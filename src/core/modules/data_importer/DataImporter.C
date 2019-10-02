// $Id$

#include "DataImporter.h"

// TiberCAD imports
#include "TiberLinearSystem.h"
#include "Messages.h"
#include "SimulationEnvironment.h"
#include "TensorGrid.h"

#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_tetgen_interface.h"
#include "libmesh/mesh_triangle_interface.h"
// unfortunately tetgen defines REAL
#undef REAL

// STL imports
#include <stdlib.h>
#include <fstream>
#include <sstream>

// HDF5 imports
//#include <hdf5.h>
//#include <hdf5_hl.h>

// ImageReader imports
//#include "planydec.h"
//#include "plpngenc.h"
//#include "planybmp.h"
//#include "Filter/plhsvconvert.h"

// boost imports
#include <boost/algorithm/string.hpp>

// Do this for shared lib
#include "TiberModule.h"

using namespace std;
using namespace libMesh;

//! List of compatible file types
const std::string DataImporter::valid_filetypes[] = {
  "hdf5",
  "vtk_structured",
  "vtk_unstructured",
  "csv1d_col",
  "csv1d_row",
  "csv2d",
  "csv3d",
  "image",
  "comsol"
};

//! Number of compatible file types
const int DataImporter::num_valid_filetypes = 9;

DataImporter*
DataImporter::_this = NULL;

DataImporter::DataImporter(const ModelOptions& options):
  SimulationInterface(options),
  _dims(0),
  _delimiter("\t ,"),
  _comment_chars({'#', '%', '!', '/'})
{
  //is_task(true);
}

DataImporter::~DataImporter(void)
{
  // Move along, nothing to see here...
}

DataImporter*
DataImporter::create(const ModelOptions& options)
{
  // Just copying from the example here...
  return new DataImporter(options);
}

void
DataImporter::do_init(void)
{
  // create linear equation system
  create_equation_system("linear");

  // get the reference to it
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // get the module options
  _filename = get_option("filename", _filename);
  _filetype = get_option("filetype", _filetype);
  _dims = get_option("dimensions", _dims);

  get_option("variable_name", _variable_name);
  get_option("unit", _unit);
  get_option("variable_alias",_variable_alias);
  get_option("dataset_name", _dataset_name);
  get_option("sizes",_sizes);
  _delimiter = get_option("delimiter", _delimiter);

  get_option("comment_characters", _comment_chars);
  get_option("print_data",_print_data);


  // add variables and attach assembly function
  system.add_variable("Data",FIRST);
  system.init();
  parse_options();

  _read_file();

  increment_solve_sequence_number();
}


// This function is called first. If you need any options for the solution variables, get them here!
void
DataImporter::do_setup_solution_variables(void)
{
  // we only have one solution variable
  Messages::info("Setting up variables...");
  declare_solution(Data,REAL,CELL,_unit.c_str());
  _variable_name = get_options().get_option("variable_name","");
  _variable_alias = get_options().get_option("variable_alias","");
  std::stringstream MyMessage;
  MyMessage << "Adding aliases: " << _variable_name << ", " << _variable_alias << std::endl;
  Messages::info(MyMessage.str());
  add_alias(_variable_name,Data);
  add_alias(_variable_alias,Data);
}

void
DataImporter::do_solve(void)
{
}

//  Actual computation is done here!
void
DataImporter::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  int i=0;
  if (values.count(Data))
  {
    TiberLinearSystem* system = &get_equation_system<TiberLinearSystem>();
    const NumericVector<Number>& solution = system->get_solution_vector();
    const DofMap& dof_map = system->get_dof_map();

    // Get mesh dimensions
    unsigned int dim = get_mesh().mesh_dimension();
    // Set the variable ID to the one(s) you defined before
    const unsigned int varid = system->variable_number("Data");

    // This section maps 3d equidistat data to the fem grid
    FEType fe_type = system->variable_type(varid);
    UniquePtr<FEBase> fe(build_finite_element(dim, fe_type));
    vector<unsigned int> dof_indices;

    //element shape functions
    const vector<Point>& real_pts = fe->get_xyz();
    ID subdomain = elem->subdomain_id();

    fe->reinit(elem, &p);

    dof_map.dof_indices(elem, dof_indices, varid);

    const unsigned int n_dofs = dof_indices.size();
    // Loop over grid positions
    for (unsigned int n = 0; n < p.size(); n++)
    {
      double data = 0.0;
      // Get the actual position
      int pos = _tensorgrid->find_element(real_pts[n]);
      // Only do something if the point is inside the domain
      if ((pos >= 0) && (pos < _tensorgrid->num_elements()))
        values[Data][n] = _data[pos];
    }
  }
}

//void
//DataImporter::do_assemble(EquationSystems& es, const std::string& system_name)
//{
//    es;
//    system_name;
//}

void
DataImporter::_read_file(void)
{
  if (_filetype == "hdf5")
    _read_hdf5();
  else if (_filetype == "csv")
    _read_csv();
  else if (_filetype == "image")
    _read_image();
  else if (_filetype == "vtk")
    _read_vtk();
  else if (_filetype == "comsol")
    _read_comsol();
  else
    Messages::info("Unknown filetype: " + _filetype);

  // Get the bounding box for automatically calculating local resolutions
  pair<Point,Point> b_box = get_environment().get_bounding_box();
  _origin=b_box.first;
  _bound=b_box.second;

  // Generate the grid
  Messages::info("Generating grid...");
  _tensorgrid = new TensorGrid(_origin,_bound,_size_x,_size_y,_size_z);
}

void
DataImporter::_read_comsol(void)
{
}

void
DataImporter::_read_hdf5(void)
{
  /*
     hid_t file_id, dataset_id, data_type;
     int i, j ,k, data_size;
     hsize_t* dims;
     int rank;

     Messages::info("Opening file "+_filename);
  // Open HDF5 file and dataset info
  file_id = H5Fopen(_filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  dataset_id = H5Dopen2(file_id,_dataset_name.c_str(),H5P_DEFAULT);
  H5LTget_dataset_ndims(file_id,_dataset_name.c_str(),&rank);
  dims = (hsize_t*)malloc(sizeof(hsize_t)*rank);
  H5LTget_dataset_info(file_id,_dataset_name.c_str(),dims,NULL,NULL);
  // Get the data type of the dataset
  data_type = H5Dget_type(dataset_id);

  // Get dimensions and sizes automatically
  _dims = rank;
  _size_x = dims[0];
  if (_dims > 1)
  _size_y = dims[1];
  else
  _size_y = 1;
  if (_dims > 2)
  _size_z = dims[2];
  else
  _size_z = 1;
  data_size = _size_x * _size_y * _size_z;

  // Allocate and read data
  _data = (double*)malloc(sizeof(double)*data_size);
  Messages::info("Loading dataset...");
  H5LTread_dataset(file_id,_dataset_name.c_str(),data_type,_data);
  // Close HDF5 file
  H5Fclose(file_id);
  */
}

void
DataImporter::_read_image(void)
{

}


void
DataImporter::_read_csv(void)
{
  Messages::info("Reading 2d CSV file...");

  std::ifstream in_file(_filename.c_str());
  if (!in_file.is_open())
  {
    throw InitFailedException("Could not open data file " + _filename);
  }

  // this will contain the different columns
  vector<vector<double>> data;

  size_t ctr = 0;
  const size_t reserve_size = 1000;
  unsigned int n_values = 0;

  string line;
  string last_comment;

  vector<string> splitted;

  while (getline(in_file, line))
  {
    // skip commented lines
    if (line.empty() || _comment_chars.count(line[0]))
    {
      last_comment = line;
      continue;
    }

    // split up the line
    boost::split(splitted, line, boost::is_any_of(_delimiter),
                                 boost::token_compress_on);

    if (ctr == 0)
    {
      // this is the first data line
      n_values = splitted.size();
      data.resize(n_values);
      for (unsigned int i = 0; i < n_values; ++i)
      data[i].reserve(reserve_size);
    }
    else if (ctr == data.size())
    {
      for (unsigned int i = 0; i < n_values; ++i)
        data[i].reserve(data.size() + reserve_size);
    }

    if (splitted.size() != n_values)
      throw InitFailedException("Inconsistent data columns in file " + _filename);

    for (unsigned int i = 0; i < n_values; ++i)
    {
      istringstream is(splitted[i]);
      double value;
      is >> value;
      data[i].push_back(value);
    }

    ctr++;
  }
  ostringstream os;
  os << "Read " << data.size() << " datasets of size " << data[0].size();
  Messages::info(os.str());

  // get info on column content
  // last_comment could contain the info on columns
  if (!last_comment.empty())
  {

  }

  vector<double>* x = nullptr;
  vector<double>* y = nullptr;
  vector<double>* z = nullptr;
  switch (_dims)
  {
    case 3:
      z = &(data[2]);

    case 2:
      y = &(data[1]);

    case 1:
      x = &(data[0]);
  }

  // Now create mesh
  _create_mesh_from_points(x, y, z);

}


void
DataImporter::_create_mesh_from_points(const vector<double>* x,
                                       const vector<double>* y,
                                       const vector<double>* z)
{
  unsigned int dim = 3;
  (x == nullptr) && --dim;
  (y == nullptr) && --dim;
  (z == nullptr) && --dim;

  if (dim == 0)
    throw InitFailedException("Data mesh has no x,y or z coordinates associated.");

  UnstructuredMesh* mesh = new Mesh(get_solver_communicator(), dim);

  size_t n_points = 0;

  if (x != nullptr)
    n_points = x->size();
  else if (y != nullptr)
    n_points = y->size();
  else if (z != nullptr)
    n_points = z->size();

  for (size_t i = 0; i < n_points; ++i)
  {
    Point p(0, 0, 0);
    if (x != nullptr) p(0) = (*x)[i];
    if (y != nullptr) p(1) = (*y)[i];
    if (z != nullptr) p(2) = (*z)[i];
    mesh->add_point(p);
  }

  if (dim == 1)
  {

  }
  else if (dim == 2)
  {
    TriangleInterface triangleif(*mesh);
    triangleif.triangulate();
  }
  else if (dim == 3)
  {
    TetGenMeshInterface tetgenif(*mesh);
    tetgenif.triangulate_pointset();
  }

  // check for degenerate elements and assign region ids
  MeshBase::element_iterator it(mesh->elements_begin());
  while (it != mesh->elements_end())
  {
    Elem* el = *it;

    ++it;

    if (el->volume() <= 1e-12)
    {
      // eliminate all degenerate elements
      mesh.delete_elem(el);
    }
    else
    {
      Point centroid(el->centroid());
      const Elem* dev_el = MeshUtils::search_element(
          &(get_device()->get_mesh()), centroid);
      ID id = INVALID_ID;
      if (dev_el != NULL)
        id = dev_el->subdomain_id();

      el->subdomain_id() = id;

      // eliminate all elements that seem to lie outside of the structure
      if (id == INVALID_ID)
        mesh->delete_elem(el);
    }
  }

  mesh->prepare_for_use(true);
  set_mesh(mesh);
}


void
DataImporter::_read_vtk(void)
{
  Messages::info("VTK Reader not yet implemented!");
}

void
DataImporter::do_print_info(void)
{
  Messages::info("Using external source for generation rate");
}

void
DataImporter::parse_options(void)
{
}

