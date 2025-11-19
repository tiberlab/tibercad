// $Id$

#include "DataImporter.h"

// TiberCAD imports
#include "tibercad/solver/TiberLinearSystem.h"
#include "tibercad/io/Messages.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/geom/TensorGrid.h"
#include "tibercad/geom/MeshUtils.h"

#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh.h"
#include "libmesh/gmsh_io.h"
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
#include "tibercad/module/TiberModule.h"

using namespace std;
using namespace libMesh;



DataImporter::DataImporter(const ModelOptions& options):
  SimulationInterface(options),
  _dims(0),
  _delimiter("\t ,"),
  _comment_chars({'#', '%', '!', '/'}),
  _translate(0.0)
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

}


void
DataImporter::setup_mesh(void)
{

  parse_options();

  _read_file();


  increment_solve_sequence_number();
}


// This function is called first. If you need any options for the solution variables, get them here!
void
DataImporter::do_setup_solution_variables(void)
{

}

void
DataImporter::do_solve(void)
{
}

void
DataImporter::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{

  TiberLinearSystem* system = &get_equation_system<TiberLinearSystem>();
  const NumericVector<Number>& solution = system->get_solution_vector();
  const DofMap& dof_map = system->get_dof_map();

  // Get mesh dimensions
  unsigned int dim = get_mesh().mesh_dimension();
  // Set the variable ID to the one(s) you defined before
  vector<unsigned int> varids;
  system->get_all_variable_numbers(varids);

  // This section maps 3d equidistant data to the fem grid
  FEType fe_type = system->variable_type(varids[0]);
  std::unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));
  vector<vector<unsigned int>> dof_indices(varids.size());

  //element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<Point>& real_pts = fe->get_xyz();
  ID subdomain = elem->subdomain_id();

  fe->reinit(elem, &p);

  for (unsigned int i = 0; i < varids.size(); ++i)
  {
    dof_map.dof_indices(elem, dof_indices[i], varids[i]);
  }

  const unsigned int n_dofs = dof_indices[0].size();
  // Loop over grid positions
  for (unsigned int n = 0; n < p.size(); n++)
  {
    for (unsigned int i = 0; i < varids.size(); ++i)
    {
      if (values.count(i))
      {
        double value = 0;
        for (unsigned int k = 0; k < dof_indices[i].size(); ++k)
        {
          //cerr << solution(dof_indices[i][k]) << "  " << phi[k][n] << "\n";
          value += solution(dof_indices[i][k]) * phi[k][n];
        }

        values[i][n] = value;

      }
    }
  }
}



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
  Messages::info("Reading csv file...");

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

  Filetype format = unknown;

  vector<string> splitted;

  while (getline(in_file, line))
  {
    boost::trim(line);

    // skip empty lines
    if (line.empty())
    {
      continue;
    }

    // skip commented lines
    if (_comment_chars.count(line[0]))
    {
      // split up the line
      boost::split(splitted, line, boost::is_any_of(_delimiter),
                                   boost::token_compress_on);

      if (splitted.size() > 2)
      {
        // try to guess format
        if ((splitted[0] == "%") && (splitted[1] == "Version:") &&
            (splitted[2] == "COMSOL"))
        {
          format = comsol;
        }
        else if ((splitted[0] == "%") && (splitted[1] == "Dimension:"))
        {
          boost::trim(splitted[2]);
          istringstream is(splitted[2]);
          is >> _dims;
        }
        else if ((splitted[0] == "%") && (splitted[1] == "Length") &&
                 (splitted[2] == "unit:"))
        {
          boost::trim(splitted[3]);
          get_options().set_option("length_units", splitted[3]);
        }
      }

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

  unsigned int xcol = 0;
  unsigned int ycol = 1;
  unsigned int zcol = 2;

  vector<string> data_vars;
  get_option("data_variables", data_vars);
  vector<string> data_units(data_vars.size(), "");
  get_option("data_units", data_units);

  // get info on column content
  // last_comment could contain the info on columns
  if (!last_comment.empty())
  {
    boost::split(splitted, last_comment, boost::is_any_of(_delimiter),
                                         boost::token_compress_on);
    if (format == comsol)
    {
      int num_vars = n_values - _dims;
      int strings_per_var = (splitted.size() - (_dims + 1)) / num_vars;

      data_vars.resize(0);
      data_vars.reserve(num_vars);
      data_units = data_vars;
      for (unsigned int i = _dims + 1; i < (splitted.size() - 1); i += strings_per_var)
      {
        boost::trim(splitted[i]);
        data_vars.push_back(splitted[i]);
        string unit(splitted[i+1]);
        unit.erase(0, 1);
        unit.pop_back();
        data_units.push_back(unit);
      }
    }

  }

  if (data_units.size() != data_vars.size())
  {
    data_units.resize(0);
    data_units.resize(data_vars.size(), "");
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

  //
  //now create the system and variables
  //

  // create linear equation system
  create_equation_system("linear");

  // get the reference to it
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // add variables
  for (unsigned int i = 0; i < data_vars.size(); ++i)
  {
    system.add_variable(data_vars[i], FIRST);
    declare_solution_ext(data_vars[i], i, SolutionDescriptor::REAL,
        SolutionDescriptor::NODES, data_units[i]);
  }

  system.init();

  if (get_mesh().n_nodes() != data[0].size())
  {
    Messages::warning("Automatically generated data import mesh has"
        " different number of nodes than provided in file.");
  }

  for (size_t n = 0; n < get_mesh().n_nodes(); ++n)
  {

    if (get_mesh().node_ptr(n)->n_vars(system.number()) > 0)
    {
      for (unsigned int i = _dims; i < n_values; ++i)
      {
        unsigned int var = i - _dims;
        if (get_mesh().node_ptr(n)->n_comp(system.number(), var) > 0)
        {
          dof_id_type dofid = get_mesh().node_ptr(n)->dof_number(system.number(), var, 0);
          get_solution_vector().set(dofid, data[i][n]);
        }
      }
    }
  }

  get_solution_vector().close();
  system.update();
  //get_solution_vector().print_matlab("test.m");
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

  // as default we assume mesh units for all lengths
  double length_units = this->get_mesh_units();

  double m = 1;
  string units = get_option("length_units", "");
  if (units == "m") m = 1;
  else if (units == "cm") m = 0.01;
  else if (units == "mm") m = 0.001;
  else if (units == "um") m = 1e-6;
  else if (units == "nm") m = 1e-9;
  else
  {
    m = get_option("length_units", length_units);
  }

  m /= length_units;
  //cerr << "multiplier = " << m << endl;

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
    if (x != nullptr) p(0) = m * (*x)[i];
    if (y != nullptr) p(1) = m * (*y)[i];
    if (z != nullptr) p(2) = m * (*z)[i];
    mesh->add_point(p + _translate);
  }

  if (dim == 1)
  {

  }
  else if (dim == 2)
  {
    //TriangleInterface triangleif(*mesh);
    //triangleif.triangulate();
    MeshUtils::triangulate_point_set(*mesh);
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
      mesh->delete_elem(el);
    }
    else
    {
      Point centroid(el->vertex_average());
      const Elem* dev_el = MeshUtils::search_element(
          &(get_environment().get_device().get_mesh()), centroid);

      if (dev_el == nullptr)
      {
        for (unsigned int n = 0; n < el->n_nodes(); ++n)
        {
          const Elem* tmp_el = MeshUtils::search_element(
              &(get_environment().get_device().get_mesh()), el->point(n));
          if (tmp_el != nullptr)
          {
            dev_el = tmp_el;
            break;
          }
        }
      }

      ID id = INVALID_ID;
      if (dev_el != nullptr)
        id = dev_el->subdomain_id();

      el->subdomain_id() = id;

      // 24/11/2021 if we delete elements here, then orphaned nodes will be deleted
      // during mesh preparation
      // eliminate all elements that seem to lie outside of the structure
      //if (id == INVALID_ID)
      //  mesh->delete_elem(el);
    }
  }

  //mesh->print_info();
  //libMesh::GmshIO(*mesh).write("mesh.msh");

  if (mesh->n_elem() == 0)
    throw InitFailedException("Data import results in empty mesh: "
        "check units or geometry.");


  mesh->allow_renumbering(false);
  mesh->prepare_for_use();
  set_mesh(mesh);
  //cerr << "created mesh with " << mesh->n_elem()
  //    << " elements (local " << mesh->n_local_elem() << ")\n";
  //cerr << "mesh comm size = " << mesh->comm().size() << "\n";
  //cerr << "# part = " << mesh->n_partitions() << "\n";
  //mesh->print_info();
  //libMesh::GmshIO(*mesh).write("mesh.msh");
}


void
DataImporter::_read_vtk(void)
{
  Messages::info("VTK Reader not yet implemented!");
}

void
DataImporter::do_print_info(void)
{
}

void
DataImporter::parse_options(void)
{

  // get the module options
  _filename = get_option("filename", _filename);
  _filetype = get_option("filetype", _filetype);
  _dims = get_option("dimensions", _dims);

  _delimiter = get_option("delimiter", _delimiter);

  get_option("comment_characters", _comment_chars);

  get_option("coordinate_translation", _translate);

}

