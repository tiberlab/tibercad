// $Id$

#include "DataImporter.h"

// TiberCAD imports
#include "TiberLinearSystem.h"
#include "Messages.h"
#include "SimulationEnvironment.h"
#include "equation_systems.h"
#include "dof_map.h"

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
  _delimiter("\\t\\s,"),
  _comment_chars({"#", "%", "!", "/"})
{
  // Move along, nothing to see here...
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
  get_option("filename", _filename);
  get_option("filetype", _filetype);
  get_option("variable_name", _variable_name);
  get_option("unit", _unit);
  get_option("variable_alias",_variable_alias);
  get_option("dataset_name", _dataset_name);
  get_option("num_dimensions",_num_dimensions);
  get_option("sizes",_sizes);
  get_option("delimiter", _delimiter);

  get_option("comment_characters", _comment_chars);
  get_option("print_data",_print_data);
  // add variables and attach assembly function
  system.add_variable("Data",FIRST);
  system.init();
  parse_options();
  _read_file();
  _print_module_info();
  increment_solve_sequence_number();
}


// This function is called first. If you need any options for the solution variables, get them here!
  void
DataImporter::do_setup_solution_variables(void)
{
  // we only have one solution variable
  Messages::info("Setting up variables...");
  declare_solution(Data,REAL,CELL,myopts.unit.c_str());
  myopts.variable_name = get_options().get_option("variable_name","");
  myopts.variable_alias = get_options().get_option("variable_alias","");
  std::stringstream MyMessage;
  MyMessage << "Adding aliases: " << myopts.variable_name << ", " << myopts.variable_alias << std::endl;
  Messages::info(MyMessage.str());
  add_alias(myopts.variable_name,Data);
  add_alias(myopts.variable_alias,Data);
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
  if (myopts.filetype == "hdf5")
    _read_hdf5();
  else if (myopts.filetype == "csv1d_col")
    _read_csv1d();
  else if (myopts.filetype == "csv1d_row")
    _read_csv1d();
  else if (myopts.filetype == "csv2d")
    _read_csv2d();
  else if (myopts.filetype == "csv3d")
    _read_csv3d();
  else if (myopts.filetype == "image")
    _read_image();
  else if (myopts.filetype == "vtk")
    _read_vtk();
  else if (myopts.filetype == "comsol")
    _read_comsol();

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

     Messages::info("Opening file "+myopts.filename);
  // Open HDF5 file and dataset info
  file_id = H5Fopen(myopts.filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  dataset_id = H5Dopen2(file_id,myopts.dataset_name.c_str(),H5P_DEFAULT);
  H5LTget_dataset_ndims(file_id,myopts.dataset_name.c_str(),&rank);
  dims = (hsize_t*)malloc(sizeof(hsize_t)*rank);
  H5LTget_dataset_info(file_id,myopts.dataset_name.c_str(),dims,NULL,NULL);
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
  H5LTread_dataset(file_id,myopts.dataset_name.c_str(),data_type,_data);
  // Close HDF5 file
  H5Fclose(file_id);
  */
}

  void
DataImporter::_read_image(void)
{
  //    Messages::info("Image Reader not yet implemented!");
  //    // TODO
  //      ostringstream os;
  //  os << "Importing data from " << _file;
  //  Messages::newline();
  //  Messages::info(os.str());
  //
  //  PLAnyPicDecoder Decoder;
  //  PLAnyBmp Bmp;
  //
  //  Decoder.MakeBmpFromFile(_file.c_str(), &Bmp);
  //
  //  int height = Bmp.GetHeight();
  //  int width = Bmp.GetWidth();
  //
  //  double scale = get_option("data_scaling", 1.0);
  //
  //  double imwidth = get_option("width", 0.0);
  //  if (imwidth > 0)
  //  {
  //    _pix_size = imwidth / width;
  //  }
  //
  //  delete _tensorgrid;
  //  Point p1(_origin);
  //  p1(0) += _pix_size * width;
  //  p1(1) += _pix_size * height;
  //
  //  double meshunits = this->get_mesh_units();
  //  _origin /= meshunits;
  //  p1 /= meshunits;
  //
  //  _tensorgrid = new TensorGrid(_origin, p1, width, height, 1);
  //
  //  //int bpp = Bmp.GetBitsPerPixel();
  //
  //  // Assuming a 24 bpp bitmap.
  //  PLPixel24 ** pLineArray;
  //  int x,y;
  //  PLPixel24  * pLine;
  //
  //  pLineArray = Bmp.GetLineArray24();
  //
  //  _data.resize(height*width);
  //
  //  // Iterate through the lines
  //  for (y = 0; y < height; y++)
  //  {
  //    // note: the image origin is the upper left corner, as usual
  //    // in informatics!
  //    pLine = pLineArray[height - y - 1];
  //
  //    // Iterate through the pixels
  //    for (x = 0; x < width; x++)
  //    {
  //      // for now we assume gray scale image
  //      double r = static_cast<double>(pLine[x].GetR());
  //      //double g = static_cast<double>(pLine[x].GetG());
  //      //double b = static_cast<double>(pLine[x].GetB());
  //      //fp_rgb_to_hsv(&r, &g, &b);
  //
  //      _data[x + y*width] = scale * r;
  //    }
  //  }
}

  void
DataImporter::_read_csv1d(void)
{
  std::ifstream in_file(myopts.filename.c_str());
  if (!in_file.is_open())
  {
    Messages::error("Could not open file! Exiting...\n");
    exit(-1);
  }
  std::string line;
  std::vector<std::string> data_str;
  // If all data is in a single line, use this
  if (myopts.filetype.find("_col") != std::string::npos)
  {
    getline(in_file,line);
    boost::split(data_str,line,boost::is_any_of(myopts.delimiter.c_str()));
  }
  // otherwise each line is a single value
  else
  {
    while (getline(in_file,line))
    {
      data_str.push_back(line);
    }
  }
  _dims = 1;
  _size_x = data_str.size();
  _size_y = 1;
  _size_z = 1;
  _data = (double*)malloc(sizeof(double)*_size_x);
#pragma omp parallel for
  for (int i=0;i<_size_x;i++)
  {
    _data[i] = atof(data_str[i].c_str());
  }
}

void
DataImporter::_read_csv2d(void)
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
  size_t data_size = 1000;
  unsigned int n_values = 0;

  string line;

  vector<string> splitted;

  while (getline(in_file, line))
  {
    // skip commented lines
    if (line.empty() || _comment_chars.count(line[0]))
      continue;

    // split up the line
    boost::split(splitted, line, boost::is_any_of(_delimiter));

    if (ctr == 0)
    {
      // this is the first data line
      n_values = splitted.size();
      data.resize(n_values);
      for (unsigned int i = 0; i < n_values; ++i)
      data[i].reserve(data_size);
    }
    else if (ctr == data_size)
    {
      for (unsigned int i = 0; i < n_values; ++i)
        data[i].reserve(data.size() + data_size);
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
  }


}

  int
DataImporter::_at(int pos_x, int pos_y)
{
  return _size_y*pos_x + pos_y;
}

  int
DataImporter::_at(int pos_x, int pos_y, int pos_z)
{
  return _size_y*_size_z*pos_x + _size_z*pos_y + pos_z;
}

  void
DataImporter::_read_csv3d(void)
{
  std::ifstream in_file(myopts.filename.c_str());
  if (!in_file.is_open())
  {
    Messages::error("Could not open file! Exiting...\n");
    exit(-1);
  }
  std::string line;
  _size_x = 0;
  _size_y = 0;
  _size_z = 0;
  std::vector<int> pos_x_list, pos_y_list, pos_z_list;
  std::vector<double> value_list;
  // Get the sizes and point lists
  while (getline(in_file,line))
  {
    std::vector<std::string> current_split;
    boost::split(current_split,line,boost::is_any_of(myopts.delimiter.c_str()));
    int pos_x = atoi(current_split[0].c_str());
    int pos_y = atoi(current_split[0].c_str());
    int pos_z = atoi(current_split[0].c_str());
    double value = atof(current_split[0].c_str());
    if (pos_x > _size_x)
      _size_x = pos_x;
    if (pos_y > _size_y)
      _size_y = pos_y;
    if (pos_z > _size_z)
      _size_z = pos_z;
    pos_x_list.push_back(pos_x);
    pos_y_list.push_back(pos_y);
    pos_z_list.push_back(pos_z);
    value_list.push_back(value);
  }

  // Initialise data
  _data = (double*)malloc(sizeof(double)*_size_x*_size_y*_size_z);
#pragma omp parallel for
  for (int i=0;i<_size_x*_size_y*_size_z;i++)
  {
    _data[i]=0;
  }
  // Fill data array
#pragma omp parallel for
  for (int i=0;i<value_list.size();i++)
  {
    int current_x = pos_x_list[i];
    int current_y = pos_y_list[i];
    int current_z = pos_z_list[i];
    double current_value = value_list[i];
    _data[_at(current_x,current_y,current_z)] = current_value;
  }
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
  // Mandatory parameters
  myopts.filename = get_options().get_option("filename","");
  myopts.filetype = get_options().get_option("filetype","");
  myopts.variable_name = get_options().get_option("variable_name","");
  // Optional parameters
  myopts.unit = get_options().get_option("unit","");
  myopts.variable_alias = get_options().get_option("variable_alias","");
  // HDF5 parameters
  // mandatory
  myopts.dataset_name = get_options().get_option("dataset","");
  // optional
  myopts.num_dimensions = get_options().get_option("num_dimensions","-1");
  myopts.sizes = get_options().get_option("sizes","-1");
  // CSV paramaters
  // optional
  myopts.delimiter = get_options().get_option("delimiter",";");
  myopts.print_data = get_options().get_option("print_data","n");
  _check_options();
}

  void
DataImporter::_check_options(void)
{
  if (myopts.filename == "")
  {
    Messages::error("No filename for reading external data provided! Exiting...\n");
    exit(-1);
  }
  if (myopts.filetype == "")
  {
    Messages::error("No filetype for external data provided! Exiting...\n");
    exit(-1);
  }
  if (myopts.filetype == "")
  {
    Messages::error("No filetype for external data specified! Exiting...\n");
    exit(-1);
  }
  bool filetype_found = false;
  for (int i=0;i<num_valid_filetypes;i++)
    if (myopts.variable_name == "")
    {
      Messages::error("Variable name must be provided! Exiting...\n");
      exit(-1);
    }
  if (myopts.filetype == "hdf5")
  {
    if (myopts.dataset_name == "")
    {
      Messages::error("HDF5 needs a dataset name! Exiting...\n");
      exit(-1);
    }
    myopts.int_num_dimensions = std::stoi(myopts.num_dimensions);
    std::vector<std::string> sizes_split;
    if (myopts.sizes != "-1")
    {
      boost::split(sizes_split,myopts.sizes,boost::is_any_of(myopts.delimiter.c_str()));
      for (int i=0;i<sizes_split.size();i++)
      {
        myopts.int_sizes.push_back(stoi(sizes_split[i]));
      }
      for (int i=sizes_split.size();i<myopts.int_num_dimensions;i++)
      {
        myopts.int_sizes.push_back(1);
      }
    }
    else
    {
      myopts.int_sizes.push_back(-1);
    }
  }
}

  void
DataImporter::_print_module_info(void)
{
  Messages::info("### Module data_from_file:");
  Messages::info("Filename: " + myopts.filename);
  Messages::info("File type: " + myopts.filetype);
  Messages::info("Variable: " + myopts.variable_name);
  Messages::info("Unit: " + myopts.unit);
  Messages::info("Variable alias: " + myopts.variable_alias);
  if (myopts.filetype == "hdf5")
  {
    Messages::info("Dataset: " + myopts.dataset_name);
    Messages::info("Number of dimensions: ", false);
    if (myopts.int_num_dimensions < 0)
    {
      Messages::info("<get from data>");
    }
    else
    {
      Messages::info(myopts.num_dimensions);
    }
    Messages::info("Sizes of data: ",false);
    if (myopts.int_sizes[0] < 0)
    {
      Messages::info("<get from data>");
    }
    else
    {
      for (int i=0;i<myopts.int_sizes.size();i++)
      {
        char *current_size;
        printf("%d ",myopts.int_sizes[i]);
      }
      printf("\n");
    }
  }
  std::size_t found = myopts.filetype.find("csv");
  if (found != std::string::npos)
  {
    Messages::info("Delimiter: " + myopts.delimiter);
  }
  if (myopts.print_data == "y")
  {
    Messages::info("Printing data...");
    for (int i=0;i<_size_x;i++)
    {
      for (int j=0;j<_size_y;j++)
      {
        for (int k=0;k<_size_z;k++)
        {
          std::stringstream MyMessage;
          MyMessage << "X: " << i << ", Y: " << j << ", Z: " << k << ", value: " <<_data[_at(i,j,k)] << std::endl;
          Messages::info(MyMessage.str());
        }
      }
    }
  }
  Messages::info("##########################");
}
