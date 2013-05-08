// $Id$

#include "ImageReader.h"
#include "TiberLinearSystem.h"
#include "TensorGrid.h"
#include "Messages.h"

#include "mesh.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_interface.h"

#include "planydec.h"
#include "plpngenc.h"
#include "planybmp.h"
#include "Filter/plhsvconvert.h"

#include "TiberModule.h"

using namespace std;



ImageReader*
ImageReader::create(const ModelOptions& options)
{
  return new ImageReader(options);
}


ImageReader::ImageReader(const ModelOptions& options) :
    SimulationInterface(options),
    _origin(0),
    _tensorgrid(NULL)
{

}


ImageReader::~ImageReader(void)
{
  delete _tensorgrid;
}

void
ImageReader::do_setup_solution_variables(void)
{
  declare_solution(Data, REAL, NODES, "");
}


void
ImageReader::do_init(void)
{
  _file = get_option("image_file", "");

  get_option("origin", _origin);
  _pix_size = get_option("pixel_size", 1e-6);



  create_equation_system("linear");
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>(0);

  system.add_variable("data", libMeshEnums::FIRST);
  system.init();

  _import_picture();

  // tag it as "solved"
  increment_solve_sequence_number();

  if (plot_solution(Data))
    plot();
}



void
ImageReader::parse_options(void)
{

}


void
ImageReader::do_solve(void)
{

}

void
ImageReader::get_solution_secure(const Elem* elem,
    map<ID, vector<double> >& values, const vector<Point>& p)
{
  if (values.count(Data))
  {
    TiberLinearSystem* system = &get_equation_system<TiberLinearSystem>();
    const NumericVector<Number>& solution = system->get_solution_vector();
    const DofMap& dof_map = system->get_dof_map();

    unsigned int dim = get_mesh().mesh_dimension();

    const unsigned int varid = system->variable_number("data");

    FEType fe_type = system->variable_type(varid);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

    vector<unsigned int> dof_indices;

    //element shape functions
    const vector<vector<Real> >& phi = fe->get_phi();
    //const vector<vector<RealGradient> >& dphi = fe->get_dphi();
    const vector<Point>& real_pts = fe->get_xyz();

    ID subdomain = elem->subdomain_id();

    fe->reinit(elem, &p);

    dof_map.dof_indices(elem, dof_indices, varid);

    const unsigned int n_dofs = dof_indices.size();

    for (unsigned int n = 0; n < p.size(); n++)
    {
      double data = 0.0;
      //for (unsigned int i = 0; i < n_dofs; i++)
      //  data += phi[i][n] * solution(dof_indices[i]);

      int pos = _tensorgrid->find_element(real_pts[n]);
      if ((pos >= 0) && (pos < _tensorgrid->num_elements()))
        data = _data[pos];

      values[Data][n] = data;
    }
  }
}

void
ImageReader::_import_picture(void)
{
  if (_file.empty()) return;

  ostringstream os;
  os << "Importing data from " << _file;
  Messages::newline();
  Messages::info(os.str());

  PLAnyPicDecoder Decoder;
  PLAnyBmp Bmp;

  Decoder.MakeBmpFromFile(_file.c_str(), &Bmp);

  int height = Bmp.GetHeight();
  int width = Bmp.GetWidth();

  double imwidth = get_option("width", 0.0);
  if (imwidth > 0)
  {
    _pix_size = imwidth / width;
  }

  delete _tensorgrid;
  Point p1(_origin);
  p1(0) += _pix_size * width;
  p1(1) += _pix_size * height;

  double meshunits = this->get_mesh_units();
  _origin /= meshunits;
  p1 /= meshunits;

  _tensorgrid = new TensorGrid(_origin, p1, width, height, 1);

  //int bpp = Bmp.GetBitsPerPixel();

  // Assuming a 24 bpp bitmap.
  PLPixel24 ** pLineArray;
  int x,y;
  PLPixel24  * pLine;

  pLineArray = Bmp.GetLineArray24();

  _data.resize(height*width);

  // Iterate through the lines
  for (y = 0; y < height; y++)
  {
    // note: the image origin is the upper left corner, as usual
    // in informatics!
    pLine = pLineArray[height - y - 1];

    // Iterate through the pixels
    for (x = 0; x < width; x++)
    {
      // for now we assume gray scale image
      double r = static_cast<double>(pLine[x].GetR());
      //double g = static_cast<double>(pLine[x].GetG());
      //double b = static_cast<double>(pLine[x].GetB());
      //fp_rgb_to_hsv(&r, &g, &b);

      _data[x + y*width] = r;
    }
  }
}
