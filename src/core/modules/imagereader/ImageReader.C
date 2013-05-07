// $Id$

#include "ImageReader.h"
#include "TiberLinearSystem.h"
#include "TensorGrid.h"
#include "Messages.h"

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
    _origin(0)
{

}


ImageReader::~ImageReader(void)
{

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
  double _pix_size = get_option("pixel_size", 1e-6);


  create_equation_system("linear");
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>(0);

  system.add_variable("data", libMeshEnums::FIRST);

  _import_picture();
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

  //int bpp = Bmp.GetBitsPerPixel();

  // Assuming a 24 bpp bitmap.
  PLPixel24 ** pLineArray;
  int x,y;
  PLPixel24  * pLine;

  pLineArray = Bmp.GetLineArray24();

  int height = Bmp.GetHeight();
  int width = Bmp.GetWidth();

  _data.resize(height*width);

  // Iterate through the lines
  for (y = 0; y < height; y++)
  {
    pLine = pLineArray[y];

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
