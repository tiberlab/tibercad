// $Id$

#include "TmmBoundaryModel.h"
#include "MaterialBoundary.h"
#include "Tmm.h"

using namespace std;

TmmBoundaryModel*
TmmBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{

  std::string type = options.get_option("type", "incidentwave");
  TmmBoundaryModel* mod =
      PhysicalModel::create<TmmBoundaryModel>("tmm_bnd_" + type,
          boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "Tmm boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}

std::string
TmmBoundaryModel::read_type(void){
  return typer;
}

void
TmmBoundaryModel::write_type(std::string str){
  typer = str;
}


void
TmmBoundaryModel::set_elements(double a0, double a1, double a2, double a3){
  _mmm00 = a0;
  _mmm01 = a1;
  _mmm10 = a2;
  _mmm11 = a3;
}



double
TmmBoundaryModel::get_element(int elm){
  switch(elm)
  {
    case 0:
      return(_mmm00);
      break;
    case 1:
      return(_mmm01);
      break;
    case 2:
      return(_mmm10);
      break;
    case 3:
      return(_mmm11);
      break;

  }
}

