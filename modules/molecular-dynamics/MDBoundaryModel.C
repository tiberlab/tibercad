// $Id: PoissonBoundaryModel.C 1902 2010-04-09 14:45:44Z maufder $

#include "MDBoundaryModel.h"

using namespace std;

MDBoundaryModel*
MDBoundaryModel::create(const MaterialBoundary* boundary,const ModelOptions& options)
{

  std::string type = options.get_option("type", "heat_reservoir");
  MDBoundaryModel* mod =
      PhysicalModel::create<MDBoundaryModel>("md_bnd_" + type,
          boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "MD boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}

