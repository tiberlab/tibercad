// $Id$

#include "TmmBoundaryModel.h"
#include "MaterialBoundary.h"

using namespace std;

TmmBoundaryModel*
TmmBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{


  std::string type = options.get_option("type", "incident_wave");
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


