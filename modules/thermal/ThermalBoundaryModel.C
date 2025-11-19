// $Id: ThermalBoundaryModel.C 2362 2011-02-21 01:02:31Z gromano $

#include "ThermalBoundaryModel.h"
#include "tibercad/physics/MaterialBoundary.h"

using namespace std;

ThermalBoundaryModel*
ThermalBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{


  std::string type = options.get_option("type", "heat_reservoir");
  ThermalBoundaryModel* mod =
      PhysicalModel::create<ThermalBoundaryModel>("thermal_bnd_" + type,
          boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "Thermal boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}


