// $Id: ThermalBoundaryModel.C 2362 2011-02-21 01:02:31Z gromano $

#include "ThermalBoundaryModel.h"

using namespace std;

ThermalBoundaryModel*
ThermalBoundaryModel::create(const ModelOptions& options)
{


  std::string type = options.get_option("type", "heat_reservoir");
  ThermalBoundaryModel* mod = dynamic_cast<ThermalBoundaryModel*>(
      PhysicalModelInterface::create("thermal_bnd_" + type, options));

  if (mod == NULL)
  {
    ostringstream os;
    os << "Thermal boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}


