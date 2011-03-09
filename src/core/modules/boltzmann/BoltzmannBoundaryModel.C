// $Id: ThermalBoundaryModel.C 2362 2011-02-21 01:02:31Z gromano $

#include "BoltzmannBoundaryModel.h"

using namespace std;

BoltzmannBoundaryModel*
BoltzmannBoundaryModel::create(const ModelOptions& options)
{


  std::string type = options.get_option("type", "heat_reservoir");
  BoltzmannBoundaryModel* mod = dynamic_cast<BoltzmannBoundaryModel*>(
      PhysicalModelInterface::create("boltzmann_bnd_" + type, options));

  if (mod == NULL)
  {
    ostringstream os;
    os << "Boltzmann boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}


