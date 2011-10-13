// $Id: ThermalBoundaryModel.C 2362 2011-02-21 01:02:31Z gromano $

#include "BoltzmannBoundaryModel.h"
#include "MaterialBoundary.h"
using namespace std;

BoltzmannBoundaryModel*
BoltzmannBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{


  std::string type = options.get_option("type", "heat_reservoir");
  BoltzmannBoundaryModel* mod =
      PhysicalModelInterface::create<BoltzmannBoundaryModel>("boltzmann_bnd_" + type,
          boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "Boltzmann boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}



