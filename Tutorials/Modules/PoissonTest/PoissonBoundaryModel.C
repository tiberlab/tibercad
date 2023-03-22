// $Id$

#include "PoissonBoundaryModel.h"
#include "MaterialBoundary.h"

using namespace std;

PoissonBoundaryModel*
PoissonBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{
  std::string type = options.get_option("type", "dirichlet");
  PoissonBoundaryModel* mod = 
      PhysicalModel::create<PoissonBoundaryModel>("contact_" + type, boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "Poisson boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}


