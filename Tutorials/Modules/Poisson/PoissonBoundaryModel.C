// $Id$

#include "PoissonBoundaryModel.h"

using namespace std;

PoissonBoundaryModel*
PoissonBoundaryModel::create(const ModelOptions& options)
{
  std::string type = options.get_option("type", "dirichlet");
  PoissonBoundaryModel* mod = dynamic_cast<PoissonBoundaryModel*>(
      PhysicalModelInterface::create("mypoisson_" + type, options));

  if (mod == NULL)
  {
    ostringstream os;
    os << "Poisson boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}


