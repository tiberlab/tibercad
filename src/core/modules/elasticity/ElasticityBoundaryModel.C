// $Id$

#include "ElasticityBoundaryModel.h"

using namespace std;


ElasticityBoundaryModel*
ElasticityBoundaryModel::create(const ModelOptions& options)
{
 

  std::string type = options.get_option("type", "clamp");

  PhysicalModelInterface* pmod = PhysicalModelInterface::create("ebnd_" + type, options);
  ElasticityBoundaryModel* mod = dynamic_cast<ElasticityBoundaryModel*>(pmod);

  if (mod == NULL)
  {
    ostringstream os;
    os << "Elasticity boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;

}


