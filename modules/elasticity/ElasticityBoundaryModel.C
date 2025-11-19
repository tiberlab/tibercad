// $Id$

#include "ElasticityBoundaryModel.h"
#include "tibercad/physics/MaterialBoundary.h"

using namespace std;


ElasticityBoundaryModel*
ElasticityBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{
 

  std::string type = options.get_option("type", "clamp");

  ElasticityBoundaryModel* mod = PhysicalModel::create<ElasticityBoundaryModel>("ebnd_" + type,
      boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "Elasticity boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;

}


