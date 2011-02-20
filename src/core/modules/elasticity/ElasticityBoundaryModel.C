// $Id$

#include "ElasticityBoundaryModel.h"

using namespace std;


ElasticityBoundaryModel*
ElasticityBoundaryModel::create(const ModelOptions& options)
{
 

   std::string type = options.get_option("type", "clamp");

   std::cout<<"ebnd_" + type<<std::endl;
  ElasticityBoundaryModel* mod = dynamic_cast<ElasticityBoundaryModel*>(
      PhysicalModelInterface::create("ebnd_" + type, options));

  if (mod == NULL)
  {
    ostringstream os;
    os << "Elasticity boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;

}


