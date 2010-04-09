// $Id$

#include "MyPoissonModel.h"


using namespace std;



MyPoissonModel*
MyPoissonModel::create(const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);
  type = "bulk_" + type;

  return dynamic_cast<MyPoissonModel*>(
      PhysicalModelInterface::create(type, options));
}
