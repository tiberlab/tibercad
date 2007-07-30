#include "MaxwellEquations.h"
#include "MaxwellPhysicalModel.h"

using namespace std;

//=======================================================================================================//
PhysicalModel*  MaxwellEquations::create_physical_model(const ModelOptions& options) const throw (ModelErrorException)
{
  MaxwellPhysicalModel* model = dynamic_cast<MaxwellPhysicalModel*> ( PhysicalModelInterface::create("maxwell", options) );
 
  if (model == NULL)
    throw ModelErrorException("MaxwellEquations: cannot create MaxwellPhysicalMode");
 
  return(model);

}

//=======================================================================================================//
