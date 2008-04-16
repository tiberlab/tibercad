#ifndef _STRAINSIMULATION_H_
#define _STRAINSIMULATION_H_
 
#include "SimulationInterface.h"

class StrainSimulation:  public SimulationInterface
{
 public:

  enum Variables
  {
    EPS_XX = 0,
    EPS_YY = 1,
    EPS_ZZ = 2,
    EPS_XY = 3,
    EPS_XZ = 4,
    EPS_YZ = 5,
    P_X = 6,
    P_Y = 7,
    P_Z = 8
  };



  /*! \copydoc SimulationInterface::create_physical_model() */
  virtual PhysicalModel*

    create_physical_model(const ModelOptions& options,
			  const Material* mat) const
    throw (ModelErrorException);



 


 protected:
  
  


 private:


};

#endif
