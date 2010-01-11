// $Id$

#ifndef _BTEFOURIER_H_
#define _BTEFOURIER_H_

#include "ThermalContact.h"
#include "vector_value.h"
#include "elem.h"
#include "point.h"
#include "SimulationInterface.h"

class BTEFourier : public ThermalContact
{
 public:
  //!Constructor	 
   BTEFourier(const ModelOptions& options);
  //!Destructor 
  virtual ~BTEFourier(){};
  //!Return the electrons_resistivity of the contact

  RealGradient get_heat_flux(const Elem* elem, const Point& p);

  //!Create a Reservoir object and return its pointer
  static  BTEFourier* create(const ModelOptions& options);

 protected:
  //!Initialize the model
  virtual void 	do_init (void);


 private:

  //Flux Variables
  enum flux_variable
    {
      JQX = 0,
      JQY,  
      JQZ
    };
  std::set<ID> ID_set;
  std::map<ID,ID> var_map;

  //!Pointer to thermal simulation
  SimulationInterface* _simul;
};


inline
BTEFourier* 
BTEFourier::create(const ModelOptions& options)
{ 
  return new BTEFourier(options);
}





#endif
