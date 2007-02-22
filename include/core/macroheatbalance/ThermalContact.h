#ifndef _THERMALCONTACT_H_
#define _THERMALCONTACT_H_

#include "BoundaryProperties.h"

//! A class that forwards the heat transport problem to boundary condition entailed
class ThermalContact: public BoundaryProperties
{
 public:
 //!Constructor
  ThermalContact() {};
  //!Destructor
  ~ThermalContact() {};

  static  ThermalContact* create(const std::string & name,  const ModelOptions &   options );

 protected:

  //!Initialize the model
  virtual void 	do_init (void) = 0;


 private:

}; 



#endif
