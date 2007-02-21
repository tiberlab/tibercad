#ifndef _THERMALCONTACT_H_
#define _THERMALCONTACT_H_

#include "BoundaryProperties.h"


class ThermalContact: public BoundaryProperties
{
 public:

  ThermalContact() {};

  ~ThermalContact() {};

  static  ThermalContact* create(const std::string & name,  const ModelOptions &   options );

 protected:

  virtual void 	do_init (void) = 0;


 private:

}; 



#endif
