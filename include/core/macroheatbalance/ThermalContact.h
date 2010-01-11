#ifndef _THERMALCONTACT_H_
#define _THERMALCONTACT_H_

#include "BoundaryProperties.h"



//! A class that forwards the heat transport problem to boundary condition entailed
class ThermalContact: public BoundaryProperties
{
 public:

  enum Type
  {
    Reservoir = 0,
    FluxContact = 1,
    ThermalSurfaceResistance = 2,
    FourierBTE = 3,
    BTEFourier = 4,
    ThermalSurfaceConductance = 5,
    Specular = 6,
    Diffusive = 7
  };

  //!Constructor
  ThermalContact(const ModelOptions& options) : BoundaryProperties(options) {};
  //!Destructor
  virtual ~ThermalContact() {};

  static  ThermalContact* create(const std::string& name,  const ModelOptions&   options );

  Type get_type(void) const;

 protected:

   
 

  //!Initialize the model
  virtual void 	do_init (void) = 0;

  void set_type(Type type);

 private:

  Type type;
    

}; 

inline ThermalContact::Type ThermalContact::get_type(void) const
{
  return type;
}

inline void ThermalContact::set_type(Type type_in)
{
  type = type_in;
}








#endif
