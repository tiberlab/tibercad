#ifndef _THERMALSURFACERESISTANCE_H_
#define _THERMALSURFACERESISTANCE_H_

#include "ThermalContact.h"



class ThermalSurfaceResistance : public ThermalContact
{
 public:
  //!Constructor	 
  ThermalSurfaceResistance();
  //!Destructor 
  ~ ThermalSurfaceResistance(){};
  //!Return the electrons_resistivity of the contact
  double get_thermal_surface_resistance(void) const;

  //!Return the holes_resistivity of the contact
  double get_temperature(void) const;

  //!Set the thermal resistance
  void set_thermal_surface_resistance(double rth);

  //!Set the external temperature
  void set_temperature(double text);

  //!Create a Reservoir object and return its pointer
  static   ThermalSurfaceResistance* create(void);

 protected:
  //!Initialize the model
  virtual void 	do_init (void);


 private:

  double _temp; 

  double _r_surf; 

};


inline
ThermalSurfaceResistance* 
ThermalSurfaceResistance::create()
{ 

  return new ThermalSurfaceResistance();
}

inline 
double
ThermalSurfaceResistance::get_thermal_surface_resistance(void) const
{

  return _r_surf;

}

inline 
double
ThermalSurfaceResistance::get_temperature( ) const
{

  return _temp;

}


inline 
void  
ThermalSurfaceResistance::set_thermal_surface_resistance(double r_surf)

{

   _r_surf = r_surf;

}

inline 
void  
ThermalSurfaceResistance::set_temperature(double temp) 
{

   _temp = temp;

}


#endif
