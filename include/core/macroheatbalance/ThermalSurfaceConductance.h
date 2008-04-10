#ifndef _THERMALSURFACECONDUCTANCE_H_
#define _THERMALSURFACECONDUCTANCE_H_

#include "ThermalContact.h"



class ThermalSurfaceConductance : public ThermalContact
{
 public:
  //!Constructor	 
    ThermalSurfaceConductance();
  //!Destructor 
  ~ThermalSurfaceConductance(){};
  //!Return the electrons_resistivity of the contact
  double get_thermal_surface_conductance(void) const;

  //!Return the holes_resistivity of the contact
  double get_temperature(void) const;

  //!Set the thermal resistance
  void set_thermal_surface_conductance(double rth);

  //!Set the external temperature
  void set_temperature(double text);

  //!Create a Reservoir object and return its pointer
  static  ThermalSurfaceConductance* create(void);

 protected:
  //!Initialize the model
  virtual void 	do_init (void);


 private:

  double _temp; 

  double _g_surf; 

};


inline
ThermalSurfaceConductance* 
ThermalSurfaceConductance::create()
{ 
  return new ThermalSurfaceConductance();
}

inline 
double
ThermalSurfaceConductance::get_thermal_surface_conductance(void) const
{

  return _g_surf;

}

inline 
double
ThermalSurfaceConductance::get_temperature( ) const
{

  return _temp;

}


inline 
void  
ThermalSurfaceConductance::set_thermal_surface_conductance(double g_surf)

{

   _g_surf = g_surf;

}

inline 
void  
ThermalSurfaceConductance::set_temperature(double temp) 
{

   _temp = temp;

}


#endif
