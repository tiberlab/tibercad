#ifndef _THERMALRESISTANCE_H_
#define _THERMALRESISTANCE_H_

#include "ThermalContact.h"



class ThermalResistance : public ThermalContact
{
 public:
  //!Constructor	 
    ThermalResistance();
  //!Destructor 
  ~ThermalResistance(){};
  //!Return the electrons_resistivity of the contact
  double get_thermal_resistance(void) const;

  //!Return the holes_resistivity of the contact
  double get_external_temperature(void) const;

  //!Set the thermal resistance
  void set_thermal_resistance(double rth);

  //!Set the external temperature
  void set_external_temperature(double text);

  //!Create a Reservoir object and return its pointer
  static  ThermalResistance* create(void);

 protected:
  //!Initialize the model
  virtual void 	do_init (void);


 private:

  double _text; 

  double _rth; 

};


inline
ThermalResistance* 
ThermalResistance::create()
{ 

  return new ThermalResistance();
}

inline 
double
ThermalResistance::get_thermal_resistance(void) const
{

  return _rth;

}

inline 
double
ThermalResistance::get_external_temperature( ) const
{

  return _text;

}


inline 
void  
ThermalResistance::set_thermal_resistance(double rth)

{

   _rth = rth;

}

inline 
void  
ThermalResistance::set_external_temperature(double text) 
{

   _text = text;

}


#endif
