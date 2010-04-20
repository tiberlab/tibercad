// $Id$

#ifndef _MESOTHERMALCONTACT_H_
#define _MESOTHERMALCONTACT_H_

#include "ThermalContact.h"
#include "Variable.h"


class MesoThermalContact : public ThermalContact
{
 public:
  //!Constructor	 
  MesoThermalContact(const ModelOptions& options);
  //!Destructor 
  ~MesoThermalContact(){};

  double get_temperature(void) const;

  double get_diffusivity(void) const;

  double get_reflectivity(void) const;

  double get_absorbivity(void) const;

  //!Create a Reservoir object and return its pointer
  static MesoThermalContact * create(const ModelOptions& options);

 protected:

 //!Initialize the model
  virtual void 	do_init (void);

 private:

  double _temperature; 
  double _diffusivity;
  double _reflectivity;
  double _absorbivity;


};


inline
MesoThermalContact* 
MesoThermalContact::create(const ModelOptions& options)
{
  return new  MesoThermalContact(options);
}

inline 
double
MesoThermalContact::get_temperature( ) const
{

  return _temperature;

}



inline 
double
 MesoThermalContact::get_diffusivity( ) const
{

  return _diffusivity;

}


inline 
double
MesoThermalContact::get_reflectivity( ) const
{

  return _reflectivity;

}

inline 
double
MesoThermalContact::get_absorbivity( ) const
{

  return _absorbivity;

}

#endif
