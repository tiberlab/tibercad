#ifndef _DIFFUSIVE_H_
#define _DIFFUSIVE_H_

#include "ThermalContact.h"
#include "Variable.h"


class Diffusive : public ThermalContact
{
 public:
  //!Constructor	 
  Diffusive();
  //!Destructor 
  ~Diffusive(){};
  //!Return the temperature of the contact
  double get_temperature(void) const;

  //!Set the temperature of contact
  void set_temperature(double Temp);

  //!Return the temperature of the contact
  double get_emittivity(void) const;

  //!Set the temperature of contact
  void set_emittivity(double ems);

  //!Create a Reservoir object and return its pointer
  static Diffusive* create(void);

 protected:
  //!Initialize the model
  virtual void 	do_init (void);

  /*! \copydoc Variable::set_variable_value() */
  virtual void set_variable_value(double value, ID id = 0);
  
  
  /*! \copydoc Variable::get_variable_value() */
  virtual double get_variable_value(ID id = 0);


 private:

  double _temperature; 
  double _emittivity;
};


inline
Diffusive* 
Diffusive::create()
{
  return new Diffusive();
}

inline 
double
Diffusive::get_temperature( ) const
{

  return _temperature;

}

inline 
void  
Diffusive::set_temperature(double Temp ) 
{

   _temperature = Temp;

}

inline 
double
Diffusive::get_emittivity( ) const
{

  return _emittivity;

}

inline 
void  
Diffusive::set_emittivity(double ems) 
{

   _emittivity = ems;

}

inline
void
Diffusive::set_variable_value(double value, ID id)
{
  ignore_unused_variable(id);
  set_temperature(value);
}


inline
double
Diffusive::get_variable_value(ID id)
{
  ignore_unused_variable(id);
  return get_temperature();
}

#endif
