#ifndef _RESERVOIR_H_
#define _RESERVOIR_H_

#include "ThermalContact.h"
#include "Variable.h"


class Reservoir : public ThermalContact
{
 public:
  //!Constructor	 
  Reservoir();
  //!Destructor 
  ~Reservoir(){};
  //!Return the temperature of the contact
  double get_temperature(void) const;

  //!Set the temperature of contact
  void set_temperature(double Temp);
   
  //!Create a Reservoir object and return its pointer
  static Reservoir* create(void);

 protected:
  //!Initialize the model
  virtual void 	do_init (void);

  /*! \copydoc Variable::set_variable_value() */
  virtual void set_variable_value(double value, ID id = 0);
  
  
  /*! \copydoc Variable::get_variable_value() */
  virtual double get_variable_value(ID id = 0);


 private:

  double _temperature; 

};


inline
Reservoir* 
Reservoir::create()
{
  return new Reservoir();
}

inline 
double
Reservoir::get_temperature( ) const
{

  return _temperature;

}

inline 
void  
Reservoir::set_temperature(double Temp ) 
{

   _temperature = Temp;

}


inline
void
Reservoir::set_variable_value(double value, ID id)
{
  ignore_unused_variable(id);
  set_temperature(value);
}


inline
double
Reservoir::get_variable_value(ID id)
{
  ignore_unused_variable(id);
  return get_temperature();
}

#endif
