// $Id$

#ifndef _SPECULAR_H_
#define _SPECULAR_H_

#include "ThermalContact.h"
#include "Variable.h"


class Specular : public ThermalContact
{
 public:
  //!Constructor	 
  Specular(const ModelOptions& options);
  //!Destructor 
  ~Specular(){};

   
  //!Create a Reservoir object and return its pointer
  static  Specular* create(const ModelOptions& options);

 //!Return the temperature of the contact
  double get_temperature(void) const;

 protected:
  //!Initialize the model
  virtual void 	do_init (void);

  /*! \copydoc Variable::set_variable_value() */
  virtual void set_variable_value(double value, ID id = 0){};
  
  
  /*! \copydoc Variable::get_variable_value() */
  virtual double get_variable_value(ID id = 0){};


 private:

  double _temperature; 

};


inline
Specular* 
Specular::create(const ModelOptions& options)
{
  return new  Specular(options);
}


inline 
double
Specular::get_temperature( ) const
{

  return _temperature;

}

#endif
