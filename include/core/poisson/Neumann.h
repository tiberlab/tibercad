#ifndef _NEUMANN_H_
#define _NEUMANN_H_

#include "PoissonContact.h"
#include "Variable.h"


class Neumann : public PoissonContact, public Variable
{
 public:
  //!Constructor	 
   Neumann();
  //!Destructor 
  ~Neumann(){};

  //!Return the temperature of the contact
  double get_field(void) const;

  //!Set the temperature of contact
  void set_field(double field);
   
  //!Create a Dirichlet object and return its pointer
  static  Neumann* create(void);

 protected:


  //!Initialize the model
  virtual void 	do_init (void);

  /*! \copydoc Variable::set_variable_value() */
  virtual void set_variable_value(double value, ID id = 0);
  
  
  /*! \copydoc Variable::get_variable_value() */
  virtual double get_variable_value(ID id = 0);


 private:

  double _field; 

};


inline
Neumann* 
Neumann::create()
{
  return new Neumann();
}

inline 
double
Neumann::get_field( ) const
{

  return _field;

}

inline 
void  
Neumann::set_field(double field ) 
{

   _field = field;

}


inline
void
Neumann::set_variable_value(double value, ID id)
{
  ignore_unused_variable(id);
  set_field(value);
}


inline
double
Neumann::get_variable_value(ID id)
{
  ignore_unused_variable(id);
  return get_field();
}

#endif
