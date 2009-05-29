#ifndef _NEUMANN_H_
#define _NEUMANN_H_

#include "PoissonContact.h"
#include "Variable.h"


class Neumann : public PoissonContact
{
 public:
  //!Constructor	 
   Neumann();
  //!Destructor 
  ~Neumann(){};

  //!Return the temperature of the contact
  double get_polarization(void) const;

  //!Set the temperature of contact
  void set_polarization(double field);
   
  //!Create a Dirichlet object and return its pointer
  static  Neumann* create(void);

 protected:


  //!Initialize the model
  virtual void 	do_init (void);


 private:

  double _polarization; 

};


inline
Neumann* 
Neumann::create()
{
  return new Neumann();
}

inline 
double
Neumann::get_polarization( ) const
{

  return _polarization;

}

inline 
void  
Neumann::set_polarization(double polarization ) 
{

   _polarization = polarization;

}

#endif
