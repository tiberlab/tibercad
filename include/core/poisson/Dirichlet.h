#ifndef _DIRICHLET_H_
#define _DIRICHLET_H_

#include "PoissonContact.h"
#include "Variable.h"


class Dirichlet : public PoissonContact
{
 public:
  //!Constructor	 
  Dirichlet();
  //!Destructor 
  ~Dirichlet(){};

  //!Return the temperature of the contact
  double get_potential(void) const;

  //!Set the temperature of contact
  void set_potential(double potential);
   
  //!Create a Dirichlet object and return its pointer
  static Dirichlet* create(void);

 protected:


  //!Initialize the model
  virtual void 	do_init (void);


 private:

  double _potential; 

};


inline
Dirichlet* 
Dirichlet::create()
{
  return new Dirichlet();
}

inline 
double
Dirichlet::get_potential( ) const
{

  return _potential;

}

inline 
void  
Dirichlet::set_potential(double potential ) 
{

   _potential = potential;

}

#endif
