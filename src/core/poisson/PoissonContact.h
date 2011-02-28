#ifndef _POISSONCONTACT_H_
#define _POISSONCONTACT_H_

#include "BoundaryProperties.h"



//! A class that forwards the heat transport problem to boundary condition entailed
class PoissonContact: public BoundaryProperties
{
 public:

  enum Type
  {
    Dirichlet = 0,
    Neumann = 1
  };

 //!Constructor
  PoissonContact(const ModelOptions& options) : BoundaryProperties(options) {};
  //!Destructor
  ~PoissonContact() {};

  static  PoissonContact* create(const std::string& name,  const ModelOptions&   options );

  Type get_type(void) const;

 protected:

  //!Initialize the model
  virtual void 	do_init (void) = 0;

  void set_type(Type type);


 private:

  Type type;
    

}; 

inline PoissonContact::Type PoissonContact::get_type(void) const
{
  return type;
}

inline void PoissonContact::set_type(Type type_in)
{
  type = type_in;
}








#endif
