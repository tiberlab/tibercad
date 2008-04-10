#ifndef _FLUXCONTACT_H_
#define _FLUXCONTACT_H_

#include "ThermalContact.h"



class FluxContact : public ThermalContact
{
 public:
  //!Constructor	 
   FluxContact();
  //!Destructor 
  ~FluxContact(){};
  //!Return the electrons_resistivity of the contact
  double get_heat_flux(void) const;

  //!Return the holes_resistivity of the contact
  double get_holes_resistivity(void) const;

  //!Set the electrons_resistivity of the contact
  void set_heat_flux(double rho_e);

  //!Create a Reservoir object and return its pointer
  static  FluxContact* create(void);

 protected:
  //!Initialize the model
  virtual void 	do_init (void);


 private:

  double _flux; 


};


inline
FluxContact* 
FluxContact::create()
{ 

  return new FluxContact();
}

inline 
double
FluxContact::get_heat_flux( ) const
{

  return _flux;

}

inline 
void  
FluxContact::set_heat_flux(double flux) 
{

   _flux = flux;

}


#endif
