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
  double get_electrons_resistivity(void) const;

  //!Return the holes_resistivity of the contact
  double get_holes_resistivity(void) const;

  //!Set the electrons_resistivity of the contact
  void set_electrons_resistivity(double rho_e);

    //!Set the holes_resistivity of the contact
  void set_holes_resistivity(double rho_h);

  //!Create a Reservoir object and return its pointer
  static  FluxContact* create(void);

 protected:
  //!Initialize the model
  virtual void 	do_init (void);


 private:

  double _rho_e; 

  double _rho_h; 

};


inline
FluxContact* 
FluxContact::create()
{ 

  return new FluxContact();
}

inline 
double
FluxContact::get_electrons_resistivity( ) const
{

  return _rho_e;

}

inline 
double
FluxContact::get_holes_resistivity( ) const
{

  return _rho_h;

}


inline 
void  
FluxContact::set_electrons_resistivity(double rho_e ) 
{

   _rho_e = rho_e;

}

inline 
void  
FluxContact::set_holes_resistivity(double rho_h ) 
{

   _rho_h = rho_h;

}


#endif
