#ifndef _WZ_LATTICETHERMALCONDUCTIVITY_H_
#define _WZ_LATTICETHERMALCONDUCTIVITY_H_

#include "LatticeThermalConductivity.h"
class WzLatticeThermalConductivity: public LatticeThermalConductivity
{
 public:
  //!constructor
  WzLatticeThermalConductivity() {};

  //!destructor
  ~WzLatticeThermalConductivity() {};


  inline  static WzLatticeThermalConductivity* create();


 //! Update the lattice thermal conductivity given the Temperature
   virtual void update_tensor(); 

  /*!
    \f$ k_x = k_y = \frac{1}{a_x + b_x  T + c_x  T^2 }  \f$


    \f$ k_z = \frac{1}{a_z + b_z  T + c_z  T^2 }  \f$
    If kappa_model is constant this method doesn't update the lattice thermal conductivity
  */



 private:

  std::string _kappa_model;

  double _kappa_a_x; 
  double _kappa_b_x; 
  double _kappa_c_x;
  double _kappa_x;

  double _kappa_a_z; 
  double _kappa_b_z; 
  double _kappa_c_z;
  double _kappa_z;

 protected:
 
  virtual void read_database(void);


  virtual void do_init(void);


  inline  virtual PhysicalModelInterface*  create_new (void) const;

};

WzLatticeThermalConductivity* WzLatticeThermalConductivity::create()
{
  return (new WzLatticeThermalConductivity());
}


PhysicalModelInterface* WzLatticeThermalConductivity::create_new (void) const
{
  return (new WzLatticeThermalConductivity() ); 
}

#endif
