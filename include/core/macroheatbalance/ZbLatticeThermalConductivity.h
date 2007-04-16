#ifndef _ZB_LATTICETHERMALCONDUCTIVITY_H_
#define _ZB_LATTICETHERMALCONDUCTIVITY_H_

#include "LatticeThermalConductivity.h"
class ZbLatticeThermalConductivity: public LatticeThermalConductivity
{
 public:
  //!constructor
  ZbLatticeThermalConductivity() {};

  //!destructor
  ~ZbLatticeThermalConductivity() {};

  //! Create a ZbLatticeThermalConductivity object
 inline  static ZbLatticeThermalConductivity* create();


//! Update the lattice thermal conductivity given the Temperature
  virtual void update_tensor(); 

  /*!
     \f$ k_x = k_y = k_z \frac{1}{a + b  T + c  T^2 }  \f$
    If kappa_model is constant this method doesn't update the lattice thermal conductivity
  */





 private:

  double _kappa_a; 
  double _kappa_b; 
  double _kappa_c; 
  double _kappa;

  
  std::string _kappa_model;
  


 protected:

  virtual void read_database(void);

  virtual void do_init(void);

  inline  virtual PhysicalModelInterface*  create_new (void) const;

};

inline
ZbLatticeThermalConductivity* ZbLatticeThermalConductivity::create()
{
  return (new ZbLatticeThermalConductivity());
}

inline
PhysicalModelInterface* ZbLatticeThermalConductivity::create_new (void) const
{
  return (new ZbLatticeThermalConductivity() ); 
}

#endif
