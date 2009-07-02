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
  virtual void re_init(){};



 private:

  double _kappa;


 protected:

  virtual void read_database(void);

  virtual void do_init(void);

  virtual void do_init_alloy (const PhysicalModelInterface *comp_A,
      const PhysicalModelInterface *comp_B, double xa);

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
