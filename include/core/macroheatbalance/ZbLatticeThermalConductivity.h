#ifndef _ZB_LATTICETHERMALCONDUCTIVITY_H_
#define _ZB_LATTICETHERMALCONDUCTIVITY_H_

#include "LatticeThermalConductivity.h"
class ZbLatticeThermalConductivity: public LatticeThermalConductivity
{
 public:
  //!constructor
  ZbLatticeThermalConductivity(const ModelOptions& options) :  LatticeThermalConductivity(options) {};

  //!destructor
  ~ZbLatticeThermalConductivity() {};

  //! Create a ZbLatticeThermalConductivity object
 inline  static ZbLatticeThermalConductivity* create(const ModelOptions& options);


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
ZbLatticeThermalConductivity* ZbLatticeThermalConductivity::create(const ModelOptions& options)
{
  return (new ZbLatticeThermalConductivity(options));
}

inline
PhysicalModelInterface* ZbLatticeThermalConductivity::create_new (void) const
{
  return (new ZbLatticeThermalConductivity(get_options()) );
}

#endif
