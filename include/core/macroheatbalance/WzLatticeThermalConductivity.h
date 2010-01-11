#ifndef _WZ_LATTICETHERMALCONDUCTIVITY_H_
#define _WZ_LATTICETHERMALCONDUCTIVITY_H_

#include "LatticeThermalConductivity.h"
class WzLatticeThermalConductivity: public LatticeThermalConductivity
{
 public:
  //!constructor
  WzLatticeThermalConductivity(const ModelOptions& options) : LatticeThermalConductivity(options) {};

  //!destructor
  ~WzLatticeThermalConductivity() {};


  inline  static WzLatticeThermalConductivity* create(const ModelOptions& options);


 //! Update the lattice thermal conductivity given the Temperature
  virtual void re_init(){};


 private:

  double _kappa_x;
  double _kappa_z;

 protected:

  virtual void read_database(void);

  virtual void read_database_alloy(void);

  virtual void do_init(void);


  virtual void do_init_alloy (const PhysicalModelInterface *comp_A,
      const PhysicalModelInterface *comp_B, double xa);

  inline  virtual PhysicalModelInterface*  create_new (void) const;

};

WzLatticeThermalConductivity* WzLatticeThermalConductivity::create(const ModelOptions& options)
{
  return (new WzLatticeThermalConductivity(options));
}


PhysicalModelInterface* WzLatticeThermalConductivity::create_new (void) const
{
  return (new WzLatticeThermalConductivity(get_options()) );
}

#endif
