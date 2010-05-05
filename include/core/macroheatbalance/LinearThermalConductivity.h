#ifndef _LINEARTHERMALCONDUCTIVITY_H_
#define _LINEARTHERMALCONDUCTIVITY_H_

#include "LatticeThermalConductivity.h"
class LinearThermalConductivity: public LatticeThermalConductivity
{
 public:
  //!constructor
  LinearThermalConductivity(const ModelOptions& options) :  LatticeThermalConductivity(options) {};

  //!destructor
  ~LinearThermalConductivity() {};

  //! Create a ZbLatticeThermalConductivity object
 inline  static LinearThermalConductivity* create(const ModelOptions& options);


//! Update the lattice thermal conductivity given the Temperature
  virtual void re_init(){};

  virtual void calculate(void);

 private:

  double _kappa;

   double kx0;
   double kz0;
   double mx;
   double mz;
   double z0;


 protected:

  virtual void read_database(void);

  virtual void do_init(void);

  virtual void do_init_alloy (const PhysicalModelInterface *comp_A,
      const PhysicalModelInterface *comp_B, double xa);

  inline  virtual PhysicalModelInterface*  create_new (void) const;

};

inline
LinearThermalConductivity* LinearThermalConductivity::create(const ModelOptions& options)
{
  return (new LinearThermalConductivity(options));
}

inline
PhysicalModelInterface* LinearThermalConductivity::create_new (void) const
{
  return (new LinearThermalConductivity(get_options()) );
}

#endif
