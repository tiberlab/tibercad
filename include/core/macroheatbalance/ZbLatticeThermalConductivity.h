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


  inline  static ZbLatticeThermalConductivity* create();

 private:

 protected:
  virtual void read_database(void);

  virtual void do_init(void);


  inline  virtual PhysicalModelInterface*  create_new (void) const;

};

ZbLatticeThermalConductivity* ZbLatticeThermalConductivity::create()
{
  return (new ZbLatticeThermalConductivity());
}


PhysicalModelInterface* ZbLatticeThermalConductivity::create_new (void) const
{
  return (new ZbLatticeThermalConductivity() ); 
}

#endif
