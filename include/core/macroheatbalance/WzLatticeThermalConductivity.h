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

 private:

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
