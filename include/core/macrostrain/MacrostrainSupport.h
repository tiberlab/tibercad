#ifndef _MACROSTRAINSUPPORT_H_
#define _MACROSTRAINSUPPORT_H_
#include "MacrostrainBoundaryProperties.h"

//! Mechanical support  \f$ {\bf u} \cdot {\bf n} = 0 \f$

class MacrostrainSupport: public  MacrostrainBoundaryProperties
{
 public:

  MacrostrainSupport() {};

  virtual ~MacrostrainSupport() {};

  static MacrostrainSupport* create(void);

 protected:

    virtual void 	do_init (void);

 private:

};


inline MacrostrainSupport* MacrostrainSupport::create()
{
  return new MacrostrainSupport;
}


inline void MacrostrainSupport::do_init()
{
  type = "support";
}
#endif
