#ifndef _MACROSTRAINSUPPORT_H_
#define _MACROSTRAINSUPPORT_H_
#include "MacrostrainBoundaryProperties.h"

//! Mechanical support  \f$ {\bf u} \cdot {\bf n} = 0 \f$

class MacrostrainSupport: public  MacrostrainBoundaryProperties
{
 public:

  virtual ~MacrostrainSupport() {};

  static MacrostrainSupport* create(const ModelOptions& options);

 protected:

  MacrostrainSupport(const ModelOptions& options);

    virtual void 	do_init (void);

 private:

};


inline
MacrostrainSupport::MacrostrainSupport(const ModelOptions& options)
 : MacrostrainBoundaryProperties(options)
{
}

inline MacrostrainSupport* MacrostrainSupport::create(const ModelOptions& options)
{
  return new MacrostrainSupport(options);
}


inline void MacrostrainSupport::do_init()
{
  type = "support";
}
#endif
