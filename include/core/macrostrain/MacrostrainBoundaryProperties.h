#ifndef _MACROSTRAINBOUNDARYPROPERTIES_H_
#define _MACROSTRAINBOUNDARYPROPERTIES_H_

#include "BoundaryProperties.h"
#include "Material.h"

class MacrostrainBoundaryProperties: public BoundaryProperties
{
 public:

   MacrostrainBoundaryProperties(const ModelOptions& options);

   static  MacrostrainBoundaryProperties* create(const std::string & name,  const ModelOptions &   options );

   std::string get_type(void) const {return type; };

 protected:

  virtual void 	do_init (void) = 0;

  std::string type; 

 private:

  
 

}; 


#endif
