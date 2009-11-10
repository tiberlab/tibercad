// $Id$


#ifndef _MACROSTRAINMODELINTERFACE_H_
#define _MACROSTRAINMODELINTERFACE_H_

#include "PhysicalModelInterface.h"
#include "PhysicalModel.h"



class MacrostrainModelInterface : public PhysicalModel
{
 public:

  
  static MacrostrainModelInterface* create (const std::string& name,  const ModelOptions& options = ModelOptions());


 protected:

  MacrostrainModelInterface(void) { };


  virtual void do_init(void) = 0;


  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) = 0 ;


  virtual PhysicalModelInterface* create_new(void) const = 0;

};



inline MacrostrainModelInterface* MacrostrainModelInterface::create (const std::string& name,  const ModelOptions& options)
{



  return dynamic_cast<MacrostrainModelInterface*> (PhysicalModelInterface::create( name, options)  ) ;


}

#endif
