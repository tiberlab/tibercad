#ifndef _TIGHTBINDINGMODELINTERFACE_H_
#define _TIGHTBINDINGMODELINTERFACE_H_

#include "PhysicalModelInterface.h"
#include "PhysicalModel.h"

//!Class that contains objects necessary for Atomistic Tight Binding Calculations
class TightBindingModelInterface : public PhysicalModel
{
public:

  //!Constructor
  TightBindingModelInterface(const ModelOptions& options);

  //!Destructor
  ~TightBindingModelInterface(){};

  //!Create a new object
  static TightBindingModelInterface* create(const std::string& name,  const ModelOptions& options = ModelOptions());


protected:


  virtual PhysicalModelInterface* create_new (void) const = 0;

  virtual void do_init(void) = 0;

};


inline
TightBindingModelInterface::TightBindingModelInterface(const ModelOptions& options) :
  PhysicalModel(options)
{
}


inline 
TightBindingModelInterface* TightBindingModelInterface::create(const std::string& name,  const ModelOptions& options)
{
  return dynamic_cast<TightBindingModelInterface*> ( PhysicalModelInterface::create(name,options) );
}


#endif
