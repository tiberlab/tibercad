#ifndef _TIGHTBINDINGMODEL_H_
#define _TIGHTBINDINGMODEL_H_

#include "PhysicalModel.h"
#include "elem.h"

//!Class that contains objects necessary for Atomistic Tight Binding Calculations
class TightBindingModel : public PhysicalModel
{
public:

  //!Constructor
  TightBindingModel();

  //!Destructor
  ~TightBindingModel();

  //!Create a new object
  static TightBindingModel* create();


protected:

  virtual PhysicalModelInterface* create_new (void) const;

  virtual void copy_from(const PhysicalModelInterface *rhs){};

  virtual void do_init();

};



inline 
TightBindingModel* TightBindingModel::create()
{
  return new  TightBindingModel();
}


inline
PhysicalModelInterface* TightBindingModel::create_new (void) const
{
  return new TightBindingModel();
}




#endif
