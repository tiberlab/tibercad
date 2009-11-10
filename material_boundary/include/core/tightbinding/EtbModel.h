#ifndef _ETBMODEL_H_
#define _ETBMODEL_H_

#include "PhysicalModelInterface.h"
#include "PhysicalModel.h"
#include "TightBindingModelInterface.h"
#include "elem.h"

//!Class that contains objects necessary for Atomistic Tight Binding Calculations
class ETBModel : public TightBindingModelInterface
{
public:

  //!Constructor
  ETBModel();

  //!Destructor
  ~ETBModel();

  //!Create a new object
  static ETBModel* create();


protected:

  virtual PhysicalModelInterface* create_new(void) const;

  virtual void do_init();

};



inline 
ETBModel* ETBModel::create()
{
  return new  ETBModel();
}


inline
PhysicalModelInterface* ETBModel::create_new (void) const
{
  return new ETBModel();
}




#endif
