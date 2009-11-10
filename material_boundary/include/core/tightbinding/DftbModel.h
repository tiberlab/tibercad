#ifndef _DFTBMODEL_H_
#define _DFTBMODEL_H_

#include "PhysicalModelInterface.h"
#include "PhysicalModel.h"
#include "TightBindingModelInterface.h"
#include "elem.h"

//!Class that contains objects necessary for Atomistic Tight Binding Calculations
class DftbModel : public TightBindingModelInterface
{
public:

  //!Constructor
  DftbModel();

  //!Destructor
  ~DftbModel();

  //!Create a new object
  static DftbModel* create(void);


protected:

  virtual PhysicalModelInterface* create_new (void) const;

  virtual void do_init();

};



inline 
DftbModel* DftbModel::create(void)
{
  return new  DftbModel();
}


inline 
PhysicalModelInterface* DftbModel::create_new(void) const
{
  return new  DftbModel();
}



#endif
