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
  DftbModel(const ModelOptions& options);

  //!Destructor
  ~DftbModel();

  //!Create a new object
  static DftbModel* create(const ModelOptions& options);


protected:

  virtual PhysicalModelInterface* create_new (void) const;

  virtual void do_init();

};



inline 
DftbModel* DftbModel::create(const ModelOptions& options)
{
  return new  DftbModel(options);
}


inline 
PhysicalModelInterface* DftbModel::create_new(void) const
{
  return new  DftbModel(get_options());
}



#endif
