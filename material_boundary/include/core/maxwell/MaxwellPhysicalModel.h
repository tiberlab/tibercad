// $Id$

#ifndef _MAXWELLPHYSICALMODEL_H_
#define _MAXWELLPHYSICALMODEL_H_

#include "PhysicalModel.h"
#include "OptDielectricConstant.h"

//!Model that is necessary to solve Maxwell equations 
class MaxwellPhysicalModel: public PhysicalModel
{
 public:
  MaxwellPhysicalModel();

  virtual ~MaxwellPhysicalModel();

  static  MaxwellPhysicalModel* create();

  //!return constant pointer to dielectric function model
  inline const OptDielectricConstant* get_dielectric_constant(void) const;

 protected:

  virtual PhysicalModelInterface* create_new (void) const;
 
  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);
  
  virtual void do_init();


 private:
  //! dielectric function model
  OptDielectricConstant* _epsilon_model;

  //! copy 
  MaxwellPhysicalModel (const MaxwellPhysicalModel&  t) {};

};

inline  const OptDielectricConstant* MaxwellPhysicalModel::get_dielectric_constant(void) const
{
  return _epsilon_model;
}


inline MaxwellPhysicalModel* MaxwellPhysicalModel::create()
{
  return new MaxwellPhysicalModel();
}

#endif
