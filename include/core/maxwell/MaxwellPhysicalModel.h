// $Id$

#ifndef _MAXWELLPHYSICALMODEL_H_
#define _MAXWELLPHYSICALMODEL_H_

#include "PhysicalModel.h"
#include "OptDielectricConstant.h"

//!Model that is necessary to solve Maxwell equations 
class MaxwellPhysicalModel: public PhysicalModel
{
 public:

  virtual ~MaxwellPhysicalModel();

  static  MaxwellPhysicalModel* create(const ModelOptions& options);

  //!return constant pointer to dielectric function model
  inline const OptDielectricConstant* get_dielectric_constant(void) const;

 protected:

  MaxwellPhysicalModel(const ModelOptions& options);

  virtual PhysicalModelInterface* create_new (void) const;
 
  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);
  
  virtual void do_init();


 private:
  //! dielectric function model
  OptDielectricConstant* _epsilon_model;


};

inline  const OptDielectricConstant* MaxwellPhysicalModel::get_dielectric_constant(void) const
{
  return _epsilon_model;
}


inline MaxwellPhysicalModel* MaxwellPhysicalModel::create(const ModelOptions& options)
{
  return new MaxwellPhysicalModel(options);
}

#endif
