
/*
 * ExcitonLayer.h
 *
 *  Created on: Sep 26, 2011
 *      Author: paveryan
*/

#include "TypeDefs.h"
#include "PhysicalModelInterface.h"
#include "PhysicalModel.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "OpticPropsInterface.h"

#ifndef OPTIC_PROPS_MODEL_H_
#define OPTIC_PROPS_MODEL_H_


class OpticPropsModel: public OpticPropsInterface

{
  public:
    OpticPropsModel(const ModelOptions& options);

     static OpticPropsModel* create(const ModelOptions& options);

     PhysicalModelInterface* create_new() const;

     virtual Complex get_dielectric_constant() const {
       return Complex(epsilon(1, 1), epsilon_imag(1, 1));
     }

     virtual const Tensor2Sym get_optical_epsilon() const {
       return epsilon;
     }

     virtual const Tensor2Sym get_optical_epsilon_imag() const {
       return epsilon_imag;
     }

     virtual double get_permeability_constant() const {
       return mu;
     }

     virtual double get_spml() const {
       return sPML;
     }

  protected:
     virtual void read_database();

     virtual void do_init(void);

     void do_init_alloy (const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa);

     Tensor2Sym epsilon;
     Tensor2Sym epsilon_imag;
     double mu;
     double sPML;

     void do_set(double eps, double eps_imag) {
       epsilon(1, 1) = eps;
       epsilon(2, 2) = eps;
       epsilon(3, 3) = eps;

       epsilon_imag(1, 1) = eps_imag;
       epsilon_imag(2, 2) = eps_imag;
       epsilon_imag(3, 3) = eps_imag;
     }

     void rotate_to_calculation_system(const Tensor2Gen& RotMatrix) {
       epsilon  = sym(RotMatrix * (epsilon * (RotMatrix.transpose())));
       epsilon_imag  = sym(RotMatrix * (epsilon_imag * (RotMatrix.transpose())));
     }
};

inline OpticPropsModel* OpticPropsModel::create(const ModelOptions& options)
{
  return new OpticPropsModel(options);
}

inline
OpticPropsModel::OpticPropsModel(const ModelOptions& options) : OpticPropsInterface(options) {
  epsilon = Tensor2Sym(1);
  epsilon_imag = Tensor2Sym(0);
  mu =  1;
  sPML = -1;
}

inline
PhysicalModelInterface*
OpticPropsModel::create_new(void) const
{
  return new OpticPropsModel(get_options());
}

#endif

