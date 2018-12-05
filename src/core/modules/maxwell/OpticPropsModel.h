#ifndef OPTIC_PROPS_MODEL_H_
#define OPTIC_PROPS_MODEL_H_
/*
 * ExcitonLayer.h
 *
 *  Created on: Sep 26, 2011
 *      Author: paveryan
*/

#include "TypeDefs.h"
#include "PhysicalModel.h"
#include "PhysicalModel.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "OpticPropsInterface.h"
#include <tensor_value.h>

class OpticPropsModel: public OpticPropsInterface

{
  public:
    OpticPropsModel(const ModelOptions& options);

     static OpticPropsModel* create(const ModelOptions& options);

     PhysicalModel* create_new() const;

     virtual libMesh::Complex get_dielectric_constant() const {
       return epsilon(0, 0);
     }

     virtual const TensorValue<Complex> get_optical_epsilon() const {
       return epsilon;
     }

     virtual double get_permeability_constant() const {
       return mu;
     }

     virtual double get_spml() const {
       return sPML;
     }

  protected:
     virtual void get_parameter_c(const std::string& name, Complex& var);

     virtual void read_database();

     virtual void do_init(void);

     void do_init_alloy (const PhysicalModel *comp_A,
                                                const PhysicalModel *comp_B, double xa);

     TensorValue<Complex> epsilon;
     double mu;
     double sPML;

     void do_set(double eps, double eps_imag) {
       Complex c(eps, eps_imag);

       epsilon(0, 0) = c;
       epsilon(1, 1) = c;
       epsilon(2, 2) = c;
     }

/*     void rotate_to_calculation_system(const Tensor2Gen& RotMatrix) {
       epsilon  = sym(RotMatrix * (epsilon * (RotMatrix.transpose())));
       epsilon_imag  = sym(RotMatrix * (epsilon_imag * (RotMatrix.transpose())));
     }*/
};

inline OpticPropsModel* OpticPropsModel::create(const ModelOptions& options)
{
  return new OpticPropsModel(options);
}

inline
OpticPropsModel::OpticPropsModel(const ModelOptions& options) : OpticPropsInterface(options) {
  epsilon = TensorValue<Complex>(Complex(1, 0));
  mu =  1;
  sPML = -1;
}

inline
PhysicalModel*
OpticPropsModel::create_new(void) const
{
  return new OpticPropsModel(get_options());
}

#endif

