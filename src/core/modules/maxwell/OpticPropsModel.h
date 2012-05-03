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

#ifndef OPTIC_PROPS_MODEL1_H_
#define OPTIC_PROPS_MODEL1_H_

class OpticPropsModel: public PhysicalModel

{
  public:
    //OpticPropertiesModel();
     //! The empty constructor.
    OpticPropsModel(const ModelOptions& options) : PhysicalModel(options) {
      epsilon = Tensor2Sym(1);
      epsilon_imag = Tensor2Sym(0);
      mu =  Tensor2Sym(1);
      sPML = -1;
    }

     //! A default (empty) destructor.
     virtual ~OpticPropsModel(void) {}

     static OpticPropsModel* create(const ModelOptions& options);

     PhysicalModelInterface* create_new() const
     {
       return new OpticPropsModel(get_options());
     }

     //TODO
     virtual Complex get_dielectric_constant() const {
       return Complex(epsilon(1, 1), epsilon_imag(1, 1));
     }

     //TODO
     virtual double get_permeability_constant() const {
       return mu(1, 1);
     }

     virtual double get_spml() const {
       return sPML;
     }

  protected:
     virtual void read_database() {
       const Database& db = get_database();
       db.set_section("permittivity");

       double eps = db.get("optical_epsilon", 1.0);

       double eps_imag = db.get("optical_epsilon_imag", 0.0);

       double dmu = db.get("optical_mu", 1.0);

       sPML = db.get("sPML", -1);

       do_set_simple(eps, eps_imag, dmu);
     }

     virtual void do_init(void) {
       double eps = get_dielectric_constant().real();
       double eps_imag = get_dielectric_constant().imag();
       double dmu = get_permeability_constant();

       get_parameter("optical_epsilon", eps);
       get_parameter("optical_epsilon_imag", eps_imag);
       get_parameter("optical_mu", dmu);
       get_parameter("sPML", sPML);

       do_set_simple(eps, eps_imag, dmu);
     }


   private:
     Tensor2Sym epsilon;
     Tensor2Sym epsilon_imag;
     Tensor2Sym mu;
     double sPML;

     void do_set_simple(double eps, double eps_imag, double dmu) {
       epsilon(1, 1) = eps;
       epsilon(2, 2) = eps;
       epsilon(3, 3) = eps;

       epsilon_imag(1, 1) = eps_imag;
       epsilon_imag(2, 2) = eps_imag;
       epsilon_imag(3, 3) = eps_imag;

       mu(1, 1) = dmu;
       mu(2, 2) = dmu;
       mu(3, 3) = dmu;
     }
};

inline OpticPropsModel* OpticPropsModel::create(const ModelOptions& options)
{
  return new OpticPropsModel(options);
}
/*inline
void
OptDielectricConstant::rotate_to_calculation_system(const Tensor2Gen& RotMatrix)
{

  // generates dielectric  matrix in calculation system
  _dielectric_constant_real  = sym(RotMatrix * (_dielectric_constant_real * (RotMatrix.transpose())));
  _dielectric_constant_imag  = sym(RotMatrix * (_dielectric_constant_imag * (RotMatrix.transpose())));


}*/

#endif
