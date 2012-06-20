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
      mu =  1;
      sPML = -1;
    }

     //! A default (empty) destructor.
     virtual ~OpticPropsModel(void) {}

     static OpticPropsModel* create(const ModelOptions& options);

     PhysicalModelInterface* create_new() const
     {
       return new OpticPropsModel(get_options());
     }

     virtual Complex get_dielectric_constant() const {
       return Complex(epsilon(1, 1), epsilon_imag(1, 1));
     }

     virtual const Tensor2Sym get_optical_epsilon() const {
       return epsilon;
     }

     virtual const Tensor2Sym get_optical_epsilon_imag() const {
       return epsilon_imag;
     }
/*
     virtual const Tensor2Sym& get_optical_epsilon_imag() const {
       return epsilon_imag;
     }
*/

     virtual double get_permeability_constant() const {
       return mu;
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

       mu = db.get("optical_mu", 1.0);

       sPML = db.get("sPML", -1);

       do_set(eps, eps_imag);

       if (db.has_variable("optical_epsilon00")) {
         epsilon(1, 1) = db.get("optical_epsilon00", 1);
         epsilon(2, 2) = db.get("optical_epsilon11", 1);
         epsilon(3, 3) = db.get("optical_epsilon22", 1);

         epsilon_imag(1, 1) = db.get("optical_epsilon_imag00", 0);
         epsilon_imag(2, 2) = db.get("optical_epsilon_imag11", 0);
         epsilon_imag(3, 3) = db.get("optical_epsilon_imag22", 0);
       }
     }

     virtual void do_init(void) {
       get_parameter("optical_mu", mu);
       get_parameter("sPML", sPML);

       if (has_parameter("optical_epsilon")) {
         double eps = get_dielectric_constant().real();
         double eps_imag = get_dielectric_constant().imag();

         get_parameter("optical_epsilon", eps);
         get_parameter("optical_epsilon_imag", eps_imag);

         do_set(eps, eps_imag);
       }

       if (has_parameter("optical_epsilon00", false)) {
         get_parameter("optical_epsilon00", epsilon(1, 1));
         get_parameter("optical_epsilon11", epsilon(2, 2));
         get_parameter("optical_epsilon22", epsilon(3, 3));
         get_parameter("optical_epsilon_imag00", epsilon_imag(1, 1));
         get_parameter("optical_epsilon_imag11", epsilon_imag(2, 2));
         get_parameter("optical_epsilon_imag22", epsilon_imag(3, 3));
       }

       const Material* mat = get_material();

       const RotatedCrystal&   cr = mat->get_rotated_crystal();

       rotate_to_calculation_system(cr.RotMatrix);
     }

     void do_init_alloy (const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa) {
       const OpticPropsModel* modA = dynamic_cast<const OpticPropsModel*>(comp_A);

       const OpticPropsModel* modB = dynamic_cast<const OpticPropsModel*>(comp_B);



       alloy(epsilon,modA->epsilon, modB->epsilon, xa);

       alloy(epsilon_imag,modA->epsilon_imag, modB->epsilon_imag, xa);

/*       const Material* mat = get_material();

       const RotatedCrystal&   cr = mat->get_rotated_crystal();

       rotate_to_calculation_system(cr.RotMatrix);*/
     }

   private:
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

#endif
