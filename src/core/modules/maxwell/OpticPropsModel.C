
#include "OpticPropsModel.h"

#include "SimulationInterface.h"
#include "Material.h"

void OpticPropsModel::read_database() {
  const Database& db = get_database();
  db.set_section("permittivity");

  double eps = db.get("optical_epsilon", 1.0);

  double eps_imag = db.get("optical_epsilon_imag", 0.0);

  mu = db.get("optical_mu", 1.0);

  sPML = db.get("sPML", -1);

  do_set(eps, eps_imag);


  if (db.has_variable("optical_epsilonxx")) {
    epsilon(0, 0) = Complex(db.get("optical_epsilon_xx", 1), db.get("optical_epsilon_xx_imag", 0));
    epsilon(1, 1) = Complex(db.get("optical_epsilon_yy", 1), db.get("optical_epsilon_yy_imag", 0));
    epsilon(2, 2) = Complex(db.get("optical_epsilon_zz", 1), db.get("optical_epsilon_zz_imag", 0));;

    epsilon(0, 1) = Complex(db.get("optical_epsilon_xy", 1), db.get("optical_epsilon_xy_imag", 0));
    epsilon(1, 2) = Complex(db.get("optical_epsilon_yz", 1), db.get("optical_epsilon_yz_imag", 0));
    epsilon(2, 0) = Complex(db.get("optical_epsilon_zx", 1), db.get("optical_epsilon_zx_imag", 0));;
  }
}

void OpticPropsModel::get_parameter_c(const std::string& name, Complex& var) {
  double a = var.real(), b = var.imag();

  get_parameter(name, a);
  get_parameter(name + "_imag", b);

  var = Complex(a, b);
}

void OpticPropsModel::do_init(void) {
  get_parameter("optical_mu", mu);
  get_parameter("sPML", sPML);

  if (has_parameter("optical_epsilon")) {
    double eps = get_dielectric_constant().real();
    double eps_imag = get_dielectric_constant().imag();

    get_parameter("optical_epsilon", eps);
    get_parameter("optical_epsilon_imag", eps_imag);

    do_set(eps, eps_imag);
  }

  if (has_parameter("optical_epsilonxx", false)) {
    get_parameter_c("optical_epsilon_xx", epsilon(0, 0));
    get_parameter_c("optical_epsilon_yy", epsilon(1, 1));
    get_parameter_c("optical_epsilon_zz", epsilon(2, 2));

    get_parameter_c("optical_epsilon_xy", epsilon(0, 1));
    get_parameter_c("optical_epsilon_yz", epsilon(1, 2));
    get_parameter_c("optical_epsilon_zx", epsilon(2, 0));
  }

  epsilon(1, 0) = epsilon(0, 1);
  epsilon(2, 1) = epsilon(1, 2);
  epsilon(0, 2) = epsilon(2, 0);

  const Material* mat = get_material();

  const RotatedCrystal&   cr = mat->get_rotated_crystal();

  //rotate_to_calculation_system(cr.RotMatrix); TODO
}

void OpticPropsModel::do_init_alloy (const PhysicalModelInterface *comp_A,
                                           const PhysicalModelInterface *comp_B, double xa) {
/*
  const OpticPropsModel* modA = dynamic_cast<const OpticPropsModel*>(comp_A);

  const OpticPropsModel* modB = dynamic_cast<const OpticPropsModel*>(comp_B);

  epsilon =  (1 - xa) * modB->epsilon + xa * modA->epsilon;
  mu = alloy(modA->mu, modB->mu, xa);*/
}

