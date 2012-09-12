
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

  if (db.has_variable("optical_epsilon00")) {
    epsilon(1, 1) = db.get("optical_epsilon00", 1);
    epsilon(2, 2) = db.get("optical_epsilon11", 1);
    epsilon(3, 3) = db.get("optical_epsilon22", 1);

    epsilon_imag(1, 1) = db.get("optical_epsilon_imag00", 0);
    epsilon_imag(2, 2) = db.get("optical_epsilon_imag11", 0);
    epsilon_imag(3, 3) = db.get("optical_epsilon_imag22", 0);
  }
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

void OpticPropsModel::do_init_alloy (const PhysicalModelInterface *comp_A,
                                           const PhysicalModelInterface *comp_B, double xa) {
  const OpticPropsModel* modA = dynamic_cast<const OpticPropsModel*>(comp_A);

  const OpticPropsModel* modB = dynamic_cast<const OpticPropsModel*>(comp_B);



  alloy(epsilon,modA->epsilon, modB->epsilon, xa);

  alloy(epsilon_imag,modA->epsilon_imag, modB->epsilon_imag, xa);

}

