// $Id$

#include "FieldAssistedMobility.h"
#include "DriftDiffusionProperties.h"

#include "Database.h"
#include "Material.h"



TIBER_MODULE(FieldAssistedMobility, field_assisted)


void
FieldAssistedMobility::read_database(void)
{
  Database& db = get_database();
  db.set_section("mobility/field_assisted");

  std::vector<double> mu0(2, _mu0);

  db.get("mu0", mu0);
  _mu0 =  get_carrier_type() == 'e' ? mu0[0] : mu0[1];

  std::vector<double> E0(2, _E0);
  db.get("E0", E0);
  _E0 =  get_carrier_type() == 'e' ? E0[0] : E0[1];
}



void
FieldAssistedMobility::do_init(void)
{
  get_parameter("mu0", _mu0);
  get_parameter("E0", _E0);

//  std::string force = get_option("driving_force", "efield");
//  if (force == "efield")
//    _force = EFIELD;
//  else if (force == "grad_fermi")
//    _force = GRADFERMI;
//  else
//  {
//    std::string msg("FieldAssistedMobility: Unknown driving force '");
//    msg += force + "'.";
//    throw InitFailedException(msg);
//  }

}



double
FieldAssistedMobility::get_mobility(void)
{
  //double T = get_driftdiffusionproperties().get_lattice_temperature();
  double E = 0.0;

  /*
  const RealGradient& grad_fermi = (get_carrier_type() == 'e') ?
    get_driftdiffusionproperties().get_grad_fermi_e() :
    get_driftdiffusionproperties().get_grad_fermi_h();

  E = grad_fermi.size();

  if ((_force == EFIELD) && (E > 1e-6))
    E = grad_fermi * get_driftdiffusionproperties().get_electric_field() / E;
  */

  E = get_driftdiffusionproperties().get_electric_field().size();

  double arg = std::sqrt(E / _E0);

  double mu = _mu0 * std::exp(arg);

  return mu;
}



void
FieldAssistedMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}


void
FieldAssistedMobility::get_derivative_grad_fermi(RealGradient& dm)
{
  dm.zero();
  //double T = get_driftdiffusionproperties().get_lattice_temperature();
  double E = 0.0;
  //if (get_carrier_type() == 'e')
  //  E = get_driftdiffusionproperties().get_grad_fermi_e().size();
  //else
  //  E = get_driftdiffusionproperties().get_grad_fermi_h().size();
/*
  if ((_force == EFIELD) && (E > 1.0))
  {

    double mu = get_mobility();
    E = get_driftdiffusionproperties().get_electric_field().size();
    double tmp = 0.5 / (std::sqrt(E * _E0) * E);


    dmu *= mu * tmp * get_driftdiffusionproperties().get_electric_field();

    if (get_carrier_type() == 'e')
      dm = dmu * get_driftdiffusionproperties().get_grad_fermi_e();
    else
      dm = dmu * get_driftdiffusionproperties().get_grad_fermi_h();

  }
*/
}


void
FieldAssistedMobility::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  const FieldAssistedMobility* scA =
    dynamic_cast<const FieldAssistedMobility*>(comp_A);
  const FieldAssistedMobility* scB =
    dynamic_cast<const FieldAssistedMobility*>(comp_B);

  _force = scA->_force;

  _mu0 = alloy(scA->_mu0, scB->_mu0, xa);
  _E0 = alloy(scA->_E0, scB->_E0, xa);
}

