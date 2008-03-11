// $Id$

#include "FieldDependentMobility.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "getpot.h"


TIBER_MODULE(FieldDependentMobility, field_dependent)


void
FieldDependentMobility::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());
  
  std::string s("beta0_");
  s += get_carrier_type();
  _beta = data(s.c_str(), _beta);
  
  s = "betaexp_";
  s += get_carrier_type();
  _betaexp = data(s.c_str(), _betaexp);

  _vsat_formula = data("Vsat_Formula", _vsat_formula);

  if (_vsat_formula == 1)
  {
    s = "vsat0_";
    s += get_carrier_type();
    _vsat0 = data(s.c_str(), _vsat0);

    s = "vsatexp_";
    s += get_carrier_type();
    _vsat_b = data(s.c_str(), _vsat_b);
  }
  else
  {
    s = "A_vsat_";
    s += get_carrier_type();
    _vsat0 = data(s.c_str(), _vsat0);

    s = "B_vsat_";
    s += get_carrier_type();
    _vsat_b = data(s.c_str(), _vsat_b);

    s = "vsat_min";
    s += get_carrier_type();
    _vsat_min = data(s.c_str(), _vsat_min);
  }

}



void
FieldDependentMobility::do_init(void)
{
  _beta = get_parameter("beta0", _beta);
  _betaexp = get_parameter("betaexp", _betaexp);
  _vsat_formula = get_parameter("Vsat_Formula", _vsat_formula);
  if (_vsat_formula == 1)
  {
    _vsat0 = get_parameter("vsat0", _vsat0);
    _vsat_b = get_parameter("vsatexp", _vsat_b);
  }
  else
  {
    _vsat0 = get_parameter("A_vsat", _vsat0);
    _vsat_b = get_parameter("B_vsat", _vsat_b);
    _vsat_min = get_parameter("vsat_min", _vsat_min);
  }

  std::string low_field_model = get_parameter("low_field_model", "constant");
  _low_field_mob = MobilityModelInterface::create(low_field_model);
  if (_low_field_mob == NULL)
  {
    std::string msg("FieldDependentMobility: Could not ");
    msg += "create low-field mobility model '" + low_field_model + "'.";
    throw InitFailedException(msg);
  }

  _low_field_mob->set_driftdiffusionproperties(&get_driftdiffusionproperties());
  _low_field_mob->set_carrier_type(get_carrier_type());
  _low_field_mob->set_material(get_material());
  _low_field_mob->init();

  std::string force = get_parameter("driving_force", "efield");
  if (force == "efield")
    _force = GRADFERMI;
  else if (force == "grad_fermi")
    _force = EFIELD;
  else
  {
    std::string msg("FieldDependentMobility: Unknown driving force '");
    msg += force + "'.";
    throw InitFailedException(msg);
  }
    
}



double
FieldDependentMobility::get_mobility(void)
{
  double T = get_driftdiffusionproperties().get_lattice_temperature();
  double E = 0.0;
  if (_force == GRADFERMI)
    if (get_carrier_type() == 'e')
      E = get_driftdiffusionproperties().get_grad_fermi_e().size();
    else
      E = get_driftdiffusionproperties().get_grad_fermi_h().size();
  else
    E = get_driftdiffusionproperties().get_electric_field().size();

  double vsat;
  if (_vsat_formula == 1)
    vsat = _vsat0 * std::pow(T / T0, -_vsat_b);
  else
    vsat = std::max(_vsat0 - _vsat_b * (T / T0), _vsat_min);

  double mu_lowfield = _low_field_mob->get_mobility();

  double beta = _beta * std::pow(T / T0, _betaexp);
  double tmp = 1.0 + std::pow(mu_lowfield * E / vsat, _beta);
  double mu = mu_lowfield * std::pow(tmp, -1.0/beta);

  return mu;
}



void
FieldDependentMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}



void
FieldDependentMobility::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  const FieldDependentMobility* scA =
    dynamic_cast<const FieldDependentMobility*>(comp_A);
  const FieldDependentMobility* scB =
    dynamic_cast<const FieldDependentMobility*>(comp_B);

  _beta = alloy(scA->_beta, scB->_beta, xa);
  _betaexp = alloy(scA->_betaexp, scB->_betaexp, xa);
  _vsat0 = alloy(scA->_vsat0, scB->_vsat0, xa);
  _vsat_b = alloy(scA->_vsat_b, scB->_vsat_b, xa);
  _vsat_min = alloy(scA->_vsat_min, scB->_vsat_min, xa);

  _low_field_mob->build_alloy(scA->_low_field_mob, scB->_low_field_mob, xa);
}

