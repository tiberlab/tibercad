// $Id$

#include "FieldDependentMobility.h"
#include "DriftDiffusionProperties.h"

#include "Database.h"
#include "Material.h"



TIBER_MODULE(FieldDependentMobility, field_dependent)


void
FieldDependentMobility::read_database(void)
{
  Database& db = get_database();
  db.set_section("mobility/field_dependent");

  std::vector<double> empty(2, 0);

  std::vector<double> data(empty);
  db.get("beta0", data, true);
  _beta =  get_carrier_type() == 'e' ? data[0] : data[1];

  data = empty;
  db.get("betaexp", data, true);
  _betaexp =  get_carrier_type() == 'e' ? data[0] : data[1];

  _vsat_formula = db.get("Vsat_Formula", _vsat_formula);

  if (_vsat_formula == 1)
  {
    data = empty;
    db.get("vsat0", data, true);
    _vsat0 =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("vsatexp", data, true);
    _vsat_b =  get_carrier_type() == 'e' ? data[0] : data[1];
  }
  else
  {
    data = empty;
    db.get("A_vsat", data, true);
    _vsat0 =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("B_vsat", data, true);
    _vsat_b =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("vsat_min", data, true);
    _vsat_min =  get_carrier_type() == 'e' ? data[0] : data[1];
  }

}



void
FieldDependentMobility::do_init(void)
{
  _beta = get_option("beta0", _beta);
  _betaexp = get_option("betaexp", _betaexp);
  _vsat_formula = get_option("Vsat_Formula", _vsat_formula);
  if (_vsat_formula == 1)
  {
    _vsat0 = get_option("vsat0", _vsat0);
    _vsat_b = get_option("vsatexp", _vsat_b);
  }
  else
  {
    _vsat0 = get_option("A_vsat", _vsat0);
    _vsat_b = get_option("B_vsat", _vsat_b);
    _vsat_min = get_option("vsat_min", _vsat_min);
  }

  std::string low_field_model = get_option("low_field_model", "doping_dependent");
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

  std::string force = get_option("driving_force", "grad_fermi");
  if (force == "efield")
    _force = EFIELD;
  else if (force == "grad_fermi")
    _force = GRADFERMI;
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

  const RealGradient& grad_fermi = (get_carrier_type() == 'e') ?
    get_driftdiffusionproperties().get_grad_fermi_e() :
    get_driftdiffusionproperties().get_grad_fermi_h();

  E = grad_fermi.size();

  if ((_force == EFIELD) && (E > 1e-6))
    E = grad_fermi * get_driftdiffusionproperties().get_electric_field() / E;

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
FieldDependentMobility::get_derivative_grad_fermi(RealGradient& dm)
{
  dm.zero();
  double T = get_driftdiffusionproperties().get_lattice_temperature();
  double E = 0.0;
  if (get_carrier_type() == 'e')
    E = get_driftdiffusionproperties().get_grad_fermi_e().size();
  else
    E = get_driftdiffusionproperties().get_grad_fermi_h().size();

  if ((_force == GRADFERMI) && (E > 1.0))
  {

    double vsat;
    if (_vsat_formula == 1)
      vsat = _vsat0 * std::pow(T / T0, -_vsat_b);
    else
      vsat = std::max(_vsat0 - _vsat_b * (T / T0), _vsat_min);

    double mu_lowfield = _low_field_mob->get_mobility();

    double beta = _beta * std::pow(T / T0, _betaexp);
    double tmp = 1.0 + std::pow(mu_lowfield * E / vsat, _beta);
    double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
    dmu *= std::pow(mu_lowfield * E / vsat, _beta) / (E * E);

    if (get_carrier_type() == 'e')
      dm = dmu * get_driftdiffusionproperties().get_grad_fermi_e();
    else
      dm = dmu * get_driftdiffusionproperties().get_grad_fermi_h();

  }
}


void
FieldDependentMobility::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  const FieldDependentMobility* scA =
    dynamic_cast<const FieldDependentMobility*>(comp_A);
  const FieldDependentMobility* scB =
    dynamic_cast<const FieldDependentMobility*>(comp_B);

  if (scA->_vsat_formula != scB->_vsat_formula)
    throw InitFailedException("Field dependent mobility has to use the same "
        "formula for both components of the alloy " + get_material()->get_name());

  _vsat_formula = scA->_vsat_formula;
  _force = scA->_force;

  _beta = alloy(scA->_beta, scB->_beta, xa);
  _betaexp = alloy(scA->_betaexp, scB->_betaexp, xa);
  _vsat0 = alloy(scA->_vsat0, scB->_vsat0, xa);
  _vsat_b = alloy(scA->_vsat_b, scB->_vsat_b, xa);
  _vsat_min = alloy(scA->_vsat_min, scB->_vsat_min, xa);

  destroy(_low_field_mob);
  _low_field_mob = create_submodel_copy(scA->_low_field_mob);
  _low_field_mob->set_driftdiffusionproperties(&get_driftdiffusionproperties());
  _low_field_mob->init_alloy(scA->_low_field_mob, scB->_low_field_mob, xa);
}

