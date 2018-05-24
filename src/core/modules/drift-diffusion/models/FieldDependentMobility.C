// $Id$

#include "FieldDependentMobility.h"
#include "DriftDiffusionProperties.h"

#include "Database.h"
#include "Material.h"

#include "TiberModule.h"



void
FieldDependentMobility::read_database(void)
{
  const Database& db = get_database();
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
FieldDependentMobility::prepare_submodels(void)
{
  assert(_low_field_mob == NULL);

  ModelOptions opts;
  std::string low_field_model = get_option("low_field_model", "doping_dependent");
  opts.set_option("type", low_field_model);
  opts.set_option("particle", get_option("particle", "electron"));

  create_submodel(_low_field_mob, "mobility", opts);
  _low_field_mob->set_carrier_type(get_carrier_type());
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

  _damping = get_option("damping_parameter", _damping);

  std::string force = get_option("driving_force", "grad_fermi");
  if (force == "efield")
    _force = EFIELD;
  else if (force == "grad_fermi")
    _force = GRADFERMI;
  else if (force == "field_parameter")
    _force = FIELDPARAM;
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

  const libMesh::RealGradient& grad_fermi = (get_carrier_type() == 'e') ?
    get_driftdiffusionproperties().get_grad_fermi_e() :
    get_driftdiffusionproperties().get_grad_fermi_h();

  //E = grad_fermi.norm();

  //if ((_force == EFIELD) && (E > 1))
    //E = grad_fermi * get_driftdiffusionproperties().get_electric_field() / E;
  if (_force == GRADFERMI)
    E = grad_fermi.norm();
  else if (_force == EFIELD)
  {
    E = grad_fermi.norm();
    if (E > 1e-6)
      E = std::fabs(grad_fermi * get_driftdiffusionproperties().get_electric_field()) / E;
    else
      E = 0.0;
  }
  else if (_force == FIELDPARAM)
  {
    E = std::fabs(grad_fermi * get_driftdiffusionproperties().get_electric_field());
    E = std::sqrt(E);
  }

  double n = (get_carrier_type() == 'e') ?
      get_driftdiffusionproperties().get_electron_density() :
      get_driftdiffusionproperties().get_hole_density();

  double damp = n / (n + _damping);
  E *= damp;



  double vsat;
  if (_vsat_formula == 1)
    vsat = _vsat0 * std::pow(T / T0, -_vsat_b);
  else
    vsat = std::max(_vsat0 - _vsat_b * (T / T0), _vsat_min);

  double mu_lowfield = _low_field_mob->get_mobility();

  double beta = _beta * std::pow(T / T0, _betaexp);
  double tmp = 1.0 + std::pow(mu_lowfield * E / vsat, beta);
  double mu = mu_lowfield * std::pow(tmp, -1.0/beta);

  return mu;
}




void
FieldDependentMobility::get_derivative_grad_potential(libMesh::RealGradient& dm)
{
  dm.zero();
  double T = get_driftdiffusionproperties().get_lattice_temperature();
  double E = 0.0;
  if (get_carrier_type() == 'e')
    E = get_driftdiffusionproperties().get_grad_fermi_e().norm();
  else
    E = get_driftdiffusionproperties().get_grad_fermi_h().norm();

  double n = (get_carrier_type() == 'e') ?
      get_driftdiffusionproperties().get_electron_density() :
      get_driftdiffusionproperties().get_hole_density();

  double damp = n / (n + _damping);

  double vsat;
  if (_vsat_formula == 1)
    vsat = _vsat0 * std::pow(T / T0, -_vsat_b);
  else
    vsat = std::max(_vsat0 - _vsat_b * (T / T0), _vsat_min);

  double mu_lowfield = _low_field_mob->get_mobility();

  double beta = _beta * std::pow(T / T0, _betaexp);

  if (0) // ((_force == FIELDPARAM))
  {
    libMesh::RealGradient grad_fermi = (get_carrier_type() == 'e') ?
        get_driftdiffusionproperties().get_grad_fermi_e() :
        get_driftdiffusionproperties().get_grad_fermi_h();

    double F = grad_fermi * get_driftdiffusionproperties().get_electric_field();
    double sign = (F < 0) ? -1 : 1;
    F = damp * std::sqrt(std::fabs(F));

    if (F > 1)
    {
      double tmp = 1.0 + std::pow(mu_lowfield * F / vsat, beta);
      double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
      dmu *= std::pow(mu_lowfield * F / vsat, beta) / F;

      dm = (0.5 * dmu * sign * damp * damp / F) * grad_fermi;
    }
  }
  else if (0) //((_force == EFIELD) && (E > 1e-6))
  {
    libMesh::RealGradient grad_fermi = (get_carrier_type() == 'e') ?
        get_driftdiffusionproperties().get_grad_fermi_e() :
        get_driftdiffusionproperties().get_grad_fermi_h();

    double F = grad_fermi * get_driftdiffusionproperties().get_electric_field();
    double sign = (F < 0) ? -1 : 1;
    F = damp * std::fabs(F) / E;

    if (std::pow(mu_lowfield * F / vsat, beta) > 1e-3)
    {
      double tmp = 1.0 + std::pow(mu_lowfield * F / vsat, beta);
      double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
      dmu *= std::pow(mu_lowfield * F / vsat, beta) / F;

      dm = dmu * sign * grad_fermi / E * damp;

    }

  }
}

void
FieldDependentMobility::get_derivative_grad_fermi(libMesh::RealGradient& dm)
{
  dm.zero();
  double T = get_driftdiffusionproperties().get_lattice_temperature();
  double E = 0.0;
  if (get_carrier_type() == 'e')
    E = get_driftdiffusionproperties().get_grad_fermi_e().norm();
  else
    E = get_driftdiffusionproperties().get_grad_fermi_h().norm();

  double n = (get_carrier_type() == 'e') ?
      get_driftdiffusionproperties().get_electron_density() :
      get_driftdiffusionproperties().get_hole_density();

  double damp = n / (n + _damping);


  double vsat;
  if (_vsat_formula == 1)
    vsat = _vsat0 * std::pow(T / T0, -_vsat_b);
  else
    vsat = std::max(_vsat0 - _vsat_b * (T / T0), _vsat_min);

  double mu_lowfield = _low_field_mob->get_mobility();

  double beta = _beta * std::pow(T / T0, _betaexp);


  if ((_force == GRADFERMI) && (E > 1e-6))
  {
    E *= damp;

    double tmp = 1.0 + std::pow(mu_lowfield * E / vsat, beta);
    double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
    dmu *= std::pow(mu_lowfield * E / vsat, beta) / (E * E);

    if (get_carrier_type() == 'e')
      dm = dmu * get_driftdiffusionproperties().get_grad_fermi_e() * damp * damp;
    else
      dm = dmu * get_driftdiffusionproperties().get_grad_fermi_h() * damp * damp;

  }
  else if (0) //((_force == FIELDPARAM))
  {
    libMesh::RealGradient grad_fermi = (get_carrier_type() == 'e') ?
        get_driftdiffusionproperties().get_grad_fermi_e() :
        get_driftdiffusionproperties().get_grad_fermi_h();

    double F = grad_fermi * get_driftdiffusionproperties().get_electric_field();
    double sign = (F < 0) ? -1 : 1;
    F = damp * std::sqrt(std::fabs(F));

    if (F > 1)
    {
      double tmp = 1.0 + std::pow(mu_lowfield * F / vsat, beta);
      double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
      dmu *= std::pow(mu_lowfield * F / vsat, beta) / F;

      dm = get_driftdiffusionproperties().get_electric_field();
      dm *= 0.5 * dmu * sign * damp * damp / F;
    }
  }
  else if (0) //((_force == EFIELD) && (E > 1e-6))
  {
    libMesh::RealGradient grad_fermi = (get_carrier_type() == 'e') ?
        get_driftdiffusionproperties().get_grad_fermi_e() :
        get_driftdiffusionproperties().get_grad_fermi_h();

    double F = grad_fermi * get_driftdiffusionproperties().get_electric_field();
    double sign = (F < 0) ? -1 : 1;
    F = damp * std::fabs(F) / E;

    //if (F > 1e-6)
    {
      double tmp = 1.0 + std::pow(mu_lowfield * F / vsat, _beta);
      double dmu = -mu_lowfield * std::pow(tmp, -1.0 / beta - 1);
      dmu *= std::pow(mu_lowfield * F / vsat, beta) / F;

      dm = dmu * sign * get_driftdiffusionproperties().get_electric_field() / E * damp;

      dm -= dmu * F * grad_fermi / (E*E);
    }

  }
}
