// $Id$

#include "SRHRecombination.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"



TIBER_MODULE(SRHRecombination, recombination, srh)


SRHRecombination::TrapAssisted::TrapAssisted(void) :
  m_trap(0.25)
{
}

double
SRHRecombination::TrapAssisted::get_gamma(double F, double T, double Et)
{
  double Et_kT = Et / T;
  double Kn_ref = 2.0/3.0 * Et_kT;
  double Kn = 4.0/3.0 * sqrt(2 * m_trap * Constants::electron_mass *
      Constants::e * Et * Et * Et) /
      (3 * Constants::hbar * F);

  double gamma = 0;

  if (Kn > Kn_ref)
  {
    double fg = sqrt(24 * Constants::electron_mass * m_trap * Constants::e * T * T * T) / Constants::hbar;
    gamma = 2 * sqrt(3*M_PI) * F / fg * exp(F*F/(fg*fg));
  }
  else
  {
    const double a1 = 0.3480242;
    const double a2 = -0.0958798;
    const double a3 = 0.7478556;
    const double p1 = 0.47047;

    double ttn = 1.0 / (1 + p1 * (0.5 * Et_kT - Kn*0.75) / sqrt(Kn * 0.375));
    gamma = Et_kT * sqrt(2 * M_PI / (3 * Kn)) *
        (a1 * ttn + a2 * ttn * ttn + a3 * ttn * ttn * ttn) * exp(Et_kT - Kn);
  }

  return gamma;
}


SRHRecombination::~SRHRecombination(void)
{
  delete _tat;
}


void
SRHRecombination::read_database(void)
{
  const Database& db = get_database();

  if (get_option("trap", false))
    db.set_section("recombination/trap");
  else
    db.set_section("recombination/SRH");

  _E_t = db.get("Etrap", _E_t);

  std::vector<double> data(2, 0);

  db.get("Talpha", data);
  _Talpha_e = data[0];
  _Talpha_h = data[1];

  data = std::vector<double>(2, 0);
  db.get("Tcoeff", data);
  _Tcoeff_e = data[0];
  _Tcoeff_h = data[1];

  data[0] = data[1] = 1e-12;
  db.get("taumin", data);
  double taumin_e = data[0];
  double taumin_h = data[1];
  data[0] = data[1] = 1e-9;
  db.get("taumax", data);
  double taumax_e = data[0];
  double taumax_h = data[1];
  data[0] = data[1] = 1e16;
  db.get("Nref", data);
  double Nref_e = data[0];
  double Nref_h = data[1];
  data[0] = data[1] = 2;
  db.get("gamma", data);
  double g_e = data[0];
  double g_h = data[1];

  double N = get_material()->get_total_doping_density();

  // electrons
  double denom = 1.0 + std::pow(N / Nref_e, g_e);
  _tau_n = taumin_e + (taumax_e - taumin_e) / denom;

  // holes
  denom = 1.0 + std::pow(N / Nref_h, g_h);
  _tau_p = taumin_h + (taumax_h - taumin_h) / denom;
}




void
SRHRecombination::read_interface_database(void)
{
  const Database& db = get_database();

  db.set_section("recombination/surface_rec");

  _E_t = db.get("Etrap", _E_t);

  std::vector<double> data(2, 1e7);
  db.get("rec_velocity", data);

  _tau_n = 1.0 / data[0];
  _tau_p = 1.0 / data[1];

}


void
SRHRecombination::do_init(void)
{
  if (get_option("trap", false))
  {
    get_parameter("sigma_n", _sigma_n);
    get_parameter("sigma_p", _sigma_p);
    _trap = true;
  }

  get_parameter("tau_n", _tau_n);
  get_parameter("tau_p", _tau_p);

  get_parameter("Et", _E_t);
  get_parameter("Nt", _density);

  std::string tmp("m");
  tmp[0] = _energy_reference;
  tmp =  get_option("reference", tmp);
  _energy_reference = tmp[0];

  ModelOptions::submodel_iterator it(get_options().submodels_begin("trap_assisted_tunneling"));
  if (get_option("trap_assisted_tunneling", false) ||
      it != get_options().submodels_end("trap_assisted_tunneling"))
  {
    _tat = new TrapAssisted();

    if (it != get_options().submodels_end("trap_assisted_tunneling"))
      _tat->m_trap = it->second.get_option("tunneling_mass", _tat->m_trap);
  }
}



void
SRHRecombination::do_init_interface(const Material* comp_A,
    const Material* comp_B)
{
  get_parameter("rec_velocity_n", _tau_n, true, new Invert());
  get_parameter("rec_velocity_p", _tau_p, true, new Invert());

  do_init();

}


double
SRHRecombination::get_trap_level(void)
{
  double ref;
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  double Ec = dd.get_conduction_band_edge();
  double Ev = dd.get_valence_band_edge();
  switch (_energy_reference)
  {
    case 'v':
      ref = Ev + _E_t;
      break;

    case 'm':
      ref = 0.5 * (Ev + Ec) + _E_t;
      break;

    default:
      ref = Ec - _E_t;
      break;
  }

  return ref;
}



void
SRHRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double n0 = dd.get_equilibrium_electron_density();
  double p0 = dd.get_equilibrium_hole_density();
  double gn = dd.get_point_data().gamma_n;
  double gp = dd.get_point_data().gamma_p;
  double T = dd.get_lattice_temperature();
  // TODO should take carrier temperatures

  double Et = get_trap_level();
  double f = std::exp((Et - dd.get_equilibrium_fermi_level()) / T);

  double tau_n = _tau_n;
  double tau_p = _tau_p;
  if (_trap)
  {
    if (_density == 0.0) _density = 1e-32;
    tau_n = 1.0 / (_density * _sigma_n * dd.get_conduction_band().get_thermal_velocity(T));
    tau_p = 1.0 / (_density * _sigma_p * dd.get_valence_band().get_thermal_velocity(T));
  }
  else
  {
    tau_n *= std::pow(T / T0, _Talpha_e) * std::exp(_Tcoeff_e * (T / T0 - 1));
    tau_p *= std::pow(T / T0, _Talpha_h) * std::exp(_Tcoeff_h * (T / T0 - 1));
  }

  if (_tat != NULL)
  {
    // TODO what is the correct value for the two dE ??
    //double dE_n = dd.get_conduction_band_edge() - dd.get_electric_potential() +
    //    dd.get_electron_electro_chemical_potential();
    //double dE_p = dd.get_electric_potential() - dd.get_valence_band_edge() -
    //    dd.get_hole_electro_chemical_potential();
    double dE_n = _E_t; //dd.get_conduction_band_edge() - Et;
    double dE_p = _E_t; //Et - dd.get_valence_band_edge();
    double gamman = _tat->get_gamma(dd.get_electric_field().size() * 100, T, dE_n);
    double gammap = _tat->get_gamma(dd.get_electric_field().size() * 100, T, dE_p);
    gamman = 1.0 / (gamman + 1);
    gammap = 1.0 / (gammap + 1);
    tau_n *= gamman;
    tau_p *= gammap;
  }

  double denom = tau_p * (n + gn * n0 * f) + tau_n * (p + gp * p0 / f);
  double tmp = n * p / denom;
  recomb_e = recomb_h = tmp - gn * gp * n0 * p0 / denom;
}



void
SRHRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double n0 = dd.get_equilibrium_electron_density();
  double p0 = dd.get_equilibrium_hole_density();
  double gn = dd.get_point_data().gamma_n;
  double gp = dd.get_point_data().gamma_p;
  double T = dd.get_lattice_temperature();

  double Et = get_trap_level();
  double f = std::exp((Et - dd.get_equilibrium_fermi_level()) / T);

  double tau_n = _tau_n;
  double tau_p = _tau_p;
  if (_trap)
  {
    if (_density == 0.0) _density = 1e-32;
    tau_n = 1.0 / (_density * _sigma_n * dd.get_conduction_band().get_thermal_velocity(T));
    tau_p = 1.0 / (_density * _sigma_p * dd.get_valence_band().get_thermal_velocity(T));
  }
  else
  {
    tau_n *= std::pow(T / T0, _Talpha_e) * std::exp(_Tcoeff_e * (T / T0 - 1));
    tau_p *= std::pow(T / T0, _Talpha_h) * std::exp(_Tcoeff_h * (T / T0 - 1));
  }

  if (_tat != NULL)
  {
    //double dE_n = dd.get_conduction_band_edge() - dd.get_electric_potential() +
    //    dd.get_electron_electro_chemical_potential();
    //double dE_p = dd.get_electric_potential() - dd.get_valence_band_edge() -
    //    dd.get_hole_electro_chemical_potential();
    double dE_n = _E_t; //dd.get_conduction_band_edge() - Et;
    double dE_p = _E_t; //Et - dd.get_valence_band_edge();
    double gamman = _tat->get_gamma(dd.get_electric_field().size() * 100, T, dE_n);
    double gammap = _tat->get_gamma(dd.get_electric_field().size() * 100, T, dE_p);
    gamman = 1.0 / (gamman + 1);
    gammap = 1.0 / (gammap + 1);
    tau_n *= gamman;
    tau_p *= gammap;
  }

  long double denom = tau_p * (n + gn * n0 * f) + tau_n * (p + gp * p0 / f);
  long double tmp = n * p / denom;
  long double SRH = tmp - gn * gp * n0 * p0 / denom;

  long double a = p / denom;
  a = a - tau_p * SRH / denom;
  long double b = n / denom;
  b = b - tau_n * SRH / denom;

  if (a < 0) a = p / denom;
  if (b < 0) b = n / denom;

  recomb_e[0] = recomb_h[0] = a;
  recomb_e[1] = recomb_h[1] = b;
}


void
SRHRecombination::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const SRHRecombination* scA =
    dynamic_cast<const SRHRecombination*>(comp_A);
  const SRHRecombination* scB =
    dynamic_cast<const SRHRecombination*>(comp_B);

  _tau_n = alloy(scA->_tau_n, scB->_tau_n, xa);
  _tau_p = alloy(scA->_tau_p, scB->_tau_p, xa);
}

