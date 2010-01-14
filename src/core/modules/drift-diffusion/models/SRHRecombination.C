// $Id$

#include "SRHRecombination.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"



TIBER_MODULE(SRHRecombination, srh)



void
SRHRecombination::read_database(void)
{
  Database& db = get_database();
  db.set_section("recombination/SRH");

  _E_t = db.get("Etrap", _E_t, true);

  std::vector<double> data(2, 0);

  db.get("Talpha", data);
  _Talpha_e = data[0];
  _Talpha_h = data[1];

  data = std::vector<double>(2, 0);
  db.get("Tcoeff", data);
  _Tcoeff_e = data[0];
  _Tcoeff_h = data[1];

  db.get("taumin", data, true);
  double taumin_e = data[0];
  double taumin_h = data[1];
  db.get("taumax", data, true);
  double taumax_e = data[0];
  double taumax_h = data[1];
  db.get("Nref", data, true);
  double Nref_e = data[0];
  double Nref_h = data[1];
  db.get("gamma", data, true);
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
SRHRecombination::do_init(void)
{
  get_parameter("tau_n", _tau_n);
  get_parameter("tau_p", _tau_p);
  get_parameter("E_trap", _E_t);
}



void
SRHRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  long double n  = dd.get_electron_density();
  long double p  = dd.get_hole_density();
  long double ni = dd.get_intrinsic_density();
  double T = dd.get_lattice_temperature();

  long double f = std::exp(_E_t / T);

  long double tau_n = _tau_n * std::pow(T / T0, _Talpha_e)
    * std::exp(_Tcoeff_e * (T / T0 - 1));
  long double tau_p = _tau_p * std::pow(T / T0, _Talpha_h)
    * std::exp(_Tcoeff_h * (T / T0 - 1));

  long double denom = tau_p * (n + ni * f) + tau_n * (p + ni / f);
  long double tmp = n * p / denom;
  recomb_e = recomb_h = tmp - ni * ni / denom;
}



void
SRHRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  long double n  = dd.get_electron_density();
  long double p  = dd.get_hole_density();
  long double ni = dd.get_intrinsic_density();
  double T = dd.get_lattice_temperature();

  long double f = std::exp(_E_t / T);

  long double tau_n = _tau_n * std::pow(T / T0, _Talpha_e)
    * std::exp(_Tcoeff_e * (T / T0 - 1));
  long double tau_p = _tau_p * std::pow(T / T0, _Talpha_h)
    * std::exp(_Tcoeff_h * (T / T0 - 1));

  long double denom = tau_p * (n + ni * f) + tau_n * (p + ni / f);
  long double tmp = n * p / denom;
  long double SRH = tmp - ni * ni / denom;

  long double a = p / denom;
  a = a - tau_p * SRH / denom;
  long double b = n / denom;
  b = b - tau_n * SRH / denom;

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

