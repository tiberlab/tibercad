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

  E_t_ = db.get("Etrap", E_t_);
  Talpha_e_ = db.get("Talpha_e", Talpha_e_);
  Talpha_h_ = db.get("Talpha_h", Talpha_h_);
  
  double N = get_material()->get_total_doping_density();

  // electrons
  double taumin = db.get("taumin_e", 0.0);
  double taumax = db.get("taumax_e", 1.0e-9);
  double Nref = db.get("Nref_e", 1e16);
  double g = db.get("gamma_e", 1.0);
  double denom = 1.0 + std::pow(N / Nref, g);
  tau_n_ = taumin + (taumax - taumin) / denom; 

  // holes
  taumin = db.get("taumin_h", 0.0);
  taumax = db.get("taumax_h", 1.0e-9);
  Nref = db.get("Nref_h", 1e16);
  g = db.get("gamma_h", 1.0);
  denom = 1.0 + std::pow(N / Nref, g);
  tau_p_ = taumin + (taumax - taumin) / denom; 
}



void
SRHRecombination::do_init(void)
{
  tau_n_ = get_parameter("tau_n", tau_n_);
  tau_p_ = get_parameter("tau_p", tau_p_);
  E_t_   = get_parameter("E_trap", E_t_);
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

  long double f = std::exp(E_t_ / T);

  long double tau_n = tau_n_ * std::pow(T / T0, Talpha_e_);
  long double tau_p = tau_p_ * std::pow(T / T0, Talpha_h_);

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

  long double f = std::exp(E_t_ / T);

  long double tau_n = tau_n_ * std::pow(T / T0, Talpha_e_);
  long double tau_p = tau_p_ * std::pow(T / T0, Talpha_h_);

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

  tau_n_ = alloy(scA->tau_n_, scB->tau_n_, xa);
  tau_p_ = alloy(scA->tau_p_, scB->tau_p_, xa);
  E_t_ = alloy(scA->E_t_, scB->E_t_, xa);
  Talpha_e_ = alloy(scA->Talpha_e_, scB->Talpha_e_, xa);
  Talpha_h_ = alloy(scA->Talpha_h_, scB->Talpha_h_, xa);
}

