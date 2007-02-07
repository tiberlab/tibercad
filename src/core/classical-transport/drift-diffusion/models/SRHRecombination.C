// $Id$

#include "SRHRecombination.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "getpot.h"


void
SRHRecombination::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  _E_t = data("Etrap", 0.0);
  
  double N = mat->get_total_doping_density();

  // electrons
  double taumin = data("taumin_e", 0.0);
  double taumax = data("taumax_e", 1.0e-9);
  double Nref = data("Nref_e", 1e16);
  double g = data("gamma_e", 1.0);
  double denom = 1.0 + std::pow(N / Nref, g);
  _tau_n = taumin + (taumax - taumin) / denom; 

  // holes
  taumin = data("taumin_h", 0.0);
  taumax = data("taumax_h", 1.0e-9);
  Nref = data("Nref_h", 1e16);
  g = data("gamma_h", 1.0);
  denom = 1.0 + std::pow(N / Nref, g);
  _tau_p = taumin + (taumax - taumin) / denom; 
}



void
SRHRecombination::do_init(void)
{
  _tau_n = get_options().get_option("tau_n", _tau_n);
  _tau_p = get_options().get_option("tau_p", _tau_p);
  _E_t   = get_options().get_option("E_t", _E_t);
}



void
SRHRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();
  double kT = dd.get_lattice_temperature();

  double f = std::exp(_E_t / kT);

  double denom = _tau_p * (n + ni * f) + _tau_n * (p + ni / f);
  recomb_e = recomb_h = (n * p - ni * ni) / denom;
}



void
SRHRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();
  double kT = dd.get_lattice_temperature();

  double f = std::exp(_E_t / kT);

  double denom = _tau_p * (n + ni * f) + _tau_n * (p + ni / f);
  double SRH = (n * p - ni * ni) / denom;

  double a = (p - _tau_p * SRH) / denom;
  double b = (n - _tau_n * SRH) / denom; 

  recomb_e[0] = recomb_h[0] = a;
  recomb_e[1] = recomb_h[1] = b;
}



void
SRHRecombination::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const SRHRecombination* scA =
    dynamic_cast<const SRHRecombination*>(comp_A);
  const SRHRecombination* scB =
    dynamic_cast<const SRHRecombination*>(comp_B);

  _tau_n = alloy(scA->_tau_n, scB->_tau_n, xa);
  _tau_p = alloy(scA->_tau_p, scB->_tau_p, xa);
}

