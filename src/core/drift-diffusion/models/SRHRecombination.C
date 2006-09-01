#include "SRHRecombination.h"
#include "DriftDiffusionProperties.h"



SRHRecombination::SRHRecombination(void)
  : RecombinationModelInterface(),
    _tau_n(1e-9),
    _tau_p(1e-9)
{
}

void
SRHRecombination::set_model_options(const ModelOptions& options)
{
  ModelOptions::const_iterator it = options.find("tau_n");
  if (it != options.end())
    _tau_n = atof((it->second).c_str());
  
  it = options.find("tau_p");
  if (it != options.end())
    _tau_p = atof((it->second).c_str());
}

void
SRHRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();

  double denom = _tau_p * (n + ni) + _tau_n * (p + ni);
  recomb_e = recomb_h = (n * p - ni * ni) / denom;
  /*
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double tn  = _tau_n;
  double tp  = _tau_p;
  double ni2 = dd.get_intrinsic_density_squared();
  double ni  = std::sqrt(ni2);
  double denom = tp * (n + ni) + tn * (p + ni);
  double G = ni2 / denom;
  double R = (n * p) / denom;
  double SRH = R - G;
  recomb_e = recomb_h = SRH;
  */
}

void
SRHRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double dn  = dd.get_electron_density_derivative();
  double p  = dd.get_hole_density();
  double dp  = dd.get_hole_density_derivative();
  double ni = dd.get_intrinsic_density();

  double denom = _tau_p * (n + ni) + _tau_n * (p + ni);
  double SRH = (n * p - ni * ni) / denom;

  double a = (p - _tau_p * SRH) * dn / denom;
  double b = (n - _tau_n * SRH) * dp / denom; 
  /*
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double dn  = dd.get_electron_density_derivative();
  double dp  = dd.get_hole_density_derivative();
  double tn  = _tau_n;
  double tp  = _tau_p;
  double ni2 = dd.get_intrinsic_density_squared();
  double ni  = std::sqrt(ni2);
  double denom = tp * (n + ni) + tn * (p + ni);
  double G = ni2 / denom;
  double R = (n * p) / denom;
  double SRH = R - G;
  double a = (p - tp * SRH) * dn / denom;
  double b = (n - tn * SRH) * dp / denom;
  */

  recomb_e[0] = recomb_h[0] = a + b;
  recomb_e[1] = recomb_h[1] = a;
  recomb_e[2] = recomb_h[2] = b;
}
