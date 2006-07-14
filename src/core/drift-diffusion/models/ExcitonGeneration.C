#include "ExcitonGeneration.h"
#include "DriftDiffusionProperties.h"


#include <typeinfo>

ExcitonGeneration::ExcitonGeneration(void)
  : RecombinationModelInterface(),
    _C(1e-10)

{
}

void
ExcitonGeneration::set_model_options(const ModelOptions& options)
{
  ModelOptions::const_iterator it = options.find("C");
  if (it != options.end())
    _C = atof((it->second).c_str());
}

void
ExcitonGeneration::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();

  recomb_e = recomb_h = _C * (n * p - ni * ni);
}

void
ExcitonGeneration::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double dn  = dd.get_electron_density_derivative();
  double p  = dd.get_hole_density();
  double dp  = dd.get_hole_density_derivative();
  double ni = dd.get_intrinsic_density();

  double a = _C * dn * p;
  double b = _C * n * dp; 

  recomb_e[0] = recomb_h[0] = a + b;
  recomb_e[1] = recomb_h[1] = a;
  recomb_e[2] = recomb_h[2] = b;
}
