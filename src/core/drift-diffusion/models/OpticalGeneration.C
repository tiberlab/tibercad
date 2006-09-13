#include "OpticalGeneration.h"
#include "DriftDiffusionProperties.h"


#include <typeinfo>

OpticalGeneration::OpticalGeneration(void)
  : RecombinationModelInterface(),
    _G(1e-10)

{
}

void
OpticalGeneration::set_model_options(const ModelOptions& options)
{
 
}

void
OpticalGeneration::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  

  recomb_e = recomb_h = -_G;
}

void
OpticalGeneration::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
 

  recomb_e[0] = recomb_h[0] = 0;
  recomb_e[1] = recomb_h[1] = 0;
  recomb_e[2] = recomb_h[2] = 0;
}


const std::string
OpticalGeneration::get_name(void) const
{
  return "optical direct generation";
}
