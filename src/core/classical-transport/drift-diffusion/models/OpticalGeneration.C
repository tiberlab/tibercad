#include "OpticalGeneration.h"
#include "DriftDiffusionProperties.h"
#include "Material.h"


#include <typeinfo>



//TIBER_MODULE(OpticalGeneration, optical)




void
OpticalGeneration::do_init(void)
{
  _G = get_options().get_option("G", 1e-10);
  _G = get_material()->get_options().get_option("G", _G);
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
}


void
OpticalGeneration::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const OpticalGeneration* scA =
    dynamic_cast<const OpticalGeneration*>(comp_A);
  const OpticalGeneration* scB =
    dynamic_cast<const OpticalGeneration*>(comp_B);

  _G = alloy(scA->_G, scB->_G, xa);
}

