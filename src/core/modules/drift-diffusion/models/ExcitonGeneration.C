// $Id$

#include "ExcitonGeneration.h"

#include "ExcitonTransport.h"
#include "ExcitonProperties.h"
#include "DriftDiffusionProperties.h"


TIBER_MODULE(ExcitonGeneration, recombination, exciton_generation)


void
ExcitonGeneration::do_init(void)
{
  C_ = get_option("C", 1e-10);
}

void
ExcitonGeneration::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();

  recomb_e = recomb_h = C_ * n * p;
}

void
ExcitonGeneration::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();

  recomb_e[0] = recomb_h[0] = C_ * p;
  recomb_e[1] = recomb_h[1] = C_ * n;
}
