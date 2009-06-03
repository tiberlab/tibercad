// $Id$

#include "ExcitonGeneration.h"

#include "ExcitonTransport.h"
#include "ExcitonProperties.h"
#include "DriftDiffusionProperties.h"


TIBER_MODULE(ExcitonGeneration, exciton_generation)


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


void
ExcitonGeneration::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const ExcitonGeneration* scA =
    dynamic_cast<const ExcitonGeneration*>(comp_A);
  const ExcitonGeneration* scB =
    dynamic_cast<const ExcitonGeneration*>(comp_B);

  C_ = alloy(scA->C_, scB->C_, xa);
}

