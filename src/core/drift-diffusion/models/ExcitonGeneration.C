#include "ExcitonTransport.h"
#include "ExcitonGeneration.h"
#include "ExcitonProperties.h"
#include "DriftDiffusionProperties.h"

void
ExcitonGeneration::do_init(void)
{
  _C = get_options().get_option("C", 1e-10);
}

void
ExcitonGeneration::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();

  recomb_e = recomb_h = _C * n * p;
}

void
ExcitonGeneration::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();

  recomb_e[0] = recomb_h[0] = _C * p;
  recomb_e[1] = recomb_h[1] = _C * n;
}


void
ExcitonGeneration::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const ExcitonGeneration* scA =
    dynamic_cast<const ExcitonGeneration*>(comp_A);
  const ExcitonGeneration* scB =
    dynamic_cast<const ExcitonGeneration*>(comp_B);

  _C = alloy(scA->_C, scB->_C, xa);
}

