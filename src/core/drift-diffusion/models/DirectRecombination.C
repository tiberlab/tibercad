// $Id$

#include "DirectRecombination.h"
#include "DriftDiffusionProperties.h"


void
DirectRecombination::do_init(void)
{
  _C = get_options().get_option("C", 1e-10);
}

void
DirectRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();

  recomb_e = recomb_h = _C * (n * p - ni * ni);
}

void
DirectRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  
  double n  = dd.get_electron_density();
  double dn  = dd.get_electron_density_derivative();
  double p  = dd.get_hole_density();
  double dp  = dd.get_hole_density_derivative();

  double a = _C * dn * p;
  double b = _C * n * dp; 

  recomb_e[0] = recomb_h[0] = a + b;
  recomb_e[1] = recomb_h[1] = -a;
  recomb_e[2] = recomb_h[2] = -b;
}


void
DirectRecombination::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const DirectRecombination* scA =
    dynamic_cast<const DirectRecombination*>(comp_A);
  const DirectRecombination* scB =
    dynamic_cast<const DirectRecombination*>(comp_B);

  _C = alloy(scA->_C, scB->_C, xa);
}

