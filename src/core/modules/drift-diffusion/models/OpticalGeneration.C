// $Id$

#include "OpticalGeneration.h"
#include "Material.h"


#include <string>



TIBER_SUBMODEL(OpticalGeneration, optical)




void
OpticalGeneration::do_init(void)
{
  // G is a sweepable value, so check it!
  get_parameter("G", G_);
}

void
OpticalGeneration::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{


  recomb_e = recomb_h = -G_;
}


void
OpticalGeneration::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  recomb_e[0] = recomb_h[0] = 0;
  recomb_e[1] = recomb_h[1] = 0;
}


void
OpticalGeneration::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const OpticalGeneration* scA =
    dynamic_cast<const OpticalGeneration*>(comp_A);
  const OpticalGeneration* scB =
    dynamic_cast<const OpticalGeneration*>(comp_B);

  G_ = alloy(scA->G_, scB->G_, xa);
}

