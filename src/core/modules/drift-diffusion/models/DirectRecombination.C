// $Id$

#include "DirectRecombination.h"
#include "DriftDiffusionProperties.h"


#include "getpot.h"



TIBER_SUBMODEL(DirectRecombination, direct)




void
DirectRecombination::read_database(void)
{
  Database& db = get_database();
  db.set_section("recombination/direct");

  C_ = db.get("C", C_);

}



void
DirectRecombination::do_init(void)
{
  get_parameter("C", C_);
}



void
DirectRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();

  recomb_e = recomb_h = C_ * (n * p - ni * ni);
}



void
DirectRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();

  recomb_e[0] = recomb_h[0] = C_ * p; // dR/dn
  recomb_e[1] = recomb_h[1] = C_ * n; // dR/dp
}



void
DirectRecombination::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const DirectRecombination* scA =
    dynamic_cast<const DirectRecombination*>(comp_A);
  const DirectRecombination* scB =
    dynamic_cast<const DirectRecombination*>(comp_B);

  C_ = alloy(scA->C_, scB->C_, xa);
}

