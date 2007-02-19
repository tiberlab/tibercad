// $Id$

#include "DirectRecombination.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "getpot.h"

void
DirectRecombination::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  _C = data("C", _C);
  
}



void
DirectRecombination::do_init(void)
{
  _C = get_options().get_option("C", _C);
  _C = (get_material()->get_options()).get_option("C", _C);
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
  double p  = dd.get_hole_density();

  recomb_e[0] = recomb_h[0] = _C * p; // dR/dn
  recomb_e[1] = recomb_h[1] = _C * n; // dR/dp
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

