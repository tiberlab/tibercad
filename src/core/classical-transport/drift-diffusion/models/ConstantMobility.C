// $Id$

#include "ConstantMobility.h"


void
ConstantMobility::do_init(void)
{
  _mu0 = get_options().get_option("mu0", 1000.0);
}

double
ConstantMobility::get_mobility(void)
{
  return _mu0;
}


void
ConstantMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}


void
ConstantMobility::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  const ConstantMobility* scA =
    dynamic_cast<const ConstantMobility*>(comp_A);
  const ConstantMobility* scB =
    dynamic_cast<const ConstantMobility*>(comp_B);

  _mu0 = alloy(scA->_mu0, scB->_mu0, xa);
}

