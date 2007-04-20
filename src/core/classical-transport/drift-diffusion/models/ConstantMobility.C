// $Id$

#include "ConstantMobility.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "getpot.h"


TIBER_MODULE(ConstantMobility, constant)


void
ConstantMobility::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());
  
  std::string s("mu_max_");
  s += get_carrier_type();
  mu0_ = data(s.c_str(), mu0_);

  s = "exponent_";
  s += get_carrier_type();
  exp_ = data(s.c_str(), exp_);
}



void
ConstantMobility::do_init(void)
{
  mu0_ = get_options().get_option("mu", mu0_);
  std::string s("mu_");
  s += get_carrier_type();
  mu0_ = get_material()->get_options().get_option(s, mu0_);
}



double
ConstantMobility::get_mobility(void)
{
  double T = get_driftdiffusionproperties().get_lattice_temperature();
  return mu0_ * std::pow(T / T0, -exp_);
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

  mu0_ = alloy(scA->mu0_, scB->mu0_, xa);
  exp_ = alloy(scA->exp_, scB->exp_, xa);
}

