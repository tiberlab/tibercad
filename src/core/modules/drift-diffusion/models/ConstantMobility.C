// $Id$

#include "ConstantMobility.h"
#include "DriftDiffusionProperties.h"

#include "Database.h"



TIBER_SUBMODEL(ConstantMobility, constant)


void
ConstantMobility::read_database(void)
{
  Database& db = get_database();
  db.set_section("mobility/constant");
  
  std::vector<double> data(2, 0);
  db.get("mu_max", data, true);
  mu0_ = get_carrier_type() == 'e' ? data[0] : data[1];

  data = std::vector<double>(2, 0);
  db.get("exponent", data);
  exp_ = get_carrier_type() == 'e' ? data[0] : data[1];

}



void
ConstantMobility::do_init(void)
{
  mu0_ = get_parameter("mu", mu0_);
  // we allow also mu_e and mu_h
  mu0_ = get_parameter(std::string("mu_") + get_carrier_type(), mu0_);
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
ConstantMobility::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  const ConstantMobility* scA =
    dynamic_cast<const ConstantMobility*>(comp_A);
  const ConstantMobility* scB =
    dynamic_cast<const ConstantMobility*>(comp_B);

  mu0_ = alloy(scA->mu0_, scB->mu0_, xa);
  exp_ = alloy(scA->exp_, scB->exp_, xa);
}

