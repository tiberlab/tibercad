/*
 * HoleConcentrationDependentMobility.C
 *
 *  Created on: Oct 11, 2016
 *      Author: mpatria
 */

#include "HoleConcentrationDependentMobility.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "TiberModule.h"



void
HoleConcentrationDependentMobility::read_database(void)
{
  const Database& db = get_database();
  db.set_section("mobility/hole_concentration_dependent");

  db.get("a", a_, true);
  db.get("b", b_, true);
  db.get("c", c_, true);
  db.get("d", d_, true);
}



void
HoleConcentrationDependentMobility::do_init(void)
{
}



double
HoleConcentrationDependentMobility::get_mobility(void)
{
  double N_hc = get_driftdiffusionproperties().get_hole_density();
  return d_ + (a_ - d_)/(1 + pow(N_hc/c_,b_));
}




void
HoleConcentrationDependentMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}



