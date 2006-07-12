// $Id$

#include "SimpleSemiconductorModel.h"

#include "point.h"
#include "elem.h"

#include <iostream>

using namespace DriftDiffusionDefs;


SimpleSemiconductorModel::SimpleSemiconductorModel(void)
  : Parent(),
    _is_prepared(false)
{
}

/*
SimpleSemiconductorModel::SimpleSemiconductorModel(
    const SimpleSemiconductorModel& model)
  : Parent(model),
    _is_prepared(model._is_prepared)
{
}
*/

void
SimpleSemiconductorModel::prepare_element_data(void)
{
  if (!_is_prepared)
  {
    double kT = SimulationOptions::T * Constants::k_B;
    electron_vt = hole_vt = kT;

    calculate_equilibrium_properties(BOTH, SimulationOptions::T);

    _is_prepared =  true;
  }
}



void
SimpleSemiconductorModel::calculate_all(double potential,
    double fermi_e, double fermi_h, const Point& coord)
{
  int coupling = get_coupling_type();

  // in this simple model all temperatures are equal
  double kT = electron_vt;

  // call the method of the parent class
  Parent::calculate_all(potential, fermi_e, fermi_h, coord);

  double n = electron_density;
  double dn = electron_density_derivative;
  double p = hole_density;
  double dp = hole_density_derivative;

  // 4.) mobilities / conductivities
  // For both statistics:
  // 
  //   mu_n = e * D_n * (1 / n) * (dn / dEf_e)
  //
  // NOTE: kT := kB * T / e includes already e
  //if (coupling & DriftDiffusionDefs::ELECTRONS)
  //{
    electron_conductivity = electron_mobility * n;
    electron_conductivity_derivatives[0] = electron_mobility * dn;
    electron_conductivity_derivatives[1] = electron_mobility * dn;
  //}
  //if (coupling & DriftDiffusionDefs::HOLES)
  //{
    hole_conductivity = hole_mobility * p;
    hole_conductivity_derivatives[0] = hole_mobility * dp;
    hole_conductivity_derivatives[2] = hole_mobility * dp;
  //}
}


