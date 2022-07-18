/*
 * VariableRangeHopping.C
 *
 *  Created on: 25 Jan 2022
 *      Author: miesu
 */
#include "VariableRangeHopping.h"
#include "DriftDiffusionProperties.h"
#include "Constants.h"

#include "TiberModule.h"






void
VariableRangeHopping::do_init(void)
{
  _mu_0 = get_option("mu_0" , _mu_0);
  _hopping_distance= get_option("r" , _hopping_distance);
  _temp_hopping_conduction = get_option("Th" , _temp_hopping_conduction);
  _gamma = get_option("gamma" , _gamma);

}






double VariableRangeHopping::get_mobility(void)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double E = dd.get_electric_field().norm();

  double kT = get_driftdiffusionproperties().get_lattice_temperature();

  double Kb = Constants::kb;  //  Boltzmann constant in eV / K



  double a =  Kb * _temp_hopping_conduction;

  double b = kT + E * _hopping_distance;

  double argument_of_exp = std::pow( ( a/b ) , _gamma);

  double exp_mu = std::exp( - argument_of_exp );



  double mu = _mu_0 * exp_mu;

  return mu;
}






void
VariableRangeHopping::get_derivative_grad_potential(libMesh::RealGradient& dm)
{
  dm.zero();

  double kT = get_driftdiffusionproperties().get_lattice_temperature();

  double Kb = Constants::kb;  //  Boltzmann constant in eV / K

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  const libMesh::RealGradient& F = dd.get_electric_field();

  double E = dd.get_electric_field().norm();



  double numerator_of_A = _mu_0 * _hopping_distance * _gamma;

  double denominator_of_A = kT + _hopping_distance * E;

  double A =  numerator_of_A / denominator_of_A;



  double numerator_of_B = Kb * _temp_hopping_conduction;

  double denominator_of_B = kT + E * _hopping_distance;

  double B = numerator_of_B / denominator_of_B;



  double B_gamma = std::pow(B, _gamma);

  double exp_minus_B_gamma = std::exp( -B_gamma );



  double dmu = (E > 1e-3) ?
               ( -1.0 ) * A * B_gamma * exp_minus_B_gamma *( 1.0 / E ) : 0.0;

  dm = F * dmu;
}
