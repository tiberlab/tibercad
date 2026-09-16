/*  
 * This file is part of the tiberCAD module dd_generic.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file VariableRangeHopping.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */

/*
 * VariableRangeHopping.C
 *
 *  Created on: 25 Jan 2022
 *      Author: miesu
 */
#include "VariableRangeHopping.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/physics/Constants.h"

#include "tibercad/module/TiberModule.h"



VariableRangeHopping::VariableRangeHopping(const ModelOptions& options)
 : MobilityModelInterface(options),
   _mu_0(100.0),
  _gamma(1.0/3.0),
  _hopping_distance(1.0e-6),
  _temp_hopping_conduction(1200.0)
{
}





VariableRangeHopping::~VariableRangeHopping(void)
{
}






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
