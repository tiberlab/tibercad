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
 * \file ExcitonDecay.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "ExcitonDecay.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/io/Database.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/io/Messages.h"

#include "mesh_base.h"
#include "quadrature.h"

#include "tibercad/module/TiberModule.h"


using namespace std;


void
ExcitonDecay::read_database(void)
{
  const Database& db = get_database();
}



void
ExcitonDecay::do_init(void)
{
  RecombinationModelInterface::do_init();

  if (get_carrier_ids().size() != 1)
    throw InitFailedException("ExcitonDecay model needs exactly one carrier");

  _tau = get_option("tau", _tau);

}




void
ExcitonDecay::calculate_rate_and_derivatives(std::vector<double>& R,
    std::vector<std::vector<double>>& dPotentials)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();


  double kT = dd.get_lattice_temperature();

  double id = get_carrier_ids()[0];
  unsigned int np = dd.n_known_carriers();

  double x = dd.get_q_density(id);
  double dx = dd.get_q_density_derivative(id);
  double Efx = -dd.get_q_fermi_potential(id);

  double beta = 1.0 / kT;
  double exponential = exp(-beta * Efx);
  double stat = 1.0 - exponential;

  double rate = stat * x / _tau;

  R[id] = rate;

  dPotentials[id][id] = - (beta * x * exponential + dx * stat) / _tau;

}
