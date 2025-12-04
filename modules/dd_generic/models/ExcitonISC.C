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
 * \file ExcitonISC.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "ExcitonISC.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/io/Database.h"
#include "tibercad/io/Messages.h"

#include "tibercad/module/TiberModule.h"

 
using namespace std;



void
ExcitonISC::read_database(void)
{
  const Database& db = get_database();
}



void
ExcitonISC::do_init(void)
{
  RecombinationModelInterface::do_init();


  if (get_carrier_names().size() != 2)
    throw InitFailedException("Exciton ISC model needs exactly "
        "two recombining carriers");

  get_parameter("C", _C);

  const DriftDiffusionProperties& dd = get_bulk_driftdiffusionproperties();
  vector<double> spin;
  for (auto name : get_carrier_names())
  {
    //if (!dd.get_carrier_properties(name)->is_exciton())
    //  throw InitFailedException("Recombination '" + get_default_name() + ": carrier '" + name + "' is not an exciton");

    spin.push_back( dd.get_carrier_properties(name)->get_spin() );
  }

  if (!(((spin[0] == 1.0) && (spin[1] == 0.0)) ||
        ((spin[0] == 0.0) && (spin[1] == 1.0))))
    throw InitFailedException("Recombination '" + get_default_name() +
        ": carriers for ISC model must be a singlet and a triplet exciton");


}




void
ExcitonISC::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  ID id1 = this->get_carrier_ids()[0];
  ID id2 = this->get_carrier_ids()[1];

  double E1  = dd.get_carrier_properties(id1)->get_band_edge();
  double E2  = dd.get_carrier_properties(id2)->get_band_edge();

  if (E1 < E2)
  {
    swap(id1, id2);
    swap(E1, E2);
  }

  double kT = dd.get_lattice_temperature();
  double x1  = dd.get_q_density(id1);
  double x2  = dd.get_q_density(id2);
  double N2  = dd.get_carrier_properties(id2)->get_maximum_density();
  double dx1 = dd.get_q_density_derivative(id1);
  double dx2 = dd.get_q_density_derivative(id2);
  double f1 = dd.get_q_fermi_potential(id1);
  double f2 = dd.get_q_fermi_potential(id2);

  double beta = 1.0/kT;
  double exponential = exp( beta*(f1-f2) );
  double stat = 1.0 - exponential;

  double x2f = x2/N2;

  double rate = _C * stat * x1 * (1 + x2f);
  double der1 = - _C * (1 + x2f) * (dx1 * stat + beta * x1 * exponential);
  double der2 = - _C * x1 * (dx2/N2 * stat - beta * (1 + x2f) * exponential);

  R[id1] = rate;
  R[id2] = -rate;

  dPotentials[id1][id1] =  der1;
  dPotentials[id1][id2] =  der2;
  dPotentials[id2][id1] = -der1;
  dPotentials[id2][id2] = -der2;

}

void
ExcitonISC::do_reinit(void)
{

}
