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
 * \file LangevinRecombination.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */



#include "LangevinRecombination.h"
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
LangevinRecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("permittivity");

  _er = db.get("permittivity", _er);
}

void
LangevinRecombination::do_init(void)
{
  get_parameter("gamma", _gamma);
}



double
LangevinRecombination::calculate_rate_and_derivatives(std::vector<double>& dPotentials)
{
  const ID id1 = this->get_carrier_ids()[0];
  const ID id2 = this->get_carrier_ids()[1];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  double Efn = -dd.get_q_fermi_potential(id1);
  double Efp = -dd.get_q_fermi_potential(id2);
  double kT = dd.get_lattice_temperature();

  double n  = dd.get_q_density(id1);
  double p  = dd.get_q_density(id2);
  double dn  = dd.get_q_density_derivative(id1);
  double dp  = dd.get_q_density_derivative(id2);
  double qn = dd.get_carrier_properties(id1)->get_charge();
  double qp = dd.get_carrier_properties(id2)->get_charge();

  double mun = dd.get_q_mobility(id1);
  double mup = dd.get_q_mobility(id2);

  double exponential = exp((Efp - Efn) / kT);
  double stat_fac = 1.0 - exponential;
  double prefactor = _gamma * Constants::e * 100 / (_er * Constants::e0) * (mun + mup);

  double g = prefactor * n * p;

  double R = g * stat_fac;

  double dR1 = prefactor * p * (dn * stat_fac - 1/kT * n * exponential);
  double dR2 = prefactor * n * (dp * stat_fac + 1/kT * p * exponential);

  dPotentials[id1] = dR1;
  dPotentials[id2] = dR2;
  dPotentials[dd.n_known_carriers()] = stat_fac * prefactor * (p * dn + n * dp);

  return R;
}

