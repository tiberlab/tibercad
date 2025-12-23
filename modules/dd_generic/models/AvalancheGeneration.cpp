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
 * \file AvalancheGeneration.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "AvalancheGeneration.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/math/TiberMath.h"

#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"


#include "tibercad/module/TiberModule.h"


using namespace std;

AvalancheGeneration::AvalancheGeneration(const ModelOptions& options) :
  RecombinationModelInterface(options),
  _w0(0.05)
{
}



AvalancheGeneration::~AvalancheGeneration(void)
{
}


void
AvalancheGeneration::read_database(void)
{
  const Database& db = get_database();

  db.set_section("recombination/AvalancheGeneration");

  _w0 = db.get("w0", _w0);

  db.get("a", _a_param);
  db.get("b", _b_param);
}




void
AvalancheGeneration::do_init(void)
{
  RecombinationModelInterface::do_init();

  unsigned int nc = get_carrier_names().size();

  _a_param.resize(nc, 0);
  _b_param.resize(nc, 0);

  get_option("a", _a_param);
  get_option("b", _b_param);
  get_parameter("w0", _w0);

}



void
AvalancheGeneration::calculate_rate_and_derivatives(std::vector<double>& R,
    std::vector<std::vector<double>>& dPotentials)
{
  const vector<ID>& ids = this->get_carrier_ids();
  unsigned int nc = ids.size();

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  //double Efn = -dd.get_q_fermi_potential(id1);
  //double Efp = -dd.get_q_fermi_potential(id2);
  double kT = dd.get_lattice_temperature();

  vector<double> alpha(nc, 0.0);

  // we add a small value to be able to divide by F in any case
  double F = dd.get_electric_field().norm() + 1e-6;

  for (unsigned int i = 0; i < nc; ++i)
  {
    ID id = ids[i];

    libMesh::RealGradient grad;
    dd.get_grad_fermi(id, grad);

    F = grad.norm() + 1e-6;

    double arg = -_b_param[i] / F;
    if (arg < -15)
      alpha[i] = 0;
    else
      alpha[i] = _a_param[i] * exp(arg);


    double flux = dd.get_q_conductivity(id) * grad.norm();

    for (unsigned int j = 0; j < nc; ++j)
      R[ids[j]] -= alpha[i] * flux;
  }

  //dPotentials[id1][id1] = dR1;
  //dPotentials[id1][id2] = dR2;
  //dPotentials[id2][id1] = dR1;
  //dPotentials[id2][id2] = dR2;
  //dPotentials[id1][dd.n_known_carriers()] = dR0;
  //dPotentials[id2][dd.n_known_carriers()] = dR0;
}

