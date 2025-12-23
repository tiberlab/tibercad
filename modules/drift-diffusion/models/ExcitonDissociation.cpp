/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file ExcitonDissociation.C
 * \brief tiberCAD driftdiffusion module implementation.
 *
 * \note This file is part of module driftdiffusion.
 */


#include "ExcitonDissociation.h"

#include "tibercad/module/SimulationInterface.h"
#include "DriftDiffusionProperties.h"

#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"



void
ExcitonDissociation::do_init(void)
{
  d_ = get_option("damping", 1.0);

  std::string ex = get_option("exciton_simulation", "");

  // find the exciton simulation to use
  _exciton_sim = SimulationInterface::find_simulation(ex);

  if (_exciton_sim == NULL)
  {
    std::string msg("ExcitonDissociation: Simulation " +
        std::string(ex) + " not found");
    throw InitFailedException(msg);
  }

  _Rdiss_id = _exciton_sim->get_solution_id("dissociation");
}




void
ExcitonDissociation::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{

  DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  const Elem* el = dd.get_element();

  recomb_e = 0.0;

  // we only use the exciton simulation if it has been solved before
  if (_exciton_sim->is_solved())
  {
    double x = 0.0;
    bool succ = _exciton_sim->get_solution(el, _Rdiss_id, x, dd.get_coordinates());
    if (succ)
      recomb_e = -d_ * x;
  }

  recomb_h = recomb_e;
}




void
ExcitonDissociation::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  recomb_e[0] = recomb_h[0] = 0.0;
  recomb_e[1] = recomb_h[1] = 0.0;
}





