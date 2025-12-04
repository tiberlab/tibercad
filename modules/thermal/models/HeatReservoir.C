/*  
 * This file is part of the tiberCAD module thermal.
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
 * \file HeatReservoir.C
 * \brief tiberCAD thermal module implementation.
 *
 * \note This file is part of module thermal.
 */


#include "HeatReservoir.h"

#include "tibercad/module/TiberModule.h"
#include "tibercad/module/SimulationInterface.h"


void
HeatReservoir::do_init(void)
{
  get_parameter("temperature", _temperature);

  std::string sim_name = get_option("host_simulation", "");

  _host_sim = 0;
  if (sim_name != "") {
    SimulationInterface* sim = SimulationInterface::find_simulation(sim_name);

    if (sim == nullptr) {
      std::string msg("External thermal host simulation " + std::string(sim_name) + " not found");
      throw InitFailedException(msg);
    }
    _host_sim = sim->get_id();
  }
}

void
HeatReservoir::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  if (_host_sim > 0)
  {
    SimulationInterface* sim = SimulationInterface::get_simulation(_host_sim);
    if (sim->is_solved())
    {
      ID T_id = sim->get_solution_id("LatticeTemp");

      double T = 0.0;
      sim->get_solution(elem, T_id, T, point);
      _temperature = T;
    }
  }

  set_coefficients(1.0, 0.0, _temperature);
}

