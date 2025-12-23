/*  
 * This file is part of the tiberCAD module common.
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
 * \file DDHeatSource.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */


#include "DDHeatSource.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"





using namespace std;


DDHeatSource::DDHeatSource(const ModelOptions& options):HeatSourceModel(options)
{

}

void
DDHeatSource::do_init(void)
{

  string dd_simul_name = get_options().get_option("transport_simulation", "driftdiffusion");
  _simul = SimulationInterface::find_simulation(dd_simul_name);

  if ( _simul == NULL)
   throw InitFailedException("Could not find " + dd_simul_name);

  vector<string> heat_sources;
  get_options().get_option("heat_sources", heat_sources);

  if (heat_sources.empty())
  {
    ID_set.insert(_simul->get_solution_id("eJoule"));
    ID_set.insert(_simul->get_solution_id("hJoule"));
    ID_set.insert(_simul->get_solution_id("RecombHeat"));
    ID_set.insert(_simul->get_solution_id("ePeltier"));
    ID_set.insert(_simul->get_solution_id("hPeltier"));
  }
  else
  {
    for (auto& src : heat_sources)
    {
      ID_set.insert(_simul->get_solution_id(src));
    }
  }

  if (ID_set.count(INVALID_ID))
   throw InitFailedException("Invalid heat sources defined in heat source model");
}

void
DDHeatSource::calculate(const Elem* elem, const Point& point)
{

  std::map<ID, vector< double > > solution;
  for (auto& id : ID_set)
    solution[id].resize(0);

  double heat_source = 0.0;
  vector<Point> h_point(1);
  h_point[0] = point;

  if (_simul->get_solution(elem, solution, h_point))
  {
    for (auto& src : solution)
      heat_source += src.second[0];
  }
  set_heat_source(heat_source * 1E6);

}


 
