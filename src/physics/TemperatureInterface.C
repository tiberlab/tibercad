/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file TemperatureInterface.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/TemperatureInterface.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/base/InitFailedException.h"

#include "elem.h"



TemperatureInterface::TemperatureInterface(void)
  : _simulation(NULL),
    _variable_name("LatticeTemp"),
    _id(INVALID_ID)
{
}


TemperatureInterface::TemperatureInterface(const TemperatureInterface& rhs)
  : _simulation(rhs._simulation),
    _variable_name(rhs._variable_name),
    _id(rhs._id)
{
}

TemperatureInterface&
TemperatureInterface::operator=(const TemperatureInterface& rhs)
{
  _simulation = rhs._simulation;
  _variable_name = rhs._variable_name;
  _id = rhs._id;
  return(*this);
}



bool
TemperatureInterface::set_simulation(const std::string& name)
{
  bool answer = false;
  if (name != "")
  {
    _simulation = SimulationInterface::find_simulation(name);
    if (_simulation == NULL)
      throw InitFailedException("No such simulation found: " + name);

    _id = _simulation->get_solution_id(_variable_name);

    if (_id == INVALID_ID)
      throw InitFailedException("Simulation " + name +
          " has no variable '" + _variable_name + "'");

    answer = true;
  }

  return answer;
}


void
TemperatureInterface::get_temperature(const Elem* elem,
    std::vector<double>& temperatures)
{
  assert(elem != NULL);

  int nn = elem->n_nodes();
  std::vector<Point> nodes(nn);

  for (unsigned int i = 0; i < nn; ++i)
     nodes[i] = elem->master_point(i);

  get_temperature(elem, nodes, temperatures, true);
}


void
TemperatureInterface::get_temperature(const Elem* elem,
    const std::vector<Point>& p, std::vector<double>& temperatures, bool refcoord)
{
  assert(elem != NULL);

  int nn = p.size();

  temperatures.resize(nn);

  if (_simulation == NULL)
  {
    for (int i = 0; i < nn; i++)
      temperatures[i] = SimulationOptions::temperature;
  }
  else
  {
    std::map<ID, std::vector<double> > temp;
    temp[_id] = std::vector<double>();

    if (_simulation->get_solution(elem, temp, p, refcoord))
      temperatures = temp[_id];
    else
      for (int i = 0; i < nn; i++)
        temperatures[i] = SimulationOptions::temperature;
  }
}


double
TemperatureInterface::get_temperature(const Elem* elem, const Point& p, bool refcoord)
{
  std::vector<Point> ps(1, p);
  std::vector<double> temp(1);

  get_temperature(elem, ps, temp, refcoord);

  return temp[0];
}

bool
TemperatureInterface::is_solved(void) const
{
  return (_simulation->is_solved());
}
