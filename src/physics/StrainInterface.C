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
 * \file StrainInterface.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/StrainInterface.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/base/InitFailedException.h"
#include "tibercad/base/RuntimeException.h"
#include "tibercad/io/Messages.h"
#include "tibercad/math/Tensor2.h"

#include "libmesh/elem.h"




StrainInterface::StrainInterface(void)
  : _simulation(NULL)
{
}



bool
StrainInterface::set_simulation(const std::string& name)
{
  bool answer = false;
  
  if (name != "")
  {
    _simulation = SimulationInterface::find_simulation(name);
    if (_simulation == NULL)
      throw InitFailedException("No such simulation found: " + name);

    _strain_id = _simulation->get_solution_id("Strain");
    _strain_cryst_id = _simulation->get_solution_id("StrainCrystal");
    _stress_id = _simulation->get_solution_id("Stress");
    _stress_cryst_id = _simulation->get_solution_id("StressCrystal");

    if (_strain_id == INVALID_ID)
      Messages::warning("Simulation " + _simulation->get_name() +
          " does not have variable \'Strain\'.");

    if (_strain_cryst_id == INVALID_ID)
      Messages::warning("Simulation " + _simulation->get_name() +
          " does not have variable \'StrainCrystal\'.");

    if (_stress_id == INVALID_ID)
      Messages::warning("Simulation " + _simulation->get_name() +
          " does not have variable \'Stress\'.");

    if (_stress_cryst_id == INVALID_ID)
      Messages::warning("Simulation " + _simulation->get_name() +
          " does not have variable \'StressCrystal\'.");

    
    if ((_strain_id == INVALID_ID) && (_strain_cryst_id == INVALID_ID) &&
        (_stress_id == INVALID_ID) && (_stress_cryst_id == INVALID_ID))
      throw InitFailedException("Simulation " + name +
          " is missing strain related variables");

    answer = true;
  }

  return answer;
}



void
StrainInterface::get_strain(const Elem* elem, const Point& point, Tensor2& strain)
{
  if (_simulation)
  {
    if (_strain_id == INVALID_ID)
      throw RuntimeException("Simulation " + _simulation->get_name() +
          " does not have variable \'Strain\'.");

    _get_data(elem, point, strain, _strain_id);
  }
}


void
StrainInterface::get_crystal_strain(const Elem* elem, const Point& point, Tensor2& strain)
{
  if (_simulation)
  {
    if (_strain_id == INVALID_ID)
      throw RuntimeException("Simulation " + _simulation->get_name() +
          " does not have variable \'StrainCrystal\'.");

    _get_data(elem, point, strain, _strain_cryst_id);
  }
}


void
StrainInterface::get_stress(const Elem* elem, const Point& point, Tensor2& stress)
{
  if (_simulation)
  {
    if (_strain_id == INVALID_ID)
      throw RuntimeException("Simulation " + _simulation->get_name() +
          " does not have variable \'Stress\'.");

    _get_data(elem, point, stress, _stress_id);
  }
}


void
StrainInterface::get_crystal_stress(const Elem* elem, const Point& point, Tensor2& stress)
{
  if (_simulation)
  {
    if (_strain_id == INVALID_ID)
      throw RuntimeException("Simulation " + _simulation->get_name() +
          " does not have variable \'StressCrystal\'.");

    _get_data(elem, point, stress, _stress_cryst_id);
  }
}


inline
void
StrainInterface::_get_data(const Elem* elem, const Point& point,
    Tensor2& data, ID id)
{
  std::vector<Point> p(1);
  p[0] = point;
  std::vector<double> values(6);

  if (_simulation->get_solution(elem, id, values, p))
  {
    data(1,1) = values[0];
    data(2,2) = values[1];
    data(3,3) = values[2];
    data(2,1) = data(1,2) = values[3];
    data(3,2) = data(2,3) = values[4];
    data(3,1) = data(1,3) = values[5];
  }
}


