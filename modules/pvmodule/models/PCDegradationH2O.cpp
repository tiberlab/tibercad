/*  
 * This file is part of the tiberCAD module pvmodule.
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
 * \file PCDegradationH2O.C
 * \brief tiberCAD pvmodule module implementation.
 *
 * \note This file is part of module pvmodule.
 */

#include "PCDegradationH2O.h"
#include "tibercad/module/SimulationInterface.h"

#include "tibercad/module/TiberModule.h"

PCDegradationH2O::PCDegradationH2O(const ModelOptions& options)
  : Photocurrent(options)
{
}


void
PCDegradationH2O::do_init(void)
{
  _initial_current = get_option("initial_photocurrent", _initial_current);

  _RH_ref = get_option("reference_humidity", _RH_ref);

  _exponent = get_option("exponent", _exponent);

  std::string water_ingress_sim = get_option("relative_humidity", "");
  _humidity_model = SimulationInterface::find_solution_provider(water_ingress_sim, "RelativeHumidity");
}


double
PCDegradationH2O::degradation_factor(double humidity) const
{
  double a = 1.0 + std::pow(humidity / _RH_ref, _exponent);

  return 1.0 / a;
}


double
PCDegradationH2O::do_get_photocurrent(const libMesh::Elem* elem,
                                      const libMesh::Point& p) const
{
  double curr = _initial_current;

  double humidity = 0;
  if (_humidity_model.first != nullptr)
  {
    if (_humidity_model.first->get_solution(elem, _humidity_model.second, humidity, p, false))
      curr *= degradation_factor(humidity);
  }

  return curr;
}
