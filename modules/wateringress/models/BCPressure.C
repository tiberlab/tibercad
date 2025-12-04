/*  
 * This file is part of the tiberCAD module wateringress.
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
 * \file BCPressure.C
 * \brief tiberCAD wateringress module implementation.
 *
 * \note This file is part of module wateringress.
 */


#include "BCPressure.h"
#include "WIUtils.h"
#include "tibercad/base/SimulationOptions.h"

#include "tibercad/module/TiberModule.h"

using namespace libMesh;


void
BCPressure::do_init(void)
{
  _relative_humidity = get_option("relative_humidity", _relative_humidity);
}


void
BCPressure::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  double temp = SimulationOptions::temperature;
  double pressure = WIUtils::goff_gratch(temp);
  pressure *= _relative_humidity / 100;

  set_coefficients(1, 0, pressure);
}

