/*  
 * This file is part of the tiberCAD module boltzmann.
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
 * \brief tiberCAD boltzmann module implementation.
 *
 * \note This file is part of module boltzmann.
 */


#include "HeatReservoir.h"

#include "tibercad/module/TiberModule.h"


void
HeatReservoir::do_init(void)
{
  get_parameter("temperature", _temperature);
}

void
HeatReservoir::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  set_coefficients(1.0, 0.0, _temperature);
}

