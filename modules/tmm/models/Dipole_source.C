/*  
 * This file is part of the tiberCAD module tmm.
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
 * \file Dipole_source.C
 * \brief tiberCAD tmm module implementation.
 *
 * \note This file is part of module tmm.
 */

/*
 * incidentwave.C
 *
 *  Created on: 4 Oct 2021
 *      Author: pamiri
 */

#include "Dipole_source.h"
#include "tibercad/module/TiberModule.h"

void Dipole_source::do_init(void)
{
  write_type("Dipole Source");
  get_parameter("kr_ratio", _kr);
  get_parameter("steps", _steps);
  set_dipole_elements(_kr,_steps);
 // typer = "Incident Wave";
}
void Dipole_source::Calculate_M_Matrix(void)
{
  ;;
}

