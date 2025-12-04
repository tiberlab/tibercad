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
 * \file ConstantPermittivity.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */


#include "ConstantPermittivity.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"


using namespace std;

void
ConstantPermittivity::do_init(void)
{
 
  get_parameter("permittivity", _permittivity_diag, true);
  set_permittivity(_permittivity_diag);
 
  rotate();
}



void
ConstantPermittivity::read_database(void)
{
   const Database& db = get_database();
   db.set_section("permittivity");

   _permittivity_diag(0) = _permittivity_diag(1) = _permittivity_diag(2) = 1;
   db.get("permittivity", _permittivity_diag, true);
  
}
