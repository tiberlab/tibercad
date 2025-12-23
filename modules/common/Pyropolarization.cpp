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
 * \file Pyropolarization.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */


#include "Pyropolarization.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"



using namespace std;

void
Pyropolarization::do_init(void)
{
 
  // read the Pz componenents in the crystal system
  get_parameter("Pz", _Pz, true, initializer(&Pyropolarization::_initP));
  _initP();
  
}


void
Pyropolarization::read_database(void)
{

  // Read the Pz componenents in the crystal system
  const Database& db = get_database();
  db.set_section("pyroelectricity");
  _Pz = db.get("Pz", _Pz);

}


void
Pyropolarization::_initP(void)
{
  PolarizationModel::do_init();

  set_polarization(libMesh::RealVectorValue(0.0, 0.0, _Pz));
  rotate();
}
