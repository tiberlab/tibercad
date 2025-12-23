/*  
 * This file is part of the tiberCAD module negf.
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
 * \file EtbModel.C
 * \brief tiberCAD negf module implementation.
 *
 * \note This file is part of module negf.
 */


#include "EtbModel.h"
#include "tibercad/io/Messages.h"

#include "tibercad/module/TiberModule.h"

EtbModel::EtbModel(const ModelOptions& options)
   : HamiltonianModel(options)
{
   //_inv_mass_crys(0);
   //_inv_mass(0);
}


void
EtbModel::do_init(void)
{
  _model = "etb";

  _simulation = get_option("simulation","none");

  if (_simulation=="none")
  {
    Messages::error("in hamiltonian submodel etb must define simulation");
    exit(1);
  }

  _degeneracy = 1;
}

