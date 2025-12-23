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
 * \file TmmBoundaryModel.C
 * \brief tiberCAD tmm module implementation.
 *
 * \note This file is part of module tmm.
 */


#include "TmmBoundaryModel.h"
#include "tibercad/physics/MaterialBoundary.h"
#include "Tmm.h"

using namespace std;

TmmBoundaryModel*
TmmBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{

  std::string type = options.get_option("type", "incidentwave");
  TmmBoundaryModel* mod =
    PhysicalModel::create<TmmBoundaryModel>("tmm_bnd_" + type,
      boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "Tmm boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}

std::string
TmmBoundaryModel::read_type(void) {
  return typer;
}

void
TmmBoundaryModel::write_type(std::string str) {
  typer = str;
}


void
TmmBoundaryModel::set_elements(double a0, double a1, double a2, double a3) {
  _mmm00 = a0;
  _mmm01 = a1;
  _mmm10 = a2;
  _mmm11 = a3;
}

void
TmmBoundaryModel::set_dipole_elements(double a0, double a1) {
  __kr = a0;
  __steps = a1;
}


double
TmmBoundaryModel::get_element(int elm) {
  switch (elm)
  {
  case 0:
    return(_mmm00);
    break;
  case 1:
    return(_mmm01);
    break;
  case 2:
    return(_mmm10);
    break;
  case 3:
    return(_mmm11);
    break;
  default:
    return(0);
    break;

  }
}
double
TmmBoundaryModel::get_kr() {
  return(__kr);
}
double
TmmBoundaryModel::get_steps() {
  return(__steps);
}

