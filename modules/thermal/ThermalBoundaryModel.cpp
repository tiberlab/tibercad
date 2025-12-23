/*  
 * This file is part of the tiberCAD module thermal.
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
 * \file ThermalBoundaryModel.C
 * \brief tiberCAD thermal module implementation.
 *
 * \note This file is part of module thermal.
 */


#include "ThermalBoundaryModel.h"
#include "tibercad/physics/MaterialBoundary.h"

using namespace std;

ThermalBoundaryModel*
ThermalBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{


  std::string type = options.get_option("type", "heat_reservoir");
  ThermalBoundaryModel* mod =
      PhysicalModel::create<ThermalBoundaryModel>("thermal_bnd_" + type,
          boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "Thermal boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}


