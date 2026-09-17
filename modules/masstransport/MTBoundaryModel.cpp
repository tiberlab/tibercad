/*  
 * This file is part of the tiberCAD module masstransport.
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
 * \file MTBoundaryModel.C
 * \brief tiberCAD masstransport module implementation.
 *
 * \note This file is part of module masstransport.
 */


#include "MTBoundaryModel.h"
#include "tibercad/physics/MaterialBoundary.h"

using namespace std;

MTBoundaryModel::~MTBoundaryModel(void) = default;

MTBoundaryModel*
MTBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{
  std::string type = options.get_option("type", "pressure");
  MTBoundaryModel* mod = 
      PhysicalModel::create<MTBoundaryModel>("contact_" + type, boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "mass transport boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}


