/*  
 * This file is part of the tiberCAD module elasticity.
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
 * \file ElasticityBoundaryModel.C
 * \brief tiberCAD elasticity module implementation.
 *
 * \note This file is part of module elasticity.
 */


#include "ElasticityBoundaryModel.h"
#include "tibercad/physics/MaterialBoundary.h"

using namespace std;


ElasticityBoundaryModel*
ElasticityBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{
 

  std::string type = options.get_option("type", "clamp");

  ElasticityBoundaryModel* mod = PhysicalModel::create<ElasticityBoundaryModel>("ebnd_" + type,
      boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "Elasticity boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;

}


