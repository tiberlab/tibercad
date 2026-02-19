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
 * \file WIModel.C
 * \brief tiberCAD wateringress module implementation.
 *
 * \note This file is part of module wateringress.
 */


#include "WIModel.h"

#include "tibercad/physics/Material.h"

using namespace std;


TiberModelObject*
WIModel::_create(const ModelOptions& options)
{
  return new WIModel(options);
}


void
WIModel::_destroy(TiberModelObject* p)
{
  delete p;
}


WIModel*
WIModel::create(const Material* mat, const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  WIModel* pm = NULL;

  if (type == "default")
    // we create the default model from explicit creation method
    pm = PhysicalModel::create<WIModel>(_create, _destroy, mat, options);
  else
  {
    // there is no such model, at the moment
    type = "bulk_" + type;
    pm = PhysicalModel::create<WIModel>(type, mat, options);
  }

  return(pm);
}



void
WIModel::do_init(void)
{
  // we read it in g/m^3/Pa
  _solubility = get_option("solubility", _solubility);

  // we read it in m^2/S
  _diffusivity = get_option("diffusivity", _diffusivity);
}


void
WIModel::calculate(const Elem* elem, const Point& point)
{
 
}

void
WIModel::prepare_submodels(void)
{
}
