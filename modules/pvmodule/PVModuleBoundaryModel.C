/*  
 * This file is part of the tiberCAD module pvmodule.
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
 * \file PVModuleBoundaryModel.C
 * \brief tiberCAD pvmodule module implementation.
 *
 * \note This file is part of module pvmodule.
 */


#include "PVModuleBoundaryModel.h"
#include "tibercad/io/Messages.h"

using namespace std;


PVModuleBoundaryModel::PVModuleBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}


TiberModelObject*
PVModuleBoundaryModel::_create(const ModelOptions& options, const void*)
{
  return new PVModuleBoundaryModel(options);
}


void
PVModuleBoundaryModel::_destroy(TiberModelObject* p)
{
  delete p;
}


PVModuleBoundaryModel*
PVModuleBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  PVModuleBoundaryModel* pm = NULL;

  if (type == "default")
    // we create the default model from explicit creation method
    pm = PhysicalModel::create<PVModuleBoundaryModel>(_create, _destroy, boundary, options);
  else
  {
    // there is no such model, at the moment
    type = "contact_" + type;
    pm = PhysicalModel::create<PVModuleBoundaryModel>(type, boundary, options);
  }

  return(pm);
}



void
PVModuleBoundaryModel::do_init(void)
{
  string type = get_option("type", "ground");
  if (type == "ground") _contact_type = GND;
  if (type == "source") _contact_type = SRC;

  string layer = get_option("layer", "bottom");
  if (layer == "bottom") _contact_layer = BOTTOM;
  if (layer == "top") _contact_layer = TOP;
  if (layer == "both") _contact_layer = BOTH;
}

