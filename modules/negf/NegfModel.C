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
 * \file NegfModel.C
 * \brief tiberCAD negf module implementation.
 *
 * \note This file is part of module negf.
 */

/*
 * Negfmodel.C
 *
 *  Created on: Jan 18, 2012
 *      Author: fpalomba
 */

#include "NegfModel.h"
#include "tibercad/physics/Material.h"

NegfModel::NegfModel(const ModelOptions& options)
   : PhysicalModel(options)
{
}

NegfModel*
NegfModel::create(const Material* mat, const ModelOptions& options)
{
  return PhysicalModel::create<NegfModel>(_create,_destroy, mat, options);
}

void
NegfModel::prepare_submodels(void)
{
    ModelOptions opts;
    opts.set_option("mass",1.0);
    create_submodels(_ham_models,"hamiltonian",opts);

}

