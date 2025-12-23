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
 * \file ConstantBodyForce.C
 * \brief tiberCAD elasticity module implementation.
 *
 * \note This file is part of module elasticity.
 */


#include "ConstantBodyForce.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"


using namespace std;


ConstantBodyForce::ConstantBodyForce(const ModelOptions& options):BodyForceModel(options)
{
}

void
ConstantBodyForce::do_init(void)
{
  libMesh::RealGradient force_source(0);
  get_parameter("F", force_source);
  set_force_source(force_source);

  libMesh::RealTensor dummy_tens(0);
  set_strain_source(dummy_tens);
  set_stress_source(dummy_tens);
}




 
