/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file DynamicalMatrix.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */


#include "DynamicalMatrix.h"
#include "tibercad/physics/Material.h"
#include "RotatedCrystal.h" 
#include "tensor.h"


DynamicalMatrix::DynamicalMatrix(const ModelOptions& options)
 : PhysicalModel(options),
   _dynamical_matrix(0)
{  
  
}



//-------------------------------------------------------------------------//


void DynamicalMatrix::do_init_alloy (const PhysicalModel *comp_A, 
                                                const PhysicalModel *comp_B, double xa) 
{ 
  const DynamicalMatrix* modA = dynamic_cast<const DynamicalMatrix*>(comp_A);

  const DynamicalMatrix* modB = dynamic_cast<const DynamicalMatrix*>(comp_B);

  alloy(_dynamical_matrix,modA->_dynamical_matrix, modB->_dynamical_matrix, xa);  

  Material* mat = get_material();
  
  const RotatedCrystal&   cr = mat->get_rotated_crystal();
  
  rotate_to_calculation_system(cr.RotMatrix);
}


