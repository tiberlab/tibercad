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
 * \file RamanTensor.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */


#include "RamanTensor.h"
#include "tibercad/physics/Material.h"
#include "RotatedCrystal.h" 
#include "tensor.h"


RamanTensor::RamanTensor(const ModelOptions& options)
 : PhysicalModel(options)
{  
}


//-------------------------------------------------------------------------//


void RamanTensor::do_init_alloy (const PhysicalModel *comp_A, 
                                                const PhysicalModel *comp_B, double xa) 
{ 
  const RamanTensor* modA = dynamic_cast<const RamanTensor*>(comp_A);

  const RamanTensor* modB = dynamic_cast<const RamanTensor*>(comp_B);

  alloy(_raman_tensor[0],modA->_raman_tensor[0], modB->_raman_tensor[0], xa);  
  alloy(_raman_tensor[1],modA->_raman_tensor[1], modB->_raman_tensor[1], xa);
  alloy(_raman_tensor[2],modA->_raman_tensor[2], modB->_raman_tensor[2], xa);
 
  Material* mat = get_material();
  
  const RotatedCrystal&   cr = mat->get_rotated_crystal();
  
  rotate_to_calculation_system(cr.RotMatrix);
}


