// $Id$

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


