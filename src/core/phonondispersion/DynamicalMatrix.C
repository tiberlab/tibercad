// $Id$

#include "DynamicalMatrix.h"
#include "Material.h"
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


