#include "DynamicalMatrix.h"
#include "Material.h"
#include "RotatedCrystal.h" 
#include "tensor.h"

//class DynamicalMatrix;

DynamicalMatrix::DynamicalMatrix()
 :_dynamical_matrix(0)
{  
  
}
//-------------------------------------------------------------------------//

void DynamicalMatrix::copy_from(const PhysicalModelInterface *rhs)
{
  const DynamicalMatrix* mod = dynamic_cast<const DynamicalMatrix*> (rhs);

  _dynamical_matrix = mod->_dynamical_matrix;
}

//-------------------------------------------------------------------------//


void DynamicalMatrix::calculate_VCA (const PhysicalModelInterface *comp_A, 
                                                const PhysicalModelInterface *comp_B, double xa) 
{ 
  const DynamicalMatrix* modA = dynamic_cast<const DynamicalMatrix*>(comp_A);

  const DynamicalMatrix* modB = dynamic_cast<const DynamicalMatrix*>(comp_B);

  alloy(_dynamical_matrix,modA->_dynamical_matrix, modB->_dynamical_matrix, xa);  

  Material* mat = get_material();
  
  const RotatedCrystal&   cr = mat->get_rotated_crystal();
  
  rotate_to_calculation_system(cr.RotMatrix);
}


