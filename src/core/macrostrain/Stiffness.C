#include "Stiffness.h"

//--------------------------------------------------------------------//

Stiffness::Stiffness() : PhysicalProperties("Stiffness")
{
  C_cr = 0;
}
 
//--------------------------------------------------------------------//

void Stiffness::rotate_to_calc_system(const Tensor2Gen& RotMatrix)
{
  // generates stiffness matrix in calculation system
 
  C_calc = push_forward(C_cr, RotMatrix);
}

//--------------------------------------------------------------------//
 
void Stiffness::set_C_tensor_crystal(const Tensor4DSym&     C)
{
  C_cr = C;
}

