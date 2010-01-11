// $Id$


#include "Stiffness.h"
#include "RotatedCrystal.h"
#include "Material.h"


//--------------------------------------------------------------------//

Stiffness::Stiffness(const ModelOptions& options) : PhysicalModelInterface(options)
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

inline void Stiffness::set_C_tensor_crystal(const Tensor4DSym&     C)
{
  C_cr = C;
}

//--------------------------------------------------------------------//


