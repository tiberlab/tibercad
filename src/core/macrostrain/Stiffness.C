// $Id$


#include "Stiffness.h"
#include "ZbStiffness.h"
#include "WzStiffness.h"
#include "RotatedCrystal.h"
#include "Material.h"


//--------------------------------------------------------------------//

Stiffness::Stiffness(const ModelOptions& options) : PhysicalModel(options)
{
  C_cr = 0;



}



Stiffness* Stiffness::create(const Material* mat, const ModelOptions& options)
{
  std::string str = mat->get_structure();
  Stiffness* mod = NULL;
  if (str == "zb")
    mod = ZbStiffness::create(options);
  else if (str == "wz")
    mod = WzStiffness::create(options);

  return mod;
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


