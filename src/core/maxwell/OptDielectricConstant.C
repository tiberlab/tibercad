#include "OptDielectricConstant.h"
#include "Material.h"
#include "RotatedCrystal.h" 
#include "tensor.h"



//-------------------------------------------------------------------------//

void OptDielectricConstant::copy_from(const PhysicalModelInterface *rhs)
{
  const OptDielectricConstant* mod = dynamic_cast<const OptDielectricConstant*> (rhs);

 

  _dielectric_constant_real  = mod->_dielectric_constant_real ;

  _dielectric_constant_imag  = mod->_dielectric_constant_imag ;
}

//-------------------------------------------------------------------------//


void OptDielectricConstant::calculate_VCA (const PhysicalModelInterface *comp_A,
                                           const PhysicalModelInterface *comp_B, double xa)

{ 
  const OptDielectricConstant* modA = dynamic_cast<const OptDielectricConstant*>(comp_A);

  const OptDielectricConstant* modB = dynamic_cast<const OptDielectricConstant*>(comp_B);

 
 

  alloy(_dielectric_constant_real,modA->_dielectric_constant_real, modB->_dielectric_constant_real, xa);
  alloy(_dielectric_constant_imag,modA->_dielectric_constant_imag, modB->_dielectric_constant_imag, xa);

  Material* mat = get_material();
 
  const RotatedCrystal&   cr = mat->get_rotated_crystal();
 
  rotate_to_calculation_system(cr.RotMatrix);
}
