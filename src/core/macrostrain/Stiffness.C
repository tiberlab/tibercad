#include "Stiffness.h" 
#include "RotatedCrystal.h"
#include "Material.h"
//----------------------------------------------------------------------//
void Stiffness::copy_from (const PhysicalModelInterface *rhs)
{

  const Stiffness* mod = dynamic_cast<const Stiffness*>(rhs);

  C_cr = mod->C_cr;

  

} 
 

//----------------------------------------------------------------------//
void Stiffness:: calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
  const Stiffness* modA = dynamic_cast<const Stiffness*>(comp_A);

  const Stiffness* modB = dynamic_cast<const Stiffness*>(comp_B);
  
  alloy(C_cr,  modA->C_cr, modB->C_cr, xa);
 
  Material*   mat = get_material();

  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

  rotate_to_calc_system(cr.RotMatrix);
  

}
//--------------------------------------------------------------------//

//--------------------------------------------------------------------//

Stiffness::Stiffness() : PhysicalModelInterface( )
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

//--------------------------------------------------------------------//
