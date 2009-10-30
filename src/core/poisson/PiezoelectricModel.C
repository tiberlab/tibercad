// $Id$


#include "PiezoelectricModel.h" 
#include "RotatedCrystal.h"
#include "Material.h"
#include "elem.h"



//----------------------------------------------------------------------//
void PiezoelectricModel:: do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
 //  const PiezoelectricModel* modA = dynamic_cast<const PiezoelectricModel*>(comp_A);

//   const PiezoelectricModel* modB = dynamic_cast<const PiezoelectricModel*>(comp_B);
  
//   alloy(_P,  modA->_P, modB->_P, xa);
 
//   Material*   mat = get_material();

//   const RotatedCrystal&   cr = mat->get_rotated_crystal ();

//   rotate_to_calc_system(cr.RotMatrix);
  

}
//--------------------------------------------------------------------//

//--------------------------------------------------------------------//
PiezoelectricModel::PiezoelectricModel() : PhysicalModelInterface( )
{
  _pol = 0;
  _strain = 0;
}
 
//--------------------------------------------------------------------//

void PiezoelectricModel::rotate_to_calc_system(const Tensor2Gen& RotMatrix)
{
 
  _pol = RotMatrix * _pol;

}
