#include "LatticeThermalConductivity.h"
#include "ZbLatticeThermalConductivity.h"
#include "Material.h"
#include "RotatedCrystal.h" 
#include "tensor.h"
using namespace std;


//-------------------------------------------------------------------------//

void LatticeThermalConductivity::copy_from(const PhysicalModelInterface *rhs)
{
  const LatticeThermalConductivity* mod = dynamic_cast<const LatticeThermalConductivity*> (rhs);

  _conductivity = mod->_conductivity;
}

//-------------------------------------------------------------------------//


void LatticeThermalConductivity::calculate_VCA (const PhysicalModelInterface *comp_A, 
                                                const PhysicalModelInterface *comp_B, double xa) 
{ 
  const ZbLatticeThermalConductivity* modA = dynamic_cast<const ZbLatticeThermalConductivity*>(comp_A);

  const ZbLatticeThermalConductivity* modB = dynamic_cast<const ZbLatticeThermalConductivity*>(comp_B);


   alloy(_conductivity,modA->_conductivity, modB->_conductivity, xa);  

   Material* mat = get_material();
 
   const RotatedCrystal&   cr = mat->get_rotated_crystal();
 
   rotate_to_calculation_system(cr.RotMatrix);
}

//-------------------------------------------------------------------------//






//-------------------------------------------------------------------------//
