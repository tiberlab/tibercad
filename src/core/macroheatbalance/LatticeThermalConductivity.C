// $Id$

#include "LatticeThermalConductivity.h"
#include "Material.h"
#include "SimulationOptions.h"
#include "RotatedCrystal.h" 
#include "tensor.h"


LatticeThermalConductivity::LatticeThermalConductivity() :
  _conductivity(0)
{  
 _temperature =  SimulationOptions::temperature; 
  
}
//-------------------------------------------------------------------------//



void LatticeThermalConductivity::do_init_alloy (const PhysicalModelInterface *comp_A, 
                                                const PhysicalModelInterface *comp_B, double xa) 
{ 
  const LatticeThermalConductivity* modA = dynamic_cast<const LatticeThermalConductivity*>(comp_A);

  const LatticeThermalConductivity* modB = dynamic_cast<const LatticeThermalConductivity*>(comp_B);

  alloy(_conductivity,modA->_conductivity, modB->_conductivity, xa);  

  Material* mat = get_material();
  
  const RotatedCrystal&   cr = mat->get_rotated_crystal();
  
  rotate_to_calculation_system(cr.RotMatrix);
}

