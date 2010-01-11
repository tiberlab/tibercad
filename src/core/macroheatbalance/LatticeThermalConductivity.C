// $Id$

#include "LatticeThermalConductivity.h"
#include "Material.h"
#include "SimulationOptions.h"
#include "RotatedCrystal.h"
#include "tensor.h"


LatticeThermalConductivity::LatticeThermalConductivity(const ModelOptions& options) :
  PhysicalModelInterface(options),
  _conductivity(0)
{
 _temperature =  SimulationOptions::temperature;

}
//-------------------------------------------------------------------------//



