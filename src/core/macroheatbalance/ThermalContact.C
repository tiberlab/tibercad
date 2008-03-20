#include "ThermalContact.h"

#include "FluxContact.h" 
#include "Reservoir.h" 
#include "ThermalResistance.h"

//==================================================================================//
ThermalContact*
ThermalContact::create(const std::string & name,  const ModelOptions &   options)
{
  ThermalContact* result = NULL;

  if (name == "heat_reservoir")
    result = Reservoir::create();

  if (name == "Thermal_flux")
    result = FluxContact::create();

  if (name == "thermal_resistance")
    result = ThermalResistance::create();
 

  if (result != NULL)
  {
       
     result->set_options(options);
  }
  else
  {
     
  }

  return result;
}




