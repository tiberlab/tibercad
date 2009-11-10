#include "ThermalContact.h"

#include "FluxContact.h" 
#include "Reservoir.h" 
#include "ThermalSurfaceResistance.h"
#include "ThermalSurfaceConductance.h"

//==================================================================================//
ThermalContact*
ThermalContact::create(const std::string & name,  const ModelOptions &   options)
{
  ThermalContact* result = NULL;

  if (name == "heat_reservoir")
    result = Reservoir::create();

  if (name == "thermal_flux")
    result = FluxContact::create();

  if (name == "thermal_surface_resistance")
    result = ThermalSurfaceResistance::create();

  if (name == "thermal_surface_conductance")
    result = ThermalSurfaceConductance::create();
 

  if (result != NULL)
  {
     result->set_options(options);
  }


  return result;
}




