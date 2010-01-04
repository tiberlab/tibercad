#include "ThermalContact.h"

#include "FluxContact.h" 
#include "Reservoir.h" 
#include "FourierBTE.h" 
#include "ThermalSurfaceResistance.h"
#include "ThermalSurfaceConductance.h"
#include "BTEFourier.h"
#include "Specular.h"
#include "Diffusive.h"
//==================================================================================//
ThermalContact*
ThermalContact::create(const std::string & name,  const ModelOptions &   options)
{
  ThermalContact* result = NULL;

  if (name == "heat_reservoir")
    result = Reservoir::create();

  if (name == "thermal_flux")
    result = FluxContact::create();

  if (name == "FourierBTE")
    result = FourierBTE::create();

 if (name == "BTEFourier")
    result = BTEFourier::create();

  if (name == "thermal_surface_resistance")
    result = ThermalSurfaceResistance::create();

  if (name == "thermal_surface_conductance")
    result = ThermalSurfaceConductance::create();

  if (name == "specular")
    result = Specular::create();

  if (name == "diffusive")
    result = Diffusive::create();

  if (result != NULL)
  {
     result->set_options(options);
  }


  return result;
}




