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
    result = Reservoir::create(options);

  if (name == "thermal_flux")
    result = FluxContact::create(options);

  if (name == "FourierBTE")
    result = FourierBTE::create(options);

 if (name == "BTEFourier")
    result = BTEFourier::create(options);

  if (name == "thermal_surface_resistance")
    result = ThermalSurfaceResistance::create(options);

  if (name == "thermal_surface_conductance")
    result = ThermalSurfaceConductance::create(options);
 
  if (name == "specular")
    result = Specular::create(options);

  if (name == "diffusive")
    result = Diffusive::create(options);


  return result;
}




