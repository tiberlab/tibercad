 #include "ThermalSurfaceResistance.h"
 #include "SimulationOptions.h"



//===================================================================================//
ThermalSurfaceResistance::ThermalSurfaceResistance(const ModelOptions& options)
 : ThermalContact(options)
{
  
  set_type(ThermalContact::ThermalSurfaceResistance);

}



 //===================================================================================//
void ThermalSurfaceResistance::do_init()
{
 
  double temp = get_options().get_option("temperature",SimulationOptions::temperature);

  double r_surf = get_options().get_option("r_surf",0.0);

  set_temperature(temp);

  set_thermal_surface_resistance(r_surf);
  

}
