 #include "ThermalSurfaceConductance.h"
 #include "SimulationOptions.h"



//===================================================================================//
ThermalSurfaceConductance::ThermalSurfaceConductance()
{
  
  set_type(ThermalContact::ThermalSurfaceConductance);

}



 //===================================================================================//
void ThermalSurfaceConductance::do_init()
{
 
  double temp = get_options().get_option("temperature",SimulationOptions::temperature);

  double g_surf = get_options().get_option("g_surf",0.0);

  set_temperature(temp);

  set_thermal_surface_conductance(g_surf);
  

}
