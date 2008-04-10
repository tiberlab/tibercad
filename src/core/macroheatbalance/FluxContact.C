 #include "FluxContact.h"
 #include "SimulationOptions.h"



//===================================================================================//
FluxContact::FluxContact()
{
  
  set_type(ThermalContact::FluxContact);

}



 //===================================================================================//
void FluxContact::do_init()
{
 
  double flux = get_options().get_option("power_density",0.0);

  set_heat_flux(flux);
  

}
