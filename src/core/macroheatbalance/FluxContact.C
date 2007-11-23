 #include "FluxContact.h"
 #include "SimulationOptions.h"



//===================================================================================//
FluxContact::FluxContact()
{
  
  set_type(ThermalContact::Neumann);

}



 //===================================================================================//
void FluxContact::do_init()
{
 
  double rho_e = get_options().get_option("electrons_resistivity",0.0);

  double rho_h = get_options().get_option("holes_resistivity",0.0);

  set_electrons_resistivity(rho_e);

  set_holes_resistivity(rho_h);
  

}
