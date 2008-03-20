 #include "ThermalResistance.h"
 #include "SimulationOptions.h"



//===================================================================================//
ThermalResistance::ThermalResistance()
{
  
  set_type(ThermalContact::ThermalResistance);

}



 //===================================================================================//
void ThermalResistance::do_init()
{
 
  double Text = get_options().get_option("ext_temperature",SimulationOptions::temperature);

  double Rth = get_options().get_option("thermal_resistance",0.0);

  set_external_temperature(Text);

  set_thermal_resistance(Rth);
  

}
