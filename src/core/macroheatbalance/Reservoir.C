 #include "Reservoir.h"
 #include "SimulationOptions.h"



//===================================================================================//
Reservoir::Reservoir()
{
  
  set_type(ThermalContact::Reservoir);

}



 //===================================================================================//
void Reservoir::do_init()
{
 
  std::string s(get_options().get_option("temperature", ""));

  //Insert the default written in the input file or the temperature of th simulation
  set_temperature(check_and_register(s, SimulationOptions::temperature));


}
