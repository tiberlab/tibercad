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

  get_parameter("temperature", _temperature);

}
