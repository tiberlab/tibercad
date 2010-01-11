 #include "Reservoir.h"
 #include "SimulationOptions.h"



//===================================================================================//
Reservoir::Reservoir(const ModelOptions& options)
 : ThermalContact(options)
{

  set_type(ThermalContact::Reservoir);

}



 //===================================================================================//
void Reservoir::do_init()
{

  get_parameter("temperature", _temperature);

}
