 #include "Specular.h"
 #include "SimulationOptions.h"



//===================================================================================//
Specular::Specular()
{

  set_type(ThermalContact::Specular);

}

 //===================================================================================//
void Specular::do_init()
{

  get_parameter("temperature", _temperature);

}
