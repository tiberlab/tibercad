// $Id$

#include "Diffusive.h"
#include "SimulationOptions.h"



//===================================================================================//
Diffusive::Diffusive(const ModelOptions& options)
  : ThermalContact(options)
{

  _temperature = SimulationOptions::temperature;
  _emittivity = 1;
  set_type(ThermalContact::Diffusive);

}



 //===================================================================================//
void Diffusive::do_init()
{

  get_parameter("temperature", _temperature);
  get_parameter("emittivity", _emittivity);
  

}
