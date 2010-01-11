// $Id$

#include "Specular.h"
#include "SimulationOptions.h"



//===================================================================================//
Specular::Specular(const ModelOptions& options)
  : ThermalContact(options)
{

  set_type(ThermalContact::Specular);

}

 //===================================================================================//
void Specular::do_init()
{

  get_parameter("temperature", _temperature);

}
