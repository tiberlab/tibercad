// $Id$

#include "MesoThermalContact.h"
#include "SimulationOptions.h"



//===================================================================================//
MesoThermalContact::MesoThermalContact(const ModelOptions& options)
  : ThermalContact(options)
{

  _temperature = SimulationOptions::temperature;


  _absorbivity = 1.0;
  _temperature = SimulationOptions::temperature;
  _reflectivity = 0.0;
  _diffusivity = 0.0;

  set_type(ThermalContact::Meso);

}



 //===================================================================================//
void MesoThermalContact::do_init()
{

 //  get_parameter("absorbivity", _absorbivity);
//   get_parameter("temperature", _temperature);
//   get_parameter("reflectivity", _reflectivity);
//   get_parameter("diffusivity", _diffusivity);
 get_parameter("A", _absorbivity);
  get_parameter("T", _temperature);
  get_parameter("R", _reflectivity);
  get_parameter("D", _diffusivity);
}
