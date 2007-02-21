#include "ThermalContact.h"


#include "Reservoir.h" 

//==================================================================================//
ThermalContact*
ThermalContact::create(const std::string & name,  const ModelOptions &   options)
{
  ThermalContact* result = NULL;

  if (name == "Heat_reservoir")
    result = Reservoir::create();

  if (result != NULL)
    result->set_options(options);

  return result;
}




