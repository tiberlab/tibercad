 #include "Reservoir.h"


//===================================================================================//
Reservoir::Reservoir()
{
  
  set_type(ThermalContact::Reservoir);
}



 //===================================================================================//
void Reservoir::do_init()
{
  const ModelOptions& options =	get_options ();

  _temperature=options.get_option("Temperature",300.0);
  

}
