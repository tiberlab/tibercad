 #include "Neumann.h"
 #include "SimulationOptions.h"



//===================================================================================//
Neumann::Neumann()
{

  set_type(PoissonContact::Neumann);

}



 //===================================================================================//
void Neumann::do_init()
{

  get_parameter("Polarization", _polarization);

}
