 #include "Dirichlet.h"
 #include "SimulationOptions.h"



//===================================================================================//
Dirichlet::Dirichlet()
  : _potential(0.0)
{
  
  set_type( PoissonContact::Dirichlet);

}



 //===================================================================================//
void Dirichlet::do_init()
{
 
  get_parameter("Voltage", _potential);

}
