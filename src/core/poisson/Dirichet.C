 #include "Dirichlet.h"
 #include "SimulationOptions.h"



//===================================================================================//
Dirichlet::Dirichlet()
{
  
  set_type( PoissonContact::Dirichlet);

}



 //===================================================================================//
void Dirichlet::do_init()
{
 
  std::string s(get_options().get_option("Voltage", ""));

  set_potential(check_and_register(s,0.0));


}
