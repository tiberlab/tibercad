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
 
  std::string s(get_options().get_option("polarization", ""));

  set_polarization(check_and_register(s,0.0));


}
