 #include "FourierBTE.h"
 #include "SimulationOptions.h"



//===================================================================================//
FourierBTE::FourierBTE()
{

  set_type(ThermalContact::FourierBTE);

}


 
double
FourierBTE::get_temperature(const Elem* elem,const Point& p)
{
 
  double temperature = 0.0;

  temperature = _lattice_temp.get_temperature(elem,p);
 
  if (~_lattice_temp.is_solved() & _global_simulation)
    temperature = _global_lattice_temp.get_temperature(elem,p);
  
 

  return temperature;

}

 //===================================================================================//
void FourierBTE::do_init()
{

  std::string temp_simul = get_option("simulation_name", "no_sim");
  bool has_sim = _lattice_temp.set_simulation(temp_simul);
  _lattice_temp.set_simulation(temp_simul);


  //First Guess simulation
  _global_simulation = get_option("global_simulation",false);
  if  (_global_simulation)
  {
    std::string gsn = get_option("global_simulation_name", "no_sim");
    bool has_sim = _global_lattice_temp.set_simulation(gsn);
  } 



}
