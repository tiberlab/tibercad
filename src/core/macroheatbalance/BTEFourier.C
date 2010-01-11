// $Id$

#include "BTEFourier.h"
#include "SimulationOptions.h"



//===================================================================================//
BTEFourier::BTEFourier(const ModelOptions& options)
  : ThermalContact(options)
{
  
  set_type(ThermalContact::BTEFourier);

}


RealGradient
BTEFourier::get_heat_flux(const Elem* elem, const Point& p) 
{

  RealGradient heat_flux(0);


  if (_simul->is_solved())
  {
    std::vector<Point> h_point(1);
    h_point[0] = p;
    std::vector< std::map< ID, double > > solution;
    _simul->get_solution(elem,h_point,ID_set,solution);
    
    
    heat_flux(0) = solution[0].find(var_map[JQX])->second;
    heat_flux(1) = solution[0].find(var_map[JQY])->second;
    heat_flux(2) = solution[0].find(var_map[JQZ])->second;

   
  }       


  return heat_flux;
}

 //===================================================================================//
void BTEFourier::do_init()
{
 
  
  std::string temp_simul = get_option("simulation_name", "no_sim");
 
  _simul = SimulationInterface::find_simulation(temp_simul);
  if ( _simul == NULL)
    throw InitFailedException("Could not find " + temp_simul);

  var_map[JQX]=_simul->get_variable_id("Jqx");
  var_map[JQY]=_simul->get_variable_id("Jqy");
  var_map[JQZ]=_simul->get_variable_id("Jqz");



  std::map<ID,ID>::iterator      it(var_map.begin());
  std::map<ID,ID>::iterator      end(var_map.end());
  for(; it!=end; ++it)
    ID_set.insert(it->second);



}

