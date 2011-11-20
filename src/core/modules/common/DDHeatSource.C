// $Id: DDHeatSource.C 2418 2011-02-28 03:40:59Z gromano $

#include "DDHeatSource.h"
#include "Material.h"



// The first string is the class name, the second one
// is the type of the model (here it is a bulk model),
// the third one is the specific model implementation.
// The library name will then be bulk_default.so

TIBER_MODULE(DDHeatSource, heat_source, joule)

using namespace std;


DDHeatSource::DDHeatSource(const ModelOptions& options):HeatSourceModel(options)
{

}

void
DDHeatSource::do_init(void)
{

  
  string dd_simul_name = get_options().get_option("transport_simulation", "driftdiffusion");
  _simul = SimulationInterface::find_simulation(dd_simul_name);

  if ( _simul == NULL)
   throw InitFailedException("Could not find " + dd_simul_name);

  var_map.clear();
  var_map[EJOULE]=_simul->get_solution_id("eJoule");
  var_map[HJOULE]=_simul->get_solution_id("hJoule");
  var_map[RECHEAT]=_simul->get_solution_id("RecombHeat");
  var_map[EPELTH]=_simul->get_solution_id("ePeltier");
  var_map[HPELTH]=_simul->get_solution_id("hPeltier");
  var_map[WNX]=_simul->get_solution_id("ePowerFlux");
  var_map[WPX]=_simul->get_solution_id("hPowerFlux");

  std::map<ID,ID>::iterator      it(var_map.begin());
  std::map<ID,ID>::iterator      end(var_map.end());
  for(; it!=end; ++it)
    ID_set.insert(it->second);
  

}

void
DDHeatSource::calculate(const Elem* elem, const Point& point)
{

  std::map<ID, vector< double > > solution;
  solution[var_map[EJOULE]].resize(0);
  solution[var_map[HJOULE]].resize(0);
  solution[var_map[RECHEAT]].resize(0);
  solution[var_map[EPELTH]].resize(0);
  solution[var_map[HPELTH]].resize(0);

  double heat_source = 0.0;
  vector<Point> h_point(1);
  h_point[0] = point;

  if  (_simul->get_solution(elem, solution, h_point))
  {
    heat_source += solution[var_map[EJOULE]][0];
    heat_source += solution[var_map[HJOULE]][0];
    heat_source += solution[var_map[RECHEAT]][0];
    heat_source += solution[var_map[EPELTH]][0];
    heat_source += solution[var_map[HPELTH]][0];
  }
  set_heat_source(heat_source * 10E6);   

}


 
