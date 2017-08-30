// $Id: DDHeatSource.C 2418 2011-02-28 03:40:59Z gromano $

#include "DDHeatSource.h"
#include "Material.h"

#include "TiberModule.h"





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

  vector<string> heat_sources;
  get_options().get_option("heat_sources", heat_sources);

  if (heat_sources.empty())
  {
    ID_set.insert(_simul->get_solution_id("eJoule"));
    ID_set.insert(_simul->get_solution_id("hJoule"));
    ID_set.insert(_simul->get_solution_id("RecombHeat"));
    ID_set.insert(_simul->get_solution_id("ePeltier"));
    ID_set.insert(_simul->get_solution_id("hPeltier"));
  }
  else
  {
    for (auto& src : heat_sources)
    {
      ID_set.insert(_simul->get_solution_id(src));
    }
  }
}

void
DDHeatSource::calculate(const Elem* elem, const Point& point)
{

  std::map<ID, vector< double > > solution;
  for (auto& id : ID_set)
    solution[id].resize(0);

  double heat_source = 0.0;
  vector<Point> h_point(1);
  h_point[0] = point;

  if (_simul->get_solution(elem, solution, h_point))
  {
    for (auto& src : solution)
      heat_source += src.second[0];
  }
  set_heat_source(heat_source * 1E6);

}


 
