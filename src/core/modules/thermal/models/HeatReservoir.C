 // $Id: HeatReservoir.C 2362 2011-02-21 01:02:31Z gromano $

#include "HeatReservoir.h"

#include "TiberModule.h"
#include "SimulationInterface.h"


void
HeatReservoir::do_init(void)
{
  get_parameter("temperature", _temperature);

  std::string sim_name = get_option("host_simulation", "");

  _host_sim = 0;
  if (sim_name != "") {
    SimulationInterface* sim = SimulationInterface::find_simulation(sim_name);

    if (sim == nullptr) {
      std::string msg("External thermal host simulation " + std::string(sim_name) + " not found");
      throw InitFailedException(msg);
    }
    _host_sim = sim->get_id();
  }
}

void
HeatReservoir::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  if (_host_sim > 0)
  {
    SimulationInterface* sim = SimulationInterface::get_simulation(_host_sim);
    if (sim->is_solved())
    {
      ID T_id = sim->get_solution_id("LatticeTemp");

      double T = 0.0;
      sim->get_solution(elem, T_id, T, point);
      _temperature = T;
    }
  }

  set_coefficients(1.0, 0.0, _temperature);
}

