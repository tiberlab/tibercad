// $Id$

#include "ExcitonDissociation.h"

#include "SimulationInterface.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"


TIBER_MODULE(ExcitonDissociation, recombination, exciton_dissociation)



void
ExcitonDissociation::do_init(void)
{
  d_ = get_option("damping", 1.0);

  std::string ex = get_option("exciton_simulation", "");

  // find the exciton simulation to use
  _exciton_sim = SimulationInterface::find_simulation(ex);

  if (_exciton_sim == NULL)
  {
    std::string msg("ExcitonDissociation: Simulation " +
        std::string(ex) + " not found");
    throw InitFailedException(msg);
  }

  _Rdiss_id = _exciton_sim->get_solution_id("dissociation");
}




void
ExcitonDissociation::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{

  DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  const Elem* el = dd.get_element();

  recomb_e = 0.0;

  // we only use the exciton simulation if it has been solved before
  if (_exciton_sim->is_solved())
  {
    double x = 0.0;
    bool succ = _exciton_sim->get_solution(el, dd.get_coordinates(), _Rdiss_id, x);
    if (succ)
      recomb_e = -d_ * x;
  }

  recomb_h = recomb_e;
}




void
ExcitonDissociation::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  recomb_e[0] = recomb_h[0] = 0.0;
  recomb_e[1] = recomb_h[1] = 0.0;
}





