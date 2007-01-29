#include "SelfconsistentSolver.h"


void
SelfconsistentSolver::do_init(void)
{
  ModelOptions& opts = get_options();
  
  const std::string& sim1 = opts.get_option("simulation1", "");
  _simulation1 = SimulationInterface::find_simulation(sim1);
  
  const std::string& sim2 = opts.get_option("simulation2", "");
  _simulation2 = SimulationInterface::find_simulation(sim2);

  // we don't tolerate NULL pointers...
  if (_simulation1 == NULL)
    throw InitFailedException("Simulation " + sim1 + " not found");
  
  if (_simulation2 == NULL)
    throw InitFailedException("Simulation " + sim2 + " not found");
}


void
SelfconsistentSolver::parse_options(void)
{
  ModelOptions& opts = get_options();

  _max_it = opts.get_option("max_iterations", _max_it);
  _rel_tol = opts.get_option("rel_tolerance", _rel_tol);
  _abs_tol = opts.get_option("abs_tolerance", _abs_tol);
  _relax = opts.get_option("relaxation_factor", _relax);

  _simulation1->set_relaxation_factor(_relax);
  _simulation2->set_relaxation_factor(_relax);
}


void
SelfconsistentSolver::do_solve(void)
{
  assert(_simulation1 != NULL);
  assert(_simulation2 != NULL);

  _simulation1->solve();
  _simulation2->solve();

  for (unsigned int i = 0; i < _max_it; i++)
  {
    _simulation1->solve();
    _simulation2->solve();

    // check for the difference between old and new solutions
  }
}


