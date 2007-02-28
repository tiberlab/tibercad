// $Id$

#include "SelfconsistentSolver.h"

using namespace std;



void
SelfconsistentSolver::do_init(void)
{
  ModelOptions& opts = get_options();
  
  const string& sim1 = opts.get_option("simulation1", "");
  _simulation1 = SimulationInterface::find_simulation(sim1);
  
  const string& sim2 = opts.get_option("simulation2", "");
  _simulation2 = SimulationInterface::find_simulation(sim2);

  // we don't tolerate NULL pointers...
  if (_simulation1 == NULL)
    throw InitFailedException("Simulation " + sim1 + " not found");
  
  if (_simulation2 == NULL)
    throw InitFailedException("Simulation " + sim2 + " not found");
  
  // we set our environment to that of the first simulation
  set_environment(&_simulation1->get_environment());
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
SelfconsistentSolver::do_equilibrium(void)
{
  
  // some sanity check
  assert(_simulation1 != NULL);
  assert(_simulation1->is_initialized());
  assert(_simulation2 != NULL);
  assert(_simulation2->is_initialized());

  _simulation1->solve_equilibrium();
  _simulation2->solve_equilibrium();
}




void
SelfconsistentSolver::do_solve(void)
{
  assert(_simulation1 != NULL);
  assert(_simulation2 != NULL);

  _simulation1->solve();
  _simulation2->solve();

  // we make a copy of the current solutions
  ID old_sol1; 
  ID old_sol2; 
  old_sol1 = _simulation1->remember_current_solution();
  old_sol2 = _simulation2->remember_current_solution();


  for (unsigned int i = 0; i < _max_it; i++)
  {
    _simulation1->solve();
    _simulation2->solve();

    // check for the difference between old and new solutions
    //double norm1 = get_norm_of_difference(*old_sol1,
    //    _simulation1->get_solution_vector());
    //double norm2 = get_norm_of_difference(*old_sol2,
    //    _simulation2->get_solution_vector());
    //cerr << "iteration " << i << ": norm1 = " << norm1 <<
    //  "  norm2 = " << norm2 << endl;

    //if ((norm1 <= _abs_tol) && (norm2 <= _abs_tol))
    //  break;
    
    _simulation1->remember_current_solution(old_sol1);
    _simulation2->remember_current_solution(old_sol2);
  }

  _simulation1->plot();
  _simulation2->plot();

  // clean up
  _simulation1->delete_remembered_solution(old_sol1);
  _simulation2->delete_remembered_solution(old_sol2);
}




double
SelfconsistentSolver::get_norm_of_difference(NumericVector<double>& vec1,
    NumericVector<double>& vec2)
{
  assert(vec1.size() == vec2.size());
  double norm = -1;

  unsigned int n = vec1.size();
  for (unsigned int i = 0; i < n; i++)
  {
    double d = fabs(vec1(i) - vec2(i));
    norm = (d > norm) ? d : norm;
  }

  return norm;
}
