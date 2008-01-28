// $Id$

#include "RelaxationMethod.h"
#include "Control.h"

using namespace std;



void
RelaxationMethod::parse_options(void)
{
  SelfconsistentSolver::parse_options();

  ModelOptions& opts = get_options();

  _relax = opts.get_option("relaxation_factor", _relax);
}




void
RelaxationMethod::do_solve(void)
{
  parse_options();

  int num_sim = get_number_of_simulations();

  SimulationIterator it(simulations_begin());
  const SimulationIterator end(simulations_end());

  initialize();

  ID old_sol_id = get_last_simulation()->remember_current_solution();

  for (unsigned int it = 0; it < get_maximum_iterations(); it++)
  {

    solve_simulations();

    double norm = get_last_simulation()->get_maximum_norm_of_difference(old_sol_id);

    bool converged = true;

    // check for the difference between old and new solutions
    cerr << "iteration " << it << ": ";

    if (norm > get_absolute_tolerance())
      converged = false;
    cerr << "norm = " << norm << endl;

    if (converged)
      break;

    get_last_simulation()->scale_solution(_relax);
    get_last_simulation()->add_scaled_remembered_solution(old_sol_id, 1.0 - _relax);
    get_last_simulation()->remember_current_solution(old_sol_id);
    
  }

  // clean up
  get_last_simulation()->delete_remembered_solution(old_sol_id);
}



