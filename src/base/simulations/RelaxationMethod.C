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

  for ( ; it != end; ++it)
    (*it)->solve();


  // we make a copy of the current solutions
  vector<ID> old_sol_ids(num_sim); 
  for (int i = 0; i < num_sim; i++)
    old_sol_ids[i] = simulation(i)->remember_current_solution();


  // for the norms of the differences
  vector<double> norms(num_sim, 0.0);

  for (unsigned int it = 0; it < get_maximum_iterations(); it++)
  {

    for (int i = 0; i < num_sim; i++)
    {
      simulation(i)->solve();
      
      // we get the norm of the difference before doing relaxation
      double norm =
        simulation(i)->get_maximum_norm_of_difference(old_sol_ids[i]);

      //if ((norms[i] > _abs_tol) && ((norm / norms[i]) < 0.5))
      //  _relax = min(_relax * 1.5, 1.0);
      //else if ((norms[i] > _abs_tol) && ((norm / norms[i]) > 1.5))
      //  _relax = max(_relax / 2.0, 0.001);
      //cerr << "relaxation: " << _relax << endl;

      norms[i] = norm;

      simulation(i)->scale_solution(_relax);
      simulation(i)->add_scaled_remembered_solution(old_sol_ids[i], 1.0 - _relax);
    }

    bool converged = true;

    // check for the difference between old and new solutions
    cerr << "iteration " << it << ": ";
    for (int i = 0; i < num_sim; i++)
    {
      if (norms[i] > get_absolute_tolerance())
        converged = false;
      
      cerr << "norm[" << i << "] = " << norms[i] << " ";
    }
    cerr << endl;

    if (converged)
      break;
    
    for (int i = 0; i < num_sim; i++)
      simulation(i)->remember_current_solution(old_sol_ids[i]);
  }

  // clean up
  for (int i = 0; i < num_sim; i++)
    simulation(i)->delete_remembered_solution(old_sol_ids[i]);
}



