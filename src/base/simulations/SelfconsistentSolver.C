// $Id$

#include "SelfconsistentSolver.h"
#include "Control.h"

using namespace std;



void
SelfconsistentSolver::do_init(void)
{
  ModelOptions& opts = get_options();
  Control& control = get_control();
  
  // get the names of the simulations to be solved
  vector<string> sims;
  opts.get_option("simulations", sims);
  int num_of_sims = sims.size();

  if (num_of_sims == 0)
    throw InitFailedException("Sweep: No simulation names provided.");
  
  _simulations.resize(num_of_sims);
  for (int i = 0; i < num_of_sims; i++)
  {
    _simulations[i] = control.find_simulation(sims[i]);
    if (_simulations[i] == NULL)
      throw InitFailedException("Sweep: Simulation " + sims[i] + " not found.");

    if (!_simulations[i]->is_initialized())
      _simulations[i]->init();
  }

  
  // we set our environment to that of the first simulation
  set_environment(&_simulations[0]->get_environment());
}


void
SelfconsistentSolver::parse_options(void)
{
  ModelOptions& opts = get_options();

  _max_it = opts.get_option("max_iterations", _max_it);
  _rel_tol = opts.get_option("rel_tolerance", _rel_tol);
  _abs_tol = opts.get_option("abs_tolerance", _abs_tol);
  _relax = opts.get_option("relaxation_factor", _relax);

  int num_sim = _simulations.size();
  for (int i = 0; i < num_sim; i++)
    _simulations[i]->set_relaxation_factor(_relax);
}


void
SelfconsistentSolver::do_equilibrium(void)
{
  int num_sim = _simulations.size();
  
  for (int i = 0; i < num_sim; i++)
  {
    // some sanity check
    assert(_simulations[i] != NULL);
    assert(_simulations[i]->is_initialized());

    _simulations[i]->solve_equilibrium();
  }
}




void
SelfconsistentSolver::do_solve(void)
{
  int num_sim = _simulations.size();

  assert(num_sim > 0);

  for (int i = 0; i < num_sim; i++)
    _simulations[i]->solve();

  // we make a copy of the current solutions
  vector<ID> old_sol_ids(num_sim); 
  for (int i = 0; i < num_sim; i++)
    old_sol_ids[i] = _simulations[i]->remember_current_solution();


  // for the norms of the differences
  vector<double> norms(num_sim);

  for (unsigned int it = 0; it < _max_it; it++)
  {

    for (int i = 0; i < num_sim; i++)
      _simulations[i]->solve();

    bool converged = true;

    // check for the difference between old and new solutions
    cerr << "iteration " << it << ": ";
    for (int i = 0; i < num_sim; i++)
    {
      norms[i] =
        get_norm_of_difference(_simulations[i]->get_solution_vector(),
            *(_simulations[i]->get_remembered_solution(old_sol_ids[i])));

      if (norms[i] > _abs_tol)
        converged = false;
      
      cerr << "norm[" << i << "] = " << norms[i] << " ";
    }
    cerr << endl;

    if (converged)
      break;
    
    for (int i = 0; i < num_sim; i++)
      _simulations[i]->remember_current_solution(old_sol_ids[i]);
  }

  // clean up
  for (int i = 0; i < num_sim; i++)
    _simulations[i]->delete_remembered_solution(old_sol_ids[i]);
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




ID
SelfconsistentSolver::do_remember_current_solution(ID id)
{
  int num_sim = _simulations.size();

  map<ID, vector<ID> >::iterator end(_remembered_sol_ids.end());
  map<ID, vector<ID> >::iterator it(_remembered_sol_ids.find(id));

  if (it != end)
  {
    assert((it->second).size() == num_sim);
    for (int i = 0; i < num_sim; i++)
      (it->second)[i] = _simulations[i]->remember_current_solution((it->second)[i]);
  }
  else
  {
    if (_remembered_sol_ids.begin() == end)
      id = 1;
    else
      id = (--end)->first + 1;

    vector<ID> ids(num_sim);
    for (int i = 0; i < num_sim; i++)
      ids[i] = _simulations[i]->remember_current_solution();

    _remembered_sol_ids[id] = ids;
  }

  return id;
}


void
SelfconsistentSolver::do_set_to_remembered_solution(ID id)
{
  int num_sim = _simulations.size();

  map<ID, vector<ID> >::iterator end(_remembered_sol_ids.end());
  map<ID, vector<ID> >::iterator it(_remembered_sol_ids.find(id));

  if (it != end)
    for (int i = 0; i < num_sim; i++)
      _simulations[i]->set_to_remembered_solution((it->second)[i]);
}



void
SelfconsistentSolver::do_delete_remembered_solution(ID id)
{
  int num_sim = _simulations.size();

  map<ID, vector<ID> >::iterator end(_remembered_sol_ids.end());
  map<ID, vector<ID> >::iterator it(_remembered_sol_ids.find(id));
  if (it != end)
    for (int i = 0; i < num_sim; i++)
      _simulations[i]->delete_remembered_solution((it->second)[i]);
}



void
SelfconsistentSolver::do_plot(void)
{
  int num_sim = _simulations.size();

  for (int i = 0; i < num_sim; i++)
    _simulations[i]->plot();
}
