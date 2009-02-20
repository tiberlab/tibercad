// $Id$

#include "SelfconsistentSolver.h"
#include "Control.h"
#include "XMonitor.h"

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
    throw InitFailedException(
        "SelfconsistentSolver: No simulation names provided.");
  
  _simulations.resize(num_of_sims);
  for (int i = 0; i < num_of_sims; i++)
  {
    _simulations[i] = control.find_simulation(sims[i]);
    if (_simulations[i] == NULL)
      throw InitFailedException("SelfconsistentSolver: Simulation " +
          sims[i] + " not found.");

    // If it is not already initialized, we initialize now
    if (!_simulations[i]->is_initialized())
      _simulations[i]->init();
  }

  if (!_simulations[num_of_sims - 1]->has_solution_vector())
  {
    ostringstream s;
    s << "SelfconsistentSolver: Simulation "
      << _simulations[num_of_sims - 1]->get_name() << " has no solution vector!";
    throw InitFailedException(s.str());
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
  _monitor = opts.get_option("monitor", true);

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




NumericVector<double>&
SelfconsistentSolver::do_get_solution_vector(void)
{
  int num_sim = _simulations.size();
  return _simulations[num_sim - 1]->get_solution_vector();
}



void
SelfconsistentSolver::do_set_solution_vector(
    const NumericVector<double>& new_solution)
{
  int num_sim = _simulations.size();
  _simulations[num_sim - 1]->set_solution_vector(new_solution);
}



void
SelfconsistentSolver::initialize(void)
{
  solve_simulations();
  //get_last_simulation()->solve();
}



void
SelfconsistentSolver::solve_simulations(void)
{
  int num_sim = _simulations.size();
  for (int i = 0; i < num_sim; i++)
    _simulations[i]->solve();
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




void
SelfconsistentSolver::build_integrated_quantities(const set<string>& names,
    vector<double>& values)
{
  vector<double> val;

  int num_sim = _simulations.size();

  for (int i = 0; i < num_sim; i++)
  {
    _simulations[i]->get_integrated_quantities(names, val);
    values.insert(values.end(), val.begin(), val.end());
  }
}




void
SelfconsistentSolver::build_integrated_quantities_description(
    const set<string>& names,
    vector<string>& legend,
    vector<string>& description)
{
  vector<string> leg;
  vector<string> desc;

  int num_sim = _simulations.size();

  for (int i = 0; i < num_sim; i++)
  {
    _simulations[i]->get_integrated_quantities_description(names, leg, desc);
    legend.insert(legend.end(), leg.begin(), leg.end());
    description.insert(description.end(), desc.begin(), desc.end());
  }
}




void
SelfconsistentSolver::open_xmonitor(void)
{
  if (get_options().get_option("xmonitor", false))
  {
    _xmonitor = XMonitor::create(get_name() + ": convergence monitor");
    _xmonitor->set_axis_labels("iteration nr.", "Logarithm of relative error");
  }
}


void
SelfconsistentSolver::close_xmonitor(void)
{
  delete _xmonitor;
  _xmonitor = NULL;
}


void
SelfconsistentSolver::draw_point(double iteration, double error, bool logarithm)
{
  if (_xmonitor != NULL)
  {
    if (logarithm)
      error = log10(error);

    _xmonitor->draw_point(iteration, error);
  }
}
