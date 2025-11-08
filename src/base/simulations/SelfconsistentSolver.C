// $Id$

#include "SelfconsistentSolver.h"
#include "base/common/Multiscale.h"
#include "XMonitor.h"

using namespace std;


SelfconsistentSolver::~SelfconsistentSolver(void)
{
  delete _multiscale;
}

void
SelfconsistentSolver::do_init(void)
{
  ModelOptions& opts = get_options();
  
  // get the names of the simulations to be solved
  vector<string> sims;
  opts.get_option("solve", sims);
  int num_of_sims = sims.size();

  if (num_of_sims == 0)
    throw InitFailedException(
        "SelfconsistentSolver: No simulation names provided.");
  
  _simulations.resize(num_of_sims);
  // as default, convergence check is done on the last simulation of the list
  _convergence_check_id = num_of_sims - 1;
  for (int i = 0; i < num_of_sims; i++)
  {
    _simulations[i] = find_simulation(sims[i]);
    if (_simulations[i] == NULL)
      throw InitFailedException("SelfconsistentSolver: Simulation " +
          sims[i] + " not found.");
    
    if (sims[i] == get_option("convergence_check", ""))
    {
      _convergence_check_id = i;
      cout<<"Convergence check: "<<sims[i]<<endl; 
    }

    // If it is not already initialized, we initialize now
    if (!_simulations[i]->is_initialized())
      _simulations[i]->init();
  }


  if (!_simulations[_convergence_check_id]->has_solution_vector())
  {
    ostringstream s;
    s << "SelfconsistentSolver: Simulation "
      << _simulations[_convergence_check_id]->get_name() << " has no solution vector!";
    throw InitFailedException(s.str());
  }

  // we set our environment to that of the first simulation
  //set_environment(&_simulations[0]->get_environment());

  parse_options();

  if (get_options().has_submodel("Multiscale"))
  {
    ModelOptions::const_submodel_iterator it(get_options().submodels_begin("Multiscale"));
    _multiscale = new Multiscale(it->second);
  }
}




void
SelfconsistentSolver::do_print_info(void)
{
}


void
SelfconsistentSolver::parse_options(void)
{
  ModelOptions& opts = get_options();

  _max_it = opts.get_option("max_iterations", _max_it);
  _rel_tol = opts.get_option("relative_tolerance", _rel_tol);
  _abs_tol = opts.get_option("absolute_tolerance", _abs_tol);
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




libMesh::NumericVector<double>&
SelfconsistentSolver::do_get_solution_vector(void)
{
  int num_sim = _simulations.size();
  //return _simulations[num_sim - 1]->get_solution_vector();
  return _simulations[_convergence_check_id]->get_solution_vector();
}



void
SelfconsistentSolver::do_set_solution_vector(
    const libMesh::NumericVector<double>& new_solution)
{
  int num_sim = _simulations.size();
  //_simulations[num_sim - 1]->set_solution_vector(new_solution); 
  _simulations[_convergence_check_id]->set_solution_vector(new_solution);
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
  {
    if (_multiscale != NULL)
      _multiscale->reinit(_simulations[i]);

    _simulations[i]->solve();

    // inherit the solution IDs
    // we do this here because the number of components may depend on the
    // solutions
    const IDSet& plotvars = _simulations[i]->get_plotvariable_ids();
    for (IDSet::iterator it(plotvars.begin()); it != plotvars.end(); ++it)
    {
      ID id = *it;
      const SolutionDescriptor& descr = _simulations[i]->get_solution_descriptor(id);

      // adjust ID
      id = TB_MAX_SIM * id + _simulations[i]->get_id();
      declare_solution_ext(descr.name(), id, descr.type(), descr.location(),
          descr.units(), descr.n_components());
    }
  }
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
  TiberCad::prepend_to_filename_suffix(get_name());

  int num_sim = _simulations.size();

  for (int i = 0; i < num_sim; i++)
    _simulations[i]->plot();

  TiberCad::drop_first_filename_suffix();
}



void
SelfconsistentSolver::get_solution_secure(map<ID, vector<double> >& values)
{
  //map<ID, vector<double> > orig(values);
  // we return everything we find
  values.clear();

  int num_sim = _simulations.size();
  for (int i = 0; i < num_sim; i++)
  {
    map<ID, vector<double> > tmp;
    _simulations[i]->get_solution(tmp);
    ID simid = _simulations[i]->get_id();

    map<ID, vector<double> >::iterator it = tmp.begin();
    const map<ID, vector<double> >::iterator end = tmp.end();
    for ( ; it != end; ++it)
      values[TB_MAX_SIM * it->first + simid] = it->second;
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
