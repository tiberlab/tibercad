// $Id$

#include "Sweep.h"
#include "Ramp.h"
#include "Utils.h"
#include "Variable.h"
#include "Boundary.h"
#include "BoundaryProperties.h"
#include "SimulationEnvironment.h"
#include "Control.h"
#include "Messages.h"

#include <fstream>



using namespace std;



Sweep::~Sweep(void)
{
}




void
Sweep::do_init(void)
{

  ModelOptions& opts = get_options();
  Control& control = get_control();


  // get the names of the simulations to be solved
  vector<string> sims;
  opts.get_option("simulation", sims);
  int num_of_sims = sims.size();


  _simulations.resize(num_of_sims);
  for (int i = 0; i < num_of_sims; i++)
  {
    _simulations[i] = control.find_simulation(sims[i]);
    if (_simulations[i] == NULL)
      throw InitFailedException("Sweep: Simulation " + sims[i] + " not found.");

    // If it is not already initialized, we initialize now
    if (!_simulations[i]->is_initialized())
      _simulations[i]->init();
  }

  // if user didn't provide a simulation name, we take the first available
  if (num_of_sims == 0)
  {
    _simulations.resize(1);
    _simulations[0] = control.find_simulation("");
    if (_simulations[0] == NULL)
      throw InitFailedException("Sweep: No simulation found.");
  }

  //
  // at this point we have for sure one simulation
  //

  // we set our environment to that of the first simulation
  set_environment(&_simulations[0]->get_environment());


  // Now we have to find the model to the variable
  _variable = opts.get_option("variable", "");
  if (_variable == "")
  {
    string msg("Sweep: You have to provide the name of ");
    msg += "the sweep variable.";
    throw InitFailedException(msg);
  }
}





void
Sweep::parse_options(void)
{
  ModelOptions& opts = get_options();

  _min_step = opts.get_option("min_step", _min_step);
  _max_step = opts.get_option("max_step", _max_step);

  double start = opts.get_option("start", 0.0);
  double stop = opts.get_option("stop", 0.0);
  unsigned int steps = opts.get_option("steps", 1);
  steps = (steps < 1) ? 1 : steps;

  double step = (stop - start) / steps;

  if (Utils::almost_equal::compare(start, stop))
  {
    steps = 0;
    _values.resize(1);
  }
  else
    _values.resize(steps + 1);

  for (unsigned int i = 0; i <= steps; i++)
    _values[i] = start + i * step;

  // we can also specify a vector with the values
  opts.get_option("values", _values);


  // check if our variable exists
  if (!Variable::is_variable(_variable))
  {
    ostringstream o;
    o << "Sweep: The variable " << _variable << " is not defined.";
    throw InitFailedException(o.str());
  }


  // whether to plot data or not
  _plot_data = opts.get_option("plot_data", _plot_data);
}





void
Sweep::do_equilibrium(void)
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
Sweep::do_plot(void)
{
  int num_sim = _simulations.size();

  for (int i = 0; i < num_sim; i++)
    _simulations[i]->plot();
}



void
Sweep::do_solve(void)
{
  assert(_variable != "");

  parse_options();

  int num_sim = _simulations.size();

  assert(num_sim > 0);

  vector<map<double, vector<double> > > sweep_data(num_sim);

  // we make a plot file for each simulation
  vector<ofstream*> plotfiles(num_sim);

  bool do_plotting = prepare_plot_files(plotfiles);

  try
  {
    do_sweep(_values, plotfiles, sweep_data);
  }
  catch (SolveFailedException& e)
  {
    if (prepare_plot_files(plotfiles))
      plot_data(plotfiles, sweep_data);

    throw e;
  }

  if (prepare_plot_files(plotfiles))
    plot_data(plotfiles, sweep_data);
}




bool
Sweep::prepare_plot_files(vector<ofstream*>& plotfiles)
{
  int num_sim = _simulations.size();

  assert(plotfiles.size() == num_sim);

  bool do_plotting = false;

  // the filename suffix
  const string& suffix = get_control().get_filename_suffix();
  // the output directory
  string outdir = get_control().get_output_dir();

  for (int i = 0; i < num_sim; i++)
  {
    // some sanity check
    assert(_simulations[i] != NULL);
    assert(_simulations[i]->is_initialized());

    vector<string> legend;
    vector<string> description;
    _simulations[i]->get_integrated_quantities_description(
        legend, description);

    // should we plot something?
    if (legend.size() == 0)
      plotfiles[i] = NULL;
    else
    {
      if (plotfiles[i] != NULL)
      {
        plotfiles[i]->close();
        delete plotfiles[i];
      }

      ostringstream suff;
      //suff.precision(6);
      //suff << fixed << _variable;
      suff << _variable;
      string plotfilename(outdir + "/" + get_name() + "_" +
          _simulations[i]->get_name() + suffix + ".dat");

      plotfiles[i] = new ofstream(plotfilename.c_str(), ios_base::trunc);
      ofstream& file = *plotfiles[i];

      //file.open(plotfilename.c_str());

      if (!file.good())
        throw SolveFailedException("Sweep: Could not open plotfile " +
            plotfilename);


      //
      // print some header
      //
      ostringstream s;
      s << "# Parameter sweep " << "(" << get_name() << ")" << endl;
      s << "# Simulation: " << _simulations[i]->get_name() << endl;
      file << s.str();

      // print the data description
      file << "# Data:" << endl;
      for (unsigned int j = 0; j < description.size(); j++)
        file << "#    * " << description[j] << endl;


      ostringstream l;
      l << "#" << endl << "# " << _variable << "   ";
      unsigned int n = legend.size();
      for (unsigned int j = 0; j < n; j++)
        l << "  " << legend[j];
      l << endl;
      file << l.str();

      do_plotting = true;
    }
  }

  return do_plotting;
}




void
Sweep::plot_data(vector<ofstream*>& plotfiles,
    vector<map<double, vector<double> > >& sweep_data)
{
  int num_sim = _simulations.size();

  for (int j = 0; j < num_sim; j++)
  {
    if (plotfiles[j] != NULL)
    {
      // it means we have something to plot
      map<double, vector<double> >::iterator it(sweep_data[j].begin());
      map<double, vector<double> >::iterator end(sweep_data[j].end());

      for ( ; it != end; ++it)
      {
        vector<double>& plotvalues = it->second;

        ostringstream l;
        l << setprecision(12) << it->first;
        unsigned int n = plotvalues.size();

        for (unsigned int k = 0; k < n; k++)
          l << "   " << plotvalues[k];
        l << endl;

        ofstream& file = *plotfiles[j];
        file << l.str();
      }
    }
  }

  // clean up
  for (int j = 0; j < num_sim; j++)
  {
    if (plotfiles[j] != NULL)
    {
      plotfiles[j]->close();
      delete plotfiles[j];
    }
  }

}




void
Sweep::get_inner_simulations(
    std::vector<SimulationInterface*>& sims) const
{
  int num_sim = _simulations.size();
  sims.resize(0);

  for (unsigned int i = 0; i < num_sim; i++)
  {
    if (_simulations[i]->get_type() != "sweep")
      sims.push_back(_simulations[i]);
    else
    {
      vector<SimulationInterface*> inner;
      Sweep* sw = static_cast<Sweep*>(_simulations[i]);
      sw->get_inner_simulations(inner);
      sims.insert(sims.end(), inner.begin(), inner.end());
      break;
    }
  }
}




void
Sweep::do_sweep(vector<double>& values, vector<ofstream*>& plotfiles,
    vector<map<double, vector<double> > >& sweep_data)
{
  unsigned int n = values.size();

  // if there are no values, we return immediately
  if (n == 0) return;

  int num_sim = _simulations.size();

  // our ramp will solve directly the innermost simulations,
  // including only the first nested sweep
  int first_nested = num_sim;
  for (int j = 0; j < num_sim; j++)
  {
    if (_simulations[j]->get_type() == "sweep")
    {
      first_nested = j;
      break;
    }
  }

  // for plotting
  vector<double> plotvalues;

  vector<double>::iterator values_begin(values.begin());
  vector<double>::iterator values_end(values.end());

  vector<SimulationInterface*> inner_sims;
  get_inner_simulations(inner_sims);

  // for backup of data in case of nested sweeps
  map<int, ID> old_sol;

  // for the variable ramp
  Ramp ramp(get_options(), inner_sims);

  // remember the current solution
  remember_solution();

  for (unsigned int i = 0; i < n; i++)
  {
    double goal = values[i];

    {
      ostringstream os;
      os << "Ramping to sweep value " << _variable << " = " << goal;
      Messages::info(os.str());
    }

    // reset nested sweep data
    // NOTE: Sweep automatically remembers the starting solution
    // when entering do_solve()
    if (i > 0)
    {
      for (int j = 0; j < num_sim; j++)
      {
        if (_simulations[j]->get_type() == "sweep")
          _simulations[j]->set_to_remembered_solution(old_sol[j]);
      }
    }

    // filename suffix
    ostringstream suffix;
    suffix << _variable << "_" << goal;
    get_control().prepend_to_filename_suffix(suffix.str());
    
    bool failed = false;
    try
    {
      ramp.ramp(goal);

      // the loop over the simulations
      for (int j = 0; j < num_sim; j++)
      {
        // now solve everything after and including the 
        // first nested sweep
        if (j >= first_nested)
          _simulations[j]->solve();

        // plot results if required
        if (_plot_data)
          _simulations[j]->plot();

        // update "something-vs.-sweepvariable" files
        if (plotfiles[j] != NULL)
        {
          // it means we have something to plot
          _simulations[j]->get_integrated_quantities(plotvalues);

          sweep_data[j][goal] = plotvalues;

          ostringstream l;
          l << setprecision(12) << goal;
          unsigned int n_data = plotvalues.size();
          for (unsigned int k = 0; k < n_data; k++)
            l << "   " << plotvalues[k];
          l << endl;

          ofstream& file = *plotfiles[j];
          file << l.str();
          file.flush();
        }

      }
    }
    catch (SolveFailedException& e)
    {
      failed = true;
      Messages::error(e.what());
      if (i == 0)
        throw SolveFailedException("Already first sweep step could not be solved.");
    }

    get_control().drop_first_filename_suffix();

    if (failed)
    {
      ostringstream os;
      os << "Sweep failed at " << _variable << " = " << goal;
      throw SolveFailedException(os.str());
    }
  }
}



ID
Sweep::do_remember_current_solution(ID id)
{
  // we only remember the first solution!
  //return remember_solution();
  ignore_unused_variable(id);

  return 1;
}


ID
Sweep::remember_solution(void)
{
  ID id = 1;
  int num_sim = _simulations.size();

  map<ID, vector<ID> >::iterator end(_remembered_sol_ids.end());
  map<ID, vector<ID> >::iterator it(_remembered_sol_ids.find(id));

  if (it != end)
  {
    assert((it->second).size() == num_sim);
    for (int i = 0; i < num_sim; i++)
      (it->second)[i] =
        _simulations[i]->remember_current_solution((it->second)[i]);
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

  // the current variable value
  _remembered_sweep_value =
    Variable::get_variable_value<double>(_variable);

  return id;
}


void
Sweep::do_set_to_remembered_solution(ID id)
{
  int num_sim = _simulations.size();

  map<ID, vector<ID> >::iterator end(_remembered_sol_ids.end());
  map<ID, vector<ID> >::iterator it(_remembered_sol_ids.begin());

  if (it != end)
    for (int i = 0; i < num_sim; i++)
      _simulations[i]->set_to_remembered_solution((it->second)[i]);

  Variable::set_variable_value(_variable, _remembered_sweep_value);
}



void
Sweep::do_delete_remembered_solution(ID id)
{
  int num_sim = _simulations.size();

  map<ID, vector<ID> >::iterator end(_remembered_sol_ids.end());
  map<ID, vector<ID> >::iterator it(_remembered_sol_ids.find(id));

  if (it != end)
    for (int i = 0; i < num_sim; i++)
      _simulations[i]->delete_remembered_solution((it->second)[i]);
}


void
Sweep::build_integrated_quantities(
    vector<double>& values)
{
  vector<double> vals;

  int num_sim = _simulations.size();
  for (int i = 0; i < num_sim; i++)
  {
    _simulations[i]->get_integrated_quantities(vals);
    values.insert(values.end(), vals.begin(), vals.end());
  }
}



void
Sweep::build_integrated_quantities_description(
    vector<string>& legend,
    vector<string>& description)
{
  vector<string> leg;
  vector<string> desc;

  int num_sim = _simulations.size();
  for (int i = 0; i < num_sim; i++)
  {
    _simulations[i]->get_integrated_quantities_description(leg, desc);
    cerr << _simulations[i]->get_name() << " " << leg.size() << endl;
    legend.insert(legend.end(), leg.begin(), leg.end());
    description.insert(description.end(), desc.begin(), desc.end());
  }
}
