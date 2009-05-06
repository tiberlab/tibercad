// $Id$

#include "Sweep.h"
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




const set<string>&
Sweep::get_plotvariables(void) const
{
  if (_plotvariables.size() != 0)
    return _plotvariables;

  return get_control().get_plotvariables();
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


  // read the variables we want to plot (type IV characteristic)
  vector<string> vars;
  opts.get_option("plotvariable", vars);
  for (unsigned int i = 0; i < vars.size(); i++)
    _plotvariables.insert(vars[i]);



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

  // the current filename suffix
  string suffix = get_control().get_filename_suffix();
  string outdir = get_control().get_output_dir();

  int num_sim = _simulations.size();

  assert(num_sim > 0);

  vector<map<double, vector<double> > > sweep_data(num_sim);

  // we make a plot file for each simulation
  vector<ofstream*> plotfiles(num_sim);

  bool do_plotting = prepare_plot_files(plotfiles);

  // if there are negative and positive values, we split them apart
  vector<double> pos_values;
  vector<double> neg_values;

  unsigned int n = _values.size();

  try
  {
    for (unsigned int i = 0; i < n; i++)
    {
      if (_values[i] >= 0.0)
        pos_values.push_back(_values[i]);
      else
        neg_values.push_back(_values[i]);
    }

    if (neg_values.size() == 0)
      do_sweep(pos_values, plotfiles, sweep_data);
    else
    {
      if (pos_values.size() > 1)
      {
        // in this case we will do: 0 -> max, 0 -> min
        if (neg_values[0] < neg_values[neg_values.size() - 1])
          reverse(neg_values.begin(), neg_values.end());

        vector<double> first_val(1, pos_values[0]);
        pos_values.erase(pos_values.begin());

        do_sweep(first_val, plotfiles, sweep_data);
        // to remember the first solution
        vector<ID> old_sol(num_sim);
        for (int i = 0; i < num_sim; i++)
          old_sol[i] = _simulations[i]->remember_current_solution();

        do_sweep(pos_values, plotfiles, sweep_data);

        // reset to the first solution
        for (int i = 0; i < num_sim; i++)
        {
          _simulations[i]->set_to_remembered_solution(old_sol[i]);
          _simulations[i]->delete_remembered_solution(old_sol[i]);
        }

        Variable::set_variable_value(_variable, first_val[0]);
      }
      else
        do_sweep(pos_values, plotfiles, sweep_data);


      do_sweep(neg_values, plotfiles, sweep_data);
    }
  }
  catch (SolveFailedException& e)
  {
    get_control().set_filename_suffix(suffix);
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

  // the current filename suffix
  string suffix = get_control().get_filename_suffix();
  // the output directory
  string outdir = get_control().get_output_dir();

  for (int i = 0; i < num_sim; i++)
  {
    // some sanity check
    assert(_simulations[i] != NULL);
    assert(_simulations[i]->is_initialized());

    vector<string> legend;
    vector<string> description;
    _simulations[i]->get_integrated_quantities_description(get_plotvariables(),
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
          _simulations[i]->get_name() + suffix + "_" + suff.str() + ".dat");

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
Sweep::do_sweep(vector<double>& values, vector<ofstream*>& plotfiles,
    vector<map<double, vector<double> > >& sweep_data)
{
  unsigned int n = values.size();

  // if there are no values, we return immediately
  if (n == 0) return;

  int num_sim = _simulations.size();

  // the current filename suffix
  string suffix = get_control().get_filename_suffix();


  // we make a copy of the current solutions
  // we need this in the case of a solver failure to go back
  // to an old successful solution
  vector<ID> old_sol(num_sim);
  for (int i = 0; i < num_sim; i++)
    old_sol[i] = _simulations[i]->remember_current_solution();

  // for plotting
  vector<double> plotvalues;

  vector<double>::iterator values_begin(values.begin());
  vector<double>::iterator values_end(values.end());

  // we look at voltages +/- 1e-9 V
  double eps = 1e-9;

  for (unsigned int i = 0; i < n; i++)
  {
    double goal = values[i];

    _last = Variable::get_variable_value(_variable);
    double step = goal - _last;
    double step_sign = (step < 0.0) ? -1 : 1;
    double old_step = 0.0;


    for (int j = 0; j < num_sim; j++)
      _simulations[j]->set_to_remembered_solution(old_sol[j]);


    // we iterate until we arrive at the goal
    double value;
    do
    {
      // check for step > max_step
      double absstep = abs(step);
      step = (absstep > _max_step) ? _max_step * step_sign : step;

      value = _last + step;
      double diff = step_sign * (goal - value);
      if (diff < 0.0)
        value = goal;

      Variable::set_variable_value(_variable, value);
      {
        ostringstream os;
        os << "Sweep value " << _variable << " = " << value;
        Messages::info(os.str());
      }

      try
      {

        // prepare filename suffix
        {
          ostringstream s;
          //s.precision(3);
          //s << suffix << "_" << _variable << "_" << fixed << value;
          s << suffix << "_" << _variable << "_" << value;
          get_control().set_filename_suffix(s.str());
        }

        bool plot_data = false;
        // write results, but only at desired sweep steps
        if (_plot_data)
          if (find_if(values_begin, values_end,
                bind2nd(Utils::almost_equal(), value)) != values_end)
          {
            plot_data = true;
          }


        // the loop over the simulations
        for (int j = 0; j < num_sim; j++)
        {
          _simulations[j]->solve();

          // update "something-vs.-sweepvariable" files
          // Here we do it also for intermediate steps
          if (plotfiles[j] != NULL)
          {
            // it means we have something to plot
            _simulations[j]->get_integrated_quantities(get_plotvariables(),
                plotvalues);

            sweep_data[j][value] = plotvalues;

            ostringstream l;
            l << setprecision(12) << value;
            unsigned int n_data = plotvalues.size();
            for (unsigned int k = 0; k < n_data; k++)
              l << "   " << plotvalues[k];
            l << endl;

            ofstream& file = *plotfiles[j];
            file << l.str();
            file.flush();
          }

          if (plot_data)
            _simulations[j]->plot();


          // remember the current solution
          // for sweeps we return to the former solution!
          if (_simulations[j]->get_type() != "sweep")
            _simulations[j]->remember_current_solution(old_sol[j]);
          else
            _simulations[j]->set_to_remembered_solution(old_sol[j]);
        }



        _last = value;



        // try to increase step, but only if it worked twice
        // with the old one
        if (step == old_step)
          step *= 2.0;
        old_step = step;
      }
      catch (SolveFailedException& e)
      {
        {
          ostringstream os;
          os << "Solve failed due to: " << Messages::endl
            << "   " << e.what();
          Messages::warning(os.str());
        }
        step = (value - _last) / 2.0;
        if (abs(step) < _min_step)
          throw SolveFailedException("Step size too small in sweep.");

        if (i == 0)
          throw SolveFailedException("Failure in first sweep step.");

        Messages::warning("Trying intermediate step");

        // set to the remembered solution
        for (int j = 0; j < num_sim; j++)
          _simulations[j]->set_to_remembered_solution(old_sol[j]);

        // it failed, so we go back to the old value
        // (In principle this is unnecessary, but we need it to not
        // get out of the do...while loop
        value = _last;
      }
    }
    while (step_sign * (value - goal) < -eps);


    if (i == 0)
      remember_solution();

  }

  get_control().set_filename_suffix(suffix);

  // clean up
  for (int j = 0; j < num_sim; j++)
    _simulations[j]->delete_remembered_solution(old_sol[j]);

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
Sweep::do_set_to_remembered_solution(ID id)
{
  int num_sim = _simulations.size();

  map<ID, vector<ID> >::iterator end(_remembered_sol_ids.end());
  map<ID, vector<ID> >::iterator it(_remembered_sol_ids.begin());

  if (it != end)
    for (int i = 0; i < num_sim; i++)
      _simulations[i]->set_to_remembered_solution((it->second)[i]);
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
    const set<string>& variables,
    vector<double>& values)
{
  vector<double> vals;

  int num_sim = _simulations.size();
  for (int i = 0; i < num_sim; i++)
  {
    _simulations[i]->get_integrated_quantities(variables, vals);
    values.insert(values.end(), vals.begin(), vals.end());
  }
}



void
Sweep::build_integrated_quantities_description(
    const set<string>& variables,
    vector<string>& legend,
    vector<string>& description)
{
  vector<string> leg;
  vector<string> desc;

  int num_sim = _simulations.size();
  for (int i = 0; i < num_sim; i++)
  {
    _simulations[i]->get_integrated_quantities_description(variables,
        leg, desc);
    legend.insert(legend.end(), leg.begin(), leg.end());
    description.insert(description.end(), desc.begin(), desc.end());
  }
}
