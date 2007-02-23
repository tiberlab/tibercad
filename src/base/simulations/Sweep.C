// $Id$

#include "Sweep.h"
#include "Sweepable.h"
#include "Boundary.h"
#include "BoundaryProperties.h"
#include "SimulationEnvironment.h"
#include "Control.h"

#include "fstream"

#include <boost/tokenizer.hpp>


using namespace std;

void
Sweep::do_init(void)
{

  using namespace boost;

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
  
  // get the model which contains the sweep variable
  // NOTE: for now we can sweep only over boundary values
  //
  // syntax is: 'simulationname.boundaryname'
  //         or 'boundaryname'
  string bnd = opts.get_option("boundary", "");
  typedef tokenizer<char_separator<char> > tokenizer;
  char_separator<char> sep(".");
  tokenizer tok(bnd, sep);
  tokenizer::iterator tokit(tok.begin());
  if (tokit == tok.end())
  {
    string msg("Sweep: You have to provide the name of ");
    msg += "the sweep variable (currently we support only boundary values).";
    throw InitFailedException(msg);
  }

  SimulationInterface* simulation;
  if ((++tokit) == tok.end())
  {
    // first check the sweep simulation
    SimulationEnvironment& env = get_environment();
    Boundary* boundary = env.get_boundary(bnd);
    if (boundary != NULL)
      simulation = _simulations[0];
    else
    {
      // only boundary name provided, let's hope it's unique...
      // (we just pick the first that matches)
      Control::simulation_iterator simit(control.simulations_begin());
      const Control::simulation_iterator simend(control.simulations_end());
      for ( ; simit != simend; ++simit)
      {
        simulation = *simit;
        if (simulation == this) // would not make sense...
          continue;

        SimulationEnvironment& env = simulation->get_environment();
        Boundary* boundary = env.get_boundary(bnd);
        if (boundary != NULL)
          break;
      }
    }
  }
  else
  {
    tokit = tok.begin();
    simulation = control.find_simulation(*tokit);
    if (simulation == NULL)
      throw InitFailedException("Sweep: Simulation " + *tokit + " not found.");

    bnd = *(++tokit);
    
  }


  // simulation now should have the sweepable boundary value
  SimulationEnvironment& env = simulation->get_environment();


  Boundary* boundary = env.get_boundary(bnd);
  if (boundary == NULL)
    throw InitFailedException("Sweep: Boundary " + bnd + " not found.");
  
  _variable = dynamic_cast<Sweepable*>(
      boundary->get_boundary_properties(simulation->get_id()));
  
  if (_variable == NULL)
    throw InitFailedException("Sweep: No sweepable entity on boundary " +
        bnd + " found.");
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

  _values.resize(steps + 1);
  double step = (stop - start) / steps;
  for (unsigned int i = 0; i <= steps; i++)
    _values[i] = start + i * step;

  // we can also specify a vector with the values
  opts.get_option("values", _values);


  //! read the variables we want to plot (type IV characteristic)
  vector<string> vars;
  opts.get_option("plotvariable", vars);
  for (unsigned int i = 0; i < vars.size(); i++)
    _plotvariables.insert(vars[i]);



  // unused
  //_do_output = opts.get_option("output", _do_output);
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
Sweep::do_solve(void)
{
  assert(_variable != NULL);

  parse_options();

  // the current filename suffix
  string suffix = get_control().get_filename_suffix();
  string outdir = get_control().get_output_dir();
  
  int num_sim = _simulations.size();

  assert(num_sim > 0);

  // we make a plot file for each simulation
  vector<double> values;
  bool do_plot = false;
  vector<ofstream*> plotfiles(num_sim);
  for (int i = 0; i < num_sim; i++)
  {
    // some sanity check
    assert(_simulations[i] != NULL);
    assert(_simulations[i]->is_initialized());
    
    vector<string> legend;
    vector<string> description;
    _simulations[i]->get_integrated_quantities_description(_plotvariables,
        legend, description);

    // should we plot something?
    if (legend.size() == 0)
      plotfiles[i] = NULL;
    else
    {

      string plotfilename(outdir + "/" + get_name() + "_" +
          _simulations[i]->get_name() + suffix + ".dat");

      plotfiles[i] = new ofstream;
      ofstream& file = *plotfiles[i];

      file.open(plotfilename.c_str());

      if (!file.good())
        throw SolveFailedException("Sweep: Could not open plotfile " +
            plotfilename);


      //
      // print some header
      // 
      ostringstream s;
      s << "# Parameter sweep " << "(" << ")" << endl;
      s << "# Simulation: " << _simulations[i]->get_name() << endl;
      file << s.str();

      // print the data description
      file << "# Data:" << endl;
      for (unsigned int j = 0; j < description.size(); j++)
        file << "#    * " << description[j] << endl;

      
      ostringstream l;
      l << "#" << endl << "# x   "; // x is out sweep variable
      unsigned int n = legend.size();
      for (unsigned int j = 0; j < n; j++)
        l << "  " << legend[j];
      l << endl;
      file << l.str();

      do_plot = true;
    }
  }


  // we make a copy of the current solutions
  // we need this in the case of a solver failure to go back
  // to an old successful solution
  vector<NumericVector<Real>* > old_sol(num_sim); 
  for (int i = 0; i < num_sim; i++)
    old_sol[i] = ((_simulations[i]->get_solution_vector()).clone()).release();



  //
  // now do the loop
  //

  vector<double>::iterator values_begin(_values.begin());
  vector<double>::iterator values_end(_values.end());

  unsigned int n = _values.size();
  // we look at voltages +/- 1e-9 V
  double eps = 1e-9;

  for (unsigned int i = 0; i < n; i++)
  {
    double goal = _values[i];
    double goal_sign = (goal < 0.0) ? -1 : 1;
    _last = _variable->get_current_value();
    double step = goal - _last;
    double old_step = 0.0;

    // we iterate until we arrive at the goal
    double value;
    do
    {
      // check for step > max_step
      double absstep = abs(step);
      double sign = step < 0.0 ? -1 : 1;
      step = (absstep > _max_step) ? _max_step * sign : step;

      value = _last + step;
      double diff = goal_sign * (goal - value);
      if (diff < 0.0)
        value = goal;
      
      _variable->set_new_value(value);
      cout  << "Sweep value = " << value << endl;
      
      try
      {
        for (int j = 0; j < num_sim; j++)
          _simulations[j]->solve();


        _last = _variable->get_current_value();


        // remember the current solution
        for (int j = 0; j < num_sim; j++)
          *(old_sol[j]) = _simulations[j]->get_solution_vector();


        // write results, but only at desired sweep steps
        if (find_if(values_begin, values_end,
              bind2nd(Utils::almost_equal(), _last)) != values_end)
        {
          ostringstream s;
          s << suffix << _last;
          get_control().set_filename_suffix(s.str());

          for (int j = 0; j < num_sim; j++)
            _simulations[j]->plot();
        }


        // update "something-vs.-sweepvariable" files
        // Here we do it also for intermediate steps
        if (do_plot)
        {
          for (int j = 0; j < num_sim; j++)
          {
            if (plotfiles[j] != NULL)
              // it means we have something to plot
            {
              ofstream& file = *plotfiles[j];
              _simulations[j]->get_integrated_quantities(_plotvariables, values);
              ostringstream l;
              l << setprecision(12) << _last;
              unsigned int n = values.size();
              for (unsigned int k = 0; k < n; k++)
                l << "   " << values[k];
              l << endl;
              
              file << l.str() << flush;
            }
          }
        }


        // try to increase step, but only if it worked twice
        // with the old one
        if (step == old_step)
          step *= 2.0;
        old_step = step;
      }
      catch (...)
      {
        step = (value - _last) / 2.0;
        if (abs(step) < _min_step)
          throw SolveFailedException("Sweep: step size small.");

        // set to the remembered solution
        for (int j = 0; j < num_sim; j++)
          _simulations[j]->get_solution_vector() = *(old_sol[j]);

        // it failed, so we go back to the old value
        // (In principle this is unnecessary, but we need it to not
        // get out of the do...while loop
        value = _last;
      }
    }
    while (goal_sign * (value - goal) < -eps);

  }


  // clean up
  for (int j = 0; j < num_sim; j++)
  {
    delete old_sol[j];

    if (plotfiles[j] != NULL)
    {
      plotfiles[j]->close();
      delete plotfiles[j];
    }
  }
}

