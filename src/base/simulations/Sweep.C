// $Id$

#include "Sweep.h"
#include "Sweepable.h"
#include "Boundary.h"
#include "BoundaryProperties.h"
#include "SimulationEnvironment.h"
#include "Control.h"

#include "fstream"


using namespace std;

void
Sweep::do_init(void)
{
  ModelOptions& opts = get_options();
  
  // get the simulation
  const string& sim = opts.get_option("simulation", "");
  _simulation = SimulationInterface::find_simulation(sim);
  
  // we don't tolerate NULL pointers...
  if (_simulation == NULL)
    throw InitFailedException("Sweep: Simulation " + sim + " not found.");

  if (!_simulation->is_initialized())
    _simulation->init();
  
  // get the model which contains the sweep variable
  // NOTE: for now we can sweep over boundary values
  const string& bnd = opts.get_option("boundary", "");
  SimulationEnvironment& env = _simulation->get_environment();

  // we have to set our environments to enable nested sweeps!
  set_environment(&env);

  if (!env.is_initialized())
    env.init();

  Boundary* boundary = env.get_boundary(bnd);
  if (boundary == NULL)
    throw InitFailedException("Sweep: Boundary " + bnd + " not found.");
  
  _variable = dynamic_cast<Sweepable*>(
      boundary->get_boundary_properties(_simulation->get_id()));
  
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
Sweep::do_solve(void)
{
  assert(_simulation != NULL);
  assert(_simulation->is_initialized());
  assert(_variable != NULL);

  parse_options();

  // the current filename suffix
  string suffix = get_control().get_filename_suffix();
  string outdir = get_control().get_output_dir();
  
  vector<double> values;
  bool do_plot = false;
  ofstream plotfile;
  {
    vector<string> legend;
    vector<string> description;
    _simulation->get_integrated_quantities_description(_plotvariables,
        legend, description);
    // should we plot something?
    if (legend.size() != 0)
    {

      string plotfilename(outdir + "/" + get_name() + suffix + ".dat");
      plotfile.open(plotfilename.c_str());

      if (!plotfile.good())
        throw SolveFailedException("Sweep: Could not open plotfile " +
            plotfilename);


      //
      // print some header
      // 
      ostringstream s;
      s << "# Parameter sweep " << "(" << ")" << endl;
      s << "# Simulation: " << _simulation->get_name() << endl;
      plotfile << s.str();

      // print the data description
      plotfile << "# Data:" << endl;
      for (unsigned int i = 0; i < description.size(); i++)
        plotfile << "#    " << description[i] << endl;

      
      ostringstream l;
      l << "#" << endl << "# x   ";
      unsigned int n = legend.size();
      for (unsigned int i = 0; i < n; i++)
        l << "  " << legend[i];
      l << endl;
      plotfile << l.str();

      do_plot = true;
    }
  }

  // we make a copy of the current solution
  // we need this in the case of a solver failure to go back
  // to an old successful solution
  AutoPtr<NumericVector<Real> > old_sol = 
    (_simulation->get_solution_vector()).clone();

  unsigned int n = _values.size();

  for (unsigned int i = 0; i < n; i++)
  {
    double goal = _values[i];
    double goal_sign = (goal < 0.0) ? -1 : 1;
    _last = _variable->get_current_value();
    double step = goal - _last;
    double old_step = 0.0;

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
      cerr  << "Sweep value = " << value << endl;
      
      try
      {
        _simulation->solve();

        _last = _variable->get_current_value();

        // write results
        ostringstream s;
        s << suffix << _last;
        get_control().set_filename_suffix(s.str());
        _simulation->plot();


        // update "something-vs.-sweepvariable" file
        if (do_plot)
        {
          _simulation->get_integrated_quantities(_plotvariables, values);
          ostringstream l;
          l << setprecision(12) << _last;
          unsigned int n = values.size();
          for (unsigned int i = 0; i < n; i++)
            l << "   " << values[i];
          l << endl;
          plotfile << l.str();
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

        _simulation->get_solution_vector() = *old_sol;

        // it failed, so we go back to the old value
        // (In principle this is unnecessary, but we need it to not
        // get out of the do...while loop
        value = _last;
      }
    }
    while (goal_sign * (value - goal) < 0.0);

    *old_sol = _simulation->get_solution_vector();
  }

  plotfile.close();
  plotfile.close();
}

