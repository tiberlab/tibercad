// $Id$

#include "Ramp.h"
#include "Control.h"
#include "SimulationInterface.h"
#include "Variable.h"
#include "Messages.h"


using namespace std;



Ramp::Ramp(const ModelOptions& options,
    const vector<SimulationInterface*>& simulations)
  : _simulations(simulations),
    _variable(""),
    _goal(0.0),
    _last(0.0),
    _initial_step(1.0),
    _min_step(1e-3),
    _max_step(1.0),
    _plot_data(false)
{

  if (simulations.size() == 0)
  {
    // get the names of the simulations to be solved
    vector<string> sims;
    options.get_option("simulation", sims);
    int num_of_sims = sims.size();

    Control& control = TiberCad::get_control();

    _simulations.resize(num_of_sims);
    for (int i = 0; i < num_of_sims; i++)
    {
      _simulations[i] = control.find_simulation(sims[i]);
      if (_simulations[i] == NULL)
        throw InitFailedException("Sweep: Simulation " + sims[i]
            + " not found.");

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
  }

  //
  // at this point we have for sure one simulation
  //


  // Now we have to find the model to the variable
  _variable = options.get_option("variable", "");
  if (_variable == "")
  {
    string msg("Sweep: You have to provide the name of ");
    msg += "the sweep variable.";
    throw InitFailedException(msg);
  }


  _goal = options.get_option("goal", _goal);
  _min_step = options.get_option("min_step", _min_step);
  _min_step = min(max(_min_step, 0.0), 1.0);
  _max_step = options.get_option("max_step", _max_step);
  _max_step = min(max(_max_step, 0.0), 1.0);
  _initial_step = options.get_option("initial_step", _initial_step);
  _initial_step = min(max(_initial_step, 0.0), 1.0);
  _plot_data = options.get_option("plot_data", _plot_data);

}


void
Ramp::ramp(void)
{
  int num_sim = _simulations.size();

  if (_old_sol_ids.size() == 0)
  {
    // we make a copy of the current solutions
    // we need this in the case of a solver failure to go back
    // to an old successful solution
    _old_sol_ids.resize(num_sim);
    for (int i = 0; i < num_sim; i++)
    {
      if (!_simulations[i]->is_solved())
      {
         ostringstream os;
         os << "Will presolve a currently unsolved system ("
            << _simulations[i]->get_name() << ")";
         Messages::warning(os.str());
      }
      _old_sol_ids[i] = _simulations[i]->remember_current_solution();
    }
  }
  else
    for (int i = 0; i < num_sim; i++)
      _simulations[i]->remember_current_solution(_old_sol_ids[i]);


  _last = Variable::get_variable_value<double>(_variable);

  // NOTE: we will treat the absolute value of the steps and their
  // sign independently
  //
  // the total step
  double step = _goal - _last;
  double sign = (step < 0) ? -1.0 : 1.0;
  step *= sign;

  double min_step = _min_step * step;
  double max_step = _max_step * step;
  double initial_step = _initial_step * step;

  double oldvalue = _last;
  double value = oldvalue;
  double oldstep = 0;
  double currstep = min(max_step, initial_step);

  do
  {
    if (currstep < min_step)
    {
      ostringstream os;
      os << "Ramp of " << _variable << " failed: step size too small.";
      throw SolveFailedException(os.str());
    }

    // update the value
    value = sign * min(sign * _goal, sign * value + currstep);

    Variable::set_variable_value(_variable, value);

    // we define the simulation index here so we can know in the
    // catch clause which simulation failed.
    int j;
    try
    {
      ostringstream os;
      os << "Trying " << _variable << " = " << value;
      Messages::info(os.str());

      for (j = 0; j < num_sim; j++)
        _simulations[j]->solve();

      // we have to remember the current solution
      for (int i = 0; i < num_sim; i++)
        _simulations[i]->remember_current_solution(_old_sol_ids[i]);

      oldvalue = value;
      double factor = (currstep == oldstep) ? 2.0 : 1.0;
      oldstep = currstep;
      currstep = min(factor * currstep, max_step);
    }
    catch (SolveFailedException& e)
    {
      for (int i = 0; i < num_sim; i++)
        _simulations[i]->set_to_remembered_solution(_old_sol_ids[i]);

      currstep /= 2.0;

      ostringstream os;
      os << "Ramping to " << _variable << " = " << value
        << " failed while solving \'" << _simulations[j]->get_name() << "\'.";
      Messages::warning(os.str());

      // this is a special case
      if (step == 0) throw SolveFailedException("Cannot decrease step size, is already 0");

      value = oldvalue;
    }
  }
  while (sign * (_goal - value) > 0);


  // plot results
  if (_plot_data)
    for (int i = 0; i < num_sim; i++)
      _simulations[i]->plot();

}


Ramp::~Ramp(void)
{
  int num_sim = _old_sol_ids.size();
  for (int j = 0; j < num_sim; j++)
    _simulations[j]->delete_remembered_solution(_old_sol_ids[j]);
}


