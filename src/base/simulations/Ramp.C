// $Id$

#include "Ramp.h"
#include "SimulationInterface.h"
#include "Variable.h"
#include "Messages.h"

#include <limits>

using namespace std;



Ramp::Ramp(const ModelOptions& options,
    const vector<SimulationInterface*>& simulations)
  : _simulations(simulations)
{

  if (simulations.size() == 0)
  {
    // get the names of the simulations to be solved
    vector<string> sims;
    options.get_option("simulation", sims);
    int num_of_sims = sims.size();

    _simulations.resize(num_of_sims);
    for (int i = 0; i < num_of_sims; i++)
    {
      _simulations[i] = SimulationInterface::find_simulation(sims[i]);
      if (_simulations[i] == NULL)
        throw InitFailedException("Sweep: Simulation " + sims[i]
            + " not found.");

      // If it is not already initialized, we initialize now
      if (!_simulations[i]->is_initialized())
        _simulations[i]->init();
    }

    /*
    // if user didn't provide a simulation name, we take the first available
    if (num_of_sims == 0)
    {
      _simulations.resize(1);
      _simulations[0] = SimulationInterface::find_simulation("");
      if (_simulations[0] == NULL)
        throw InitFailedException("Sweep: No simulation found.");
    }
    */
  }

  //
  // at this point we have for sure one simulation
  //


  // Now we have to find the model to the variable
  string variable = options.get_option("variable", "");
  if (variable == "")
  {
    string msg("Sweep: You have to provide the name of ");
    msg += "the sweep variable.";
    throw InitFailedException(msg);
  }

  Utils::extract_vector(variable, _variable);


  _goal = options.get_option("goal", _goal);

  _min_step = options.get_option("min_relative_step", _min_step);
  _min_step = min(max(_min_step, 0.0), 1.0);
  _min_abs_step = options.get_option("min_step", _min_abs_step);

  _max_step = options.get_option("max_relative_step", _max_step);
  _max_step = min(max(_max_step, 0.0), 1.0);
  _max_abs_step = options.get_option("max_step", _max_abs_step);

  _initial_step = options.get_option("initial_relative_step", _initial_step);
  _initial_step = min(max(_initial_step, 0.0), 1.0);
  _initial_abs_step = options.get_option("initial_step", _initial_abs_step);

  //_plot_data = options.get_option("plot_data", _plot_data);

}



void
Ramp::ramp(double goal)
{
  _goal = goal;
  ramp();
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
         //os << "INFO: I will first solve a currently unsolved system ("
         //   << _simulations[i]->get_name() << "):";
         //Messages::info(os.str());
         //_simulations[i]->solve();
         os << _simulations[i]->get_name() << " has never been solved before ramping.";
         Messages::warning(os.str());
      }
      _old_sol_ids[i] = _simulations[i]->remember_current_solution();
    }
  }
  else
    for (int i = 0; i < num_sim; i++)
      _simulations[i]->remember_current_solution(_old_sol_ids[i]);


  _last = VariableValue::get_variable_value<double>(_variable[0]);

  // NOTE: we will treat the absolute value of the steps and their
  // sign independently
  //
  // the total step
  double step = _goal - _last;
  double sign = (step < 0) ? -1.0 : 1.0;
  step *= sign;

  double min_step = max(_min_abs_step, _min_step * step);
  double max_step = min(_max_abs_step, _max_step * step);
  double initial_step = min(_initial_step * step, _initial_abs_step);
  //double min_step = _min_step * step;
  //double max_step = _max_step * step;
  //double initial_step = _initial_step * step;

  double oldvalue = _last;
  double value = oldvalue;
  double oldstep = 0;
  double currstep = min(max_step, initial_step);

  do
  {
    if ((currstep < min_step) && (step > 0))
    {
      ostringstream os;
      os << "Ramp of " << _variable[0];
      if (_variable.size() > 1)
        os << "...";
      os << " failed: step size too small.";
      throw SolveFailedException(os.str());
    }

    // update the value
    value = sign * min(sign * _goal, sign * value + currstep);

    for (auto&& var : _variable)
      VariableValue::set_variable_value(var, value);

    // we define the simulation index here so we can know in the
    // catch clause which simulation failed.
    int j;
    try
    {
      ostringstream os;
      os << "Trying " << _variable[0];
      for (unsigned int i = 1; i < _variable.size(); ++i)
        os << "," << _variable[i];
      os << " = " << value;
      Messages::info(os.str());

      for (j = 0; j < num_sim; j++)
        _simulations[j]->solve();

      // we have to remember the current solution
      for (int i = 0; i < num_sim; i++)
        _simulations[i]->remember_current_solution(_old_sol_ids[i]);

      oldvalue = value;
      double factor = 2;
      //double factor = (currstep == oldstep) ? 2.0 : 1.0;
      oldstep = currstep;
      currstep = min(factor * currstep, max_step);
    }
    catch (SolveFailedException& e)
    {
      for (int i = 0; i < num_sim; i++)
      {
        _simulations[i]->set_to_remembered_solution(_old_sol_ids[i]);
        //TiberCad::prepend_to_filename_suffix("_reloaded");
        //_simulations[i]->plot();
        //TiberCad::drop_first_filename_suffix();
      }

      currstep /= 2.0;

      ostringstream os;
      os << "Ramping to " << _variable[0];
      if (_variable.size() > 1)
        os << "...";
      os << " = " << value
         << " failed while solving \'" << _simulations[j]->get_name() << "\'.";

      if (_simulations[j]->verbose() > 1)
        os << "\n" << e.what();

      Messages::warning(os.str());

      // this is a special case
      if (step == 0) throw SolveFailedException("Cannot decrease step size, is already 0");

      value = oldvalue;
      for (auto&& var : _variable)
        VariableValue::set_variable_value(var, value);
      //for (int i = 0; i < num_sim; i++)
      //  _simulations[i]->solve();
    }
  }
  while (abs(_goal - value) > 1e9 * (abs(_goal) * numeric_limits<double>::epsilon()
      + numeric_limits<double>::min()));
  // the epsilon() here prevents from resolving two times the same sweep value
  // due to fixed point numerics

  // At the next call, we begin with the last successful step size,
  // if this is smaller than the original initial step
  //if (currstep >= min_step)
  //  _initial_abs_step = oldstep;
  //  _initial_abs_step = min(_initial_abs_step, oldstep);

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


