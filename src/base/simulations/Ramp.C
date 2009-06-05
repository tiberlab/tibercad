// $Id$

#include "Ramp.h"
#include "Control.h"
#include "SimulationInterface.h"
#include "Variable.h"
#include "Messages.h"


using namespace std;


Ramp::Ramp(const ModelOptions& options)
  : _variable(""),
    _goal(0.0),
    _last(0.0),
    _initial_step(1.0),
    _min_step(1e-3),
    _max_step(1.0)
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
}


void
Ramp::ramp(void)
{
  int num_sim = _simulations.size();

  // we make a copy of the current solutions
  // we need this in the case of a solver failure to go back
  // to an old successful solution
  vector<ID> old_sol(num_sim);

  _last = Variable::get_variable_value<double>(_variable);
  double step = _goal - _last;
  double min_step = _min_step * step;
  double max_step = _max_step * step;
  double initial_step = _initial_step * step;

  double oldvalue = _last;
  double oldstep = 0;
  double currstep = min(max_step, initial_step);
  double value = oldvalue + currstep;

  do
  {
    if (currstep < min_step)
    {
      ostringstream os;
      os << "Ramp of " << _variable << " failed: step size too small.";
      throw SolveFailedException(os.str());
    }

    for (int i = 0; i < num_sim; i++)
      old_sol[i] = _simulations[i]->remember_current_solution();

    Variable::set_variable_value(_variable, value);

    try
    {
      for (int j = 0; j < num_sim; j++)
        _simulations[j]->solve();


      oldvalue = value;
      oldstep = currstep;
      currstep = max(2 * currstep, max_step);
      value += currstep;
    }
    catch (SolveFailedException& e)
    {
      {
        ostringstream os;
        os << "Solve of failed due to: " << Messages::endl
        << "   " << e.what();
        Messages::warning(os.str());
      }

      currstep /= 2;
      value = oldvalue + currstep;
    }
  }
  while (value < _goal);

}

