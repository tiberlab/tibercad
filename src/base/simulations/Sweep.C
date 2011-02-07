// $Id$

#include "Sweep.h"
#include "Ramp.h"
#include "Utils.h"
#include "Variable.h"
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

  // get the names of the simulations to be solved
  vector<string> sims;
  opts.get_option("solve", sims);
  int num_of_sims = sims.size();


  _simulations.resize(num_of_sims);
  for (int i = 0; i < num_of_sims; i++)
  {
    _simulations[i] = find_simulation(sims[i]);
    if (_simulations[i] == NULL)
      throw InitFailedException("Sweep: Simulation " + sims[i] + " not found.");

    // If it is not already initialized, we initialize now
    if (!_simulations[i]->is_initialized())
      _simulations[i]->init();


  }

  // if user didn't provide a simulation name, we take the first available
  if (num_of_sims == 0)
  {
    SimulationIterator it = simulations_begin();
    if (it == simulations_end())
      throw InitFailedException("Sweep: No simulation found.");

    _simulations.resize(1);
    _simulations[0] = *it;
  }

  //
  // at this point we have for sure one simulation
  //

  // we set our environment to that of the first simulation
  // 02-09-2010 Sweep should not have an environment
  //set_environment(&_simulations[0]->get_environment());


  // Now we have to find the model to the variable
  _variable = opts.get_option("variable", "");
  if (_variable == "")
  {
    string msg("Sweep: You have to provide the name of ");
    msg += "the sweep variable.";
    throw InitFailedException(msg);
  }

  parse_options();
}




void
Sweep::do_print_info(void)
{
  // this is a dirty trick to not have the variables in the command
  // line output

  for (int i = 0; i < _simulations.size(); i++)
  {
    // inherit the solution IDs
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



void
Sweep::parse_options(void)
{
  ModelOptions& opts = get_options();

  //_min_step = opts.get_option("min_step", _min_step);
  //_max_step = opts.get_option("max_step", _max_step);

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

  do_sweep(_values, plotfiles, sweep_data);

  // clean up
  for (size_t i = 0; i < num_sim; ++i)
    delete plotfiles[i];
}





void
Sweep::write_global_data(SimulationInterface& simulation, ofstream*& plotfile)
{
  map<ID, vector<double> > data;
  simulation.get_solution(data);

  if (data.size() > 0)
  {
    if (plotfile == NULL)
    {
      //
      // prepare the file
      //

      // the output directory
      string outdir = get_output_directory();

      string plotfilename(outdir + "/" + get_name() + "_" +
          simulation.get_name() + _suffix + ".dat");

      plotfile = new ofstream(plotfilename.c_str(), ios_base::trunc);
      ofstream& file = *plotfile;


      if (!file.good())
        throw SolveFailedException("Sweep: Could not open plotfile " +
            plotfilename);


      //
      // print some header
      //
      ostringstream s;
      s << "# Parameter sweep " << "(" << get_name() << ")" << endl;
      s << "# Simulation    : " << simulation.get_name() << endl;
      //s << "# Device        : " << simulation.get_device()->get_name() << endl;
      s << "# Sweep variable: " << _variable << endl;
      s << "# Data:" << endl;
      file << s.str();

      int width[3] = {25, 15, 10};
      int tot_width = width[0] + width[1] + width[2];

      ostringstream line;
      line.width(tot_width);
      line.fill('-');
      line << "";
      file << "# " << line.str() << endl;

      {
        ostringstream os;
        os << "# Name";
        int w = width[0] + 2;
        os.width(w - os.tellp());
        os << "" << "Units";
        w += width[1];
        os.width(w - os.tellp());
        os << "" << "Type";
        w += width[2];
        file << os.str() << "\n# ";
        file << line.str() << "\n# ";
      }

      ostringstream l;
      l << "#" << endl << "# " << _variable << " ";
      map<ID, vector<double> >::iterator it = data.begin();
      const map<ID, vector<double> >::iterator end = data.end();
      for ( ; it != end; ++it)
      {
        const SolutionDescriptor& descr = simulation.get_solution_descriptor(it->first);
        const string& name = descr.name();

        int w = width[0];
        ostringstream os;
        os << name;
        os.width(w - os.tellp());
        os << "" << descr.units();
        w += width[1];
        os.width(w - os.tellp());
        os << "" << descr.type();
        w += width[2];
        file << os.str() << "\n# ";


        switch (descr.type())
        {
          case SolutionDescriptor::VECTOR:
            l << name << "_x "
            << name << "_y "
            << name << "_z ";
            break;

          case SolutionDescriptor::TENSOR:
            l << name << "_11 "
            << name << "_22 "
            << name << "_33 "
            << name << "_12 "
            << name << "_23 "
            << name << "_13 ";
            break;

          case SolutionDescriptor::NTUPLE:
            for (unsigned int i = 0; i < descr.n_components(); i++)
              l << name << "_" << i << " ";
            break;

          default:
            l << name << " ";
            break;
        }
      }
      file << line.str() << "\n";
      l << endl;
      file << l.str();
    } // end of preparation


    //
    // simply write all data
    //

    ofstream& file = *plotfile;

    file << _goal;

    map<ID, vector<double> >::iterator it = data.begin();
    map<ID, vector<double> >::iterator end = data.end();
    for ( ; it != end; ++it)
    {
      const vector<double>& values = it->second;
      for (size_t i = 0; i < values.size(); ++i)
        file << " " << values[i];
    }

    file << endl;

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

  // the filename suffix
  _suffix = TiberCad::get_filename_suffix();

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
    _goal = values[i];

    {
      ostringstream os;
      os << "Ramping to sweep value " << _variable << " = " << _goal;
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
    // we strip the leading '$' for the filename
    ostringstream suffix;
    suffix << _variable.substr(1, string::npos) << "_" << _goal;
    TiberCad::prepend_to_filename_suffix(suffix.str());

    bool failed = false;
    try
    {
      ramp.ramp(_goal);

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
        write_global_data(*_simulations[j], plotfiles[j]);

      }
    }
    catch (SolveFailedException& e)
    {
      failed = true;
      Messages::error(e.what());
      if (i == 0)
        throw SolveFailedException("Already first sweep step could not be solved.");
    }

    TiberCad::drop_first_filename_suffix();

    if (failed)
    {
      ostringstream os;
      os << "Sweep failed at " << _variable << " = " << _goal;
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
Sweep::get_solution_secure(map<ID, vector<double> >& values)
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


