// $Id$

#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Control.h"

#include "DriftDiffusion.h"
#include "ExcitonTransport.h"
#include "Macrostrain.h"
#include "EnvelopFunctionApprox.h"
#include "Sweep.h"


#include "Utils.h"

#include "GMVIO_cell.h"
#include "tecplot_IO_cell.h"
#include "gnuplot_io.h"

// LibMesh includes
#include "system.h"

#include <iostream>

using namespace std;

SimulationInterface::SimulationMap
SimulationInterface::_simulation_map;



SimulationInterface::SimulationInterface(void)
  : _environment(0),
    _is_initialized(false),
    _is_solved(false),
    _relaxation_factor(1.0)
{
  ID new_id = _simulation_map.size() + 1;
  _id = new_id;

  _simulation_map[new_id] = this;
}




SimulationInterface*
SimulationInterface::create(const string& type,
        const ModelOptions& options)
{
  SimulationInterface* sim = NULL;

  if (type == "driftdiffusion")
    sim = DriftDiffusion::create();
  else if (type == "excitontransport")
    sim = ExcitonTransport::create();
  else if (type == "macrostrain")
    sim = Macrostrain::create();
  else if (type == "efaschroedinger")
    sim = EnvelopFunctionApprox::create();
  else if (type == "sweep")
    sim = Sweep::create();

  if (sim != NULL)
  {
    sim->set_options(options);

    // we let it know what's its identifier
    sim->set_type(type);

    //! set the name
    string defaultname = Utils::extract_typename(typeid(*sim));
    sim->_name = sim->get_options().get_option("name", defaultname);
    sim->_options.delete_option("name");

#ifdef DEBUG
    cout << "Added simulator" << endl;
    cout << "        ID   = " << sim->get_id() << endl;
    cout << "        type = " << sim->get_type() << endl;
    cout << "        name = " << sim->get_name() <<
      " / default name = " << sim->get_default_name() << endl;
#endif
  }

  return sim;
}




void
SimulationInterface::destroy(SimulationInterface* p)
{
  // TODO better call a module internal method
  delete p;
}




void
SimulationInterface::init(void) throw (InitFailedException)
{
  if (!_is_initialized)
  {
    // build name for equation systems
    create_equation_system_name();

    if (_environment != NULL)
      _environment->prepare_for_solve();
    
    do_init();
    
  }

  _is_initialized = true;
}




void
SimulationInterface::create_equation_system_name(void)
{
  ostringstream o;
  o << get_name() << get_id();
  _eq_system_name = o.str();
}


string
SimulationInterface::get_default_name(void) const
{
  return Utils::extract_typename(typeid(*this));
}



SimulationInterface*
SimulationInterface::find_simulation(const string& name)
{
  SimulationInterface* sim = NULL;

  SimulationMap::iterator it(_simulation_map.begin());
  SimulationMap::iterator end(_simulation_map.end());

  if (it != end)
  {
    if (name == "") // we just take the first we can find ...
      sim = it->second;
    else
    {
      // look for user defined names
      for ( ; (it != end) && ((it->second)->get_name() != name); ++it);

      if (it != end)
        sim = it->second;

      if (it == end)
      {
        // name could be model identifier
        it = _simulation_map.begin();
        for ( ; (it != end) && ((it->second)->get_type() != name); ++it);

        if (it != end)
          sim = it->second;
      }

      if (it == end)
      {
        // we even look for the default name
        it = _simulation_map.begin();
        for ( ; (it != end) && ((it->second)->get_default_name() != name); ++it);

        if (it != end)
          sim = it->second;
      }

    }
  }

  return sim;
}




EquationSystems&
SimulationInterface::get_equation_systems(void) const
{
  return _environment->get_device().get_equation_systems();
}




void
SimulationInterface::solve(void) throw (SolveFailedException) 
{
 

  PerfLog perflog(get_name() + ": solve", false);
  perflog.start_event("solve");

  assert(_is_initialized);

  if (_environment != NULL)
    _environment->prepare_for_solve();

  do_solve();

  _is_solved = true;
  
  perflog.stop_event("solve");
}




NumericVector<Real>&
SimulationInterface::get_solution_vector(void)
{
  const EquationSystems& eq = get_equation_systems();
  const System& sys = eq.get_system(get_equation_system_name());

  return *sys.solution;
}



BoundaryProperties*
SimulationInterface::create_boundary_model(const ModelOptions& options) const
      throw (ModelErrorException)
{
  ignore_unused_variable(options);

  return NULL;
}


      
PhysicalModel*
SimulationInterface::create_physical_model(const ModelOptions& options) const
      throw (ModelErrorException)
{
  ignore_unused_variable(options);
  
  return NULL;
}



void
SimulationInterface::plot(void)
{
  const Device& dev = get_environment().get_device();

  string suffix = get_control().get_filename_suffix();
  string outdir = get_control().get_output_dir();
  string format = get_control().get_output_format();

  string suff;
  if (format == "gmv")
    suff = ".gmv";
  else if (format == "ise")
    suff = ".plt";

  vector<double> results;
  vector<string> names;

  // nodal values
  get_nodal_results(get_control().get_plotvariables(), results, names);
  if (names.size() > 0)
  {
    string filename(outdir + "/" + get_name() +
        "_nodal" + suffix + suff);

    if (format == "gmv")
      GMVIO(dev.get_mesh()).write_nodal_data(filename, results, names);
    else if (format == "gnuplot")
      GnuPlotIO(dev.get_mesh()).write_nodal_data(filename, results, names);
    else if (format == "ise")
      TecplotIO(dev.get_mesh()).write_nodal_data(filename, results, names);
    else
    {
      cout << "Output format not supported. Falling back to GMV." << endl;
      GMVIO(dev.get_mesh()).write_nodal_data(filename, results, names);
    }
  }

  // elemental values
  get_elemental_results(get_control().get_plotvariables(), results, names);
  if (names.size() > 0)
  {
    string filename(outdir + "/" + get_name() +
        "_elemental" + suffix + suff);

    if (format == "gmv")
      GMVIO_cell(dev.get_mesh()).write_ascii_cell_data(filename, results, names);
    else if (format == "gnuplot")
      cout << "GnuPlot does not currently support cell data." << endl;
    else if (format == "ise")
      TecplotIO_cell(dev.get_mesh()).write_cell_data(filename, results, names);
    else
    {
      cout << "Output format not supported. Falling back to GMV." << endl;
      GMVIO_cell(dev.get_mesh()).write_ascii_cell_data(filename, results, names);
    }
  }

  // integrated properties
  vector<string> description;
  get_integrated_quantities_description(get_control().get_plotvariables(),
      names, description);
  if (names.size() > 0)
  {
    string filename(outdir + "/" + get_name() + suffix + ".dat");
    ofstream file;
    file.open(filename.c_str());
    if (file.good())
    {
      // header
      file << "# Simulation: " << get_name() << endl;
      file << "# Data:" << endl;
      for (unsigned int i = 0; i < description.size(); i++)
        file << "#    * " << description[i] << endl;
      file << "#" << endl;
      
      build_integrated_quantities(get_control().get_plotvariables(), results);
      
      unsigned int nn = names.size();
      unsigned int nr = results.size();

      // if nn == nr, we print data in columns, otherwise on a row
      if (nn == nr)
      {
        ostringstream l;
        l << setprecision(12);
        for (unsigned int i = 0; i < nn; i++)
          l << names[i] << "   " << results[i] << endl;

        file << l.str();
      }
      else
      {
        // legend
        ostringstream l;
        l << setprecision(12);
        l << "# ";
        for (unsigned int i = 0; i < nn; i++)
          l << names[i] << "   ";
        l << endl;

        // data
        for (unsigned int i = 0; i < nr; i++)
          l << results[i] << "   ";
        l << endl;
        file << l.str();
      }

      file.close();
    }
  }

}
