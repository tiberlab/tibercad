// $Id$

#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Control.h"
#include "DLLoader.h"

#include "DriftDiffusion.h"
#include "ExcitonTransport.h"
#include "Macrostrain.h"
#include "EnvelopFunctionApprox.h"
#include "MacroHeatBalance.h"
#include "QuantumDensity.h"
#include "OpticsKP.h"
#include "QuantumDispersion.h"
#include "TightBinding.h"
#include "OptRecombinSpectrum.h"
#include "Poisson.h"
#include "TunnelingCurrent.h"

#include "Sweep.h"
#include "SelfconsistentSolver.h"

#include "Utils.h"

#include "GMVIO_cell.h"
#include "tecplot_IO_cell.h"
#include "gnuplot_io.h"
#include "GraceIO.h"

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
    _equilibrium_is_solved(false),
    _relaxation_factor(1.0)
{
  ID new_id = _simulation_map.size() + 1;
  _id = new_id;

  _simulation_map[new_id] = this;
}


SimulationInterface::~SimulationInterface(void)
{
  map<ID, NumericVector<double>*>::iterator it(_remembered_solutions.begin());
  map<ID, NumericVector<double>*>::iterator end(_remembered_solutions.end());
  for ( ; it != end; ++it)
    delete it->second;
}


SimulationInterface*
SimulationInterface::create(const string& type,
        const ModelOptions& options)
{
  SimulationInterface* sim = NULL;

  // First we attempt to open a shared library
  //
  DLLoader::LibraryInterface iface;
  bool success = DLLoader::open_library(type, iface);

  create_t create_fnc = (create_t) iface.create_fnc;
  destroy_t destroy_fnc = (destroy_t) iface.destroy_fnc;

  if(success)
    sim = create_fnc();
  else
  {
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
    else if (type == "thermal")
      sim = MacroHeatBalance::create();
    else if (type == "selfconsistent")
      sim = SelfconsistentSolver::create();
    else if (type == "quantumdensity")
      sim = QuantumDensity::create();
    else if (type == "opticskp")
      sim = OpticsKP::create();
    else if (type == "quantumdispersion")
      sim = QuantumDispersion::create();
    else if (type == "tunnelingcurrent")
      sim = TunnelingCurrent::create();
    else if (type == "tightbinding")
      sim = TightBinding::create();
    else if (type == "opticalspectrum")
      sim = OptRecombinSpectrum::create();
    else if (type == "poisson")
      sim = Poisson::create();
}

  if (sim != NULL)
  {
    sim->_libhandle = iface.handle;
    sim->_create = create_fnc;
    sim->_destroy = destroy_fnc;


    sim->set_options(options);

    // we let it know what's its identifier
    sim->set_type(type);

    //! set the name
    string defaultname = Utils::extract_typename(typeid(*sim));
    sim->_name = sim->get_options().get_option("name", defaultname);
    sim->_options.delete_option("name");
    
    sim->_relaxation_factor =
      sim->get_options().get_option("relaxation_factor", sim->_relaxation_factor);
    sim->_options.delete_option("relaxation_factor");

#ifdef DEBUG
    cout << "Added simulator" << endl;
    cout << "        ID   = " << sim->get_id() << endl;
    cout << "        type = " << sim->get_type() << endl;
    cout << "        name = " << sim->get_name() <<
      " / default name = " << sim->get_default_name() << endl;
    cout << "        address = " << sim << endl;
#endif
  }

  return sim;
}




void
SimulationInterface::destroy(SimulationInterface* p)
{
  if (p != NULL)
  {
#ifdef DEBUG
    cerr << "Deleted simulator (ID = " << p->get_id() <<
      " name = " << p->get_name() << " / type_id = " <<
      p->get_default_name() << ")";
    cerr << " address = " << p << endl;
#endif

    libhandle_t libhandle = p->_libhandle;
    destroy_t destroy_fnc = p->_destroy;

    if (destroy_fnc != NULL)
      destroy_fnc(p);
    else
      delete p;

    if (libhandle != NULL)
      DLLoader::close_library(libhandle);
  }
}




void
SimulationInterface::init(void) throw (InitFailedException)
{
  if (!_is_initialized)
  {
#ifdef DEBUG
    cerr << "Initialize " << get_name() << "... ";
#endif

    // build name for equation systems
    create_equation_system_name();

    if (_environment != NULL)
    {
      _environment->prepare_for_solve();
      _scaling.set_calc_mesh_units((_environment->get_device()).get_mesh_units());
    }

    do_init();
    
  }

  _is_initialized = true;

#ifdef DEBUG
    cerr << "done" << endl;
#endif
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
SimulationInterface::solve_equilibrium(void) throw (SolveFailedException) 
{
 
  if (!_equilibrium_is_solved)
  {

    PerfLog perflog(get_name() + ": solve_equilibrium", false);
    perflog.start_event("solve_equilibrium");

    assert(is_initialized());

    if (_environment != NULL)
      _environment->prepare_for_solve();


    do_equilibrium();

    _is_solved = true;

    perflog.stop_event("solve_equilibrium");

    _equilibrium_is_solved = true;
  }
}




void
SimulationInterface::solve(void) throw (SolveFailedException) 
{
 

  PerfLog perflog(get_name() + ": solve", false);
  perflog.start_event("solve");

  assert(is_initialized());

  if (_environment != NULL)
    _environment->prepare_for_solve();


  try
  {
    do_solve();
  }
  catch (PetscRuntimeError& e)
  {
    ostringstream s;
    s << "Control: Solve of " << get_name() << " failed." << endl <<
      "    Cause: " << e.what() << " : " << e.get_reason();
    throw SolveFailedException(s.str());
  }
  catch (runtime_error& e)
  {
    ostringstream s;
    s << "Control: Solve of " << get_name() << " failed." << endl <<
      "    Cause: " << e.what();
    throw SolveFailedException(s.str());
  }
  catch (...)
  {
    ostringstream s;
    s << "Control: Solve of " << get_name() << " failed." << endl <<
      "    Cause: Unknown";
    throw SolveFailedException(s.str());
  }

  _is_solved = true;
  
  perflog.stop_event("solve");
}




NumericVector<double>&
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
SimulationInterface::get_integrated_quantities(
    const std::set<std::string>& variables, std::vector<double>& values)
{

  if (_environment != NULL)
    get_environment().prepare_for_solve();

  values.resize(0);
  build_integrated_quantities(variables, values);
}



void
SimulationInterface::plot(void)
{
  if (_environment != NULL)
    get_environment().prepare_for_solve();

  do_plot();
}
  


void
SimulationInterface::do_plot(void)
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
  else if (format == "grace")
    suff = ".dat";

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
    else if (format == "grace")
      GraceIO(dev.get_mesh()).write_nodal_data(filename, results, names);
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
    else if (format == "grace")
      GraceIO(dev.get_mesh()).write_elemental_data(filename, results, names);
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




ID
SimulationInterface::do_remember_current_solution(ID id)
{
  map<ID, NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, NumericVector<double>*>::iterator it(_remembered_solutions.find(id));

  if (it != end)
    *(it->second) = get_solution_vector();
  else
  {
    if (_remembered_solutions.begin() == end)
      id = 1;
    else
      id = (--end)->first + 1;

    _remembered_solutions[id] = get_solution_vector().clone().release();
  }


  return id;
}


void
SimulationInterface::do_set_to_remembered_solution(ID id)
{
  map<ID, NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, NumericVector<double>*>::iterator it(_remembered_solutions.find(id));

  if (it != end)
    get_solution_vector() = *(it->second);
}



void
SimulationInterface::do_delete_remembered_solution(ID id)
{
  map<ID, NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, NumericVector<double>*>::iterator it(_remembered_solutions.find(id));
  if (it != end)
  {
    delete it->second;
    _remembered_solutions.erase(it);
  }
}



NumericVector<double>*
SimulationInterface::get_remembered_solution(ID id)
{
  NumericVector<double>* vec = NULL;

  map<ID, NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, NumericVector<double>*>::iterator it(_remembered_solutions.find(id));
  if (it != end)
    vec = it->second;

  return vec;
}


double
SimulationInterface::do_maximum_norm_of_difference(ID id)
{
  double norm = 0.0;

  map<ID, NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, NumericVector<double>*>::iterator it(_remembered_solutions.find(id));
  if (it != end)
  {
    NumericVector<double>& old = *(it->second);
    NumericVector<double>& current = get_solution_vector();

    assert(old.size() == current.size());

    unsigned int n = old.size();
    for (unsigned int i = 0; i < n; i++)
    {
      double d = fabs(current(i) - old(i));
      norm = (d > norm) ? d : norm;
    }
  }

  return norm;
}



AutoPtr<FEBase>
SimulationInterface::build_finite_element(unsigned int dim, FEType type,
    bool scale)
{
  assert(type.family == libMeshEnums::LAGRANGE);

  double x0 = scale ? get_scaling().get_length_scaling() : 1.0;
  double mu = get_scaling().get_calc_mesh_units();
  TiberCad::Symmetry sym = get_environment().get_device().get_symmetry();

  FEBase* fe;

  switch (dim)
  {
    case 1:
      {
        FiniteElement<1, libMeshEnums::LAGRANGE>* fem =
          new FiniteElement<1, libMeshEnums::LAGRANGE>(type);
        fem->set_symmetry(sym);
        fem->set_scaling(x0, mu);
        fe = fem;
      }
      break;
    case 2:
      {
        FiniteElement<2, libMeshEnums::LAGRANGE>* fem =
          new FiniteElement<2, libMeshEnums::LAGRANGE>(type);
        fem->set_symmetry(sym);
        fem->set_scaling(x0, mu);
        fe = fem;
      }
      break;
    case 3:
      {
        FiniteElement<3, libMeshEnums::LAGRANGE>* fem =
          new FiniteElement<3, libMeshEnums::LAGRANGE>(type);
        fem->set_symmetry(sym);
        fem->set_scaling(x0, mu);
        fe = fem;
      }
      break;
    default:
      fe = NULL;
  }
  assert(fe != NULL);
  
  return AutoPtr<FEBase>(fe);
}



void
SimulationInterface::get_elemental_results(
    const std::set<std::string>& variables,
    std::vector<double>& results, std::vector<std::string>& legend)
{
  results.resize(0);
  legend.resize(0);
  if (variables.size() > 0)
  {
    build_elemental_results(variables, results, legend);

    unsigned int n = get_environment().get_mesh().n_active_elem();

   

    if (results.size() != n * legend.size())
    {
      ostringstream s;
      s << "SimulationInterface::get_elemental_results: simulation "
        << this->get_name() << " returned broken elemental results.";
      throw runtime_error(s.str()); 
    }
  }
}


void
SimulationInterface::get_nodal_results(
    const std::set<std::string>& variables,
    std::vector<double>& results, std::vector<std::string>& legend)
{
  results.resize(0);
  legend.resize(0);
  if (variables.size() > 0)
  {
    build_nodal_results(variables, results, legend);

    unsigned int n = get_environment().get_mesh().n_nodes();
    if (results.size() != n * legend.size())
    {
      ostringstream s;
      s << "SimulationInterface::get_nodal_results: simulation "
        << this->get_name() << " returned broken nodal results.";
      throw runtime_error(s.str()); 
    }
  }
}


void
SimulationInterface::get_integrated_quantities_description(
    const std::set<std::string>& variables,
    std::vector<std::string>& legend,
    std::vector<std::string>& description)
{
  legend.resize(0);
  description.resize(0);
  if (variables.size() > 0)
    build_integrated_quantities_description(variables, legend, description);
}



