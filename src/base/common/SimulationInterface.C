// $Id$

#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Control.h"
#include "DLLoader.h"
#include "Material.h"
#include "Embracing.h"

#ifndef BUILD_TIBER_MODULES
#include "DriftDiffusion.h"
#include "ExcitonTransport.h"
#include "DSSC.h"
#endif
#include "Macrostrain.h"
#include "EnvelopFunctionApprox.h"
#include "MacroHeatBalance.h"
#include "QuantumDensity.h"
#include "OpticsKP.h"
#include "QuantumDispersion.h"
#include "Dftb.h"
#include "OptRecombinSpectrum.h"
#include "Poisson.h"
#include "TunnelingCurrent.h"
#include "MaxwellEquations.h"
#include "CrackStrain.h"


#include "Sweep.h"
#include "RelaxationMethod.h"
#include "ModifiedBroyden.h"
#include "Utils.h"

#include "DataOutput.h"


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
    _has_solution_vector(true),
    _verbosity(1)
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

  EmbracingMap::iterator embit(_embracings.begin());
  const EmbracingMap::iterator embend(_embracings.end());
  for ( ; embit != embend; ++embit)
    delete embit->second;

}


SimulationInterface*
SimulationInterface::create(const string& type,
                            const ModelOptions& options)
{
  SimulationInterface* sim = NULL;

  string type_name(type);
  string flavour = options.get_option("flavour", "");
  if (flavour.size() != 0)
    type_name += "_" + flavour;

  // First we attempt to open a shared library
  //
  DLLoader::LibraryInterface iface;
  bool success = DLLoader::open_library(type_name, iface);

  create_t create_fnc = (create_t) iface.create_fnc;
  destroy_t destroy_fnc = (destroy_t) iface.destroy_fnc;

  if(success)
    sim = create_fnc();
  else
  {
#ifndef BUILD_TIBER_MODULES
    if (type_name == "driftdiffusion")
      sim = DriftDiffusion::create();
    else if (type_name == "dssc")
      sim = DSSC::create();
    else if (type_name == "excitontransport")
      sim = ExcitonTransport::create();
    else if (type_name == "macrostrain")
#else
    if (type_name == "macrostrain")
#endif
      sim = Macrostrain::create();
    else if (type_name == "crackstrain")
      sim = CrackStrain::create();
    else if (type_name == "efaschroedinger")
      sim = EnvelopFunctionApprox::create();
    else if (type_name == "sweep")
      sim = Sweep::create();
    else if (type_name == "thermal")
      sim = MacroHeatBalance::create();
    else if (type_name == "selfconsistent")
      sim = RelaxationMethod::create();
    else if (type_name == "selfconsistent_relaxation")
      sim = RelaxationMethod::create();
    else if (type_name == "selfconsistent_broyden")
      sim = ModifiedBroyden::create();
    else if (type_name == "quantumdensity")
      sim = QuantumDensity::create();
    else if (type_name == "opticskp")
      sim = OpticsKP::create();
    else if (type_name == "quantumdispersion")
      sim = QuantumDispersion::create();
    else if (type_name == "tunnelingcurrent")
      sim = TunnelingCurrent::create();
#ifdef ENABLE_DFTB
    else if (type_name == "tightbinding")
      sim = Dftb::create();
#endif
    else if (type_name == "opticalspectrum")
      sim = OptRecombinSpectrum::create();
    else if (type_name == "poisson")
      sim = Poisson::create();
    else if (type_name == "maxwell")
      sim = MaxwellEquations::create();

  }

  if (sim != NULL)
  {
    sim->_libhandle = iface.handle;
    sim->_create = create_fnc;
    sim->_destroy = destroy_fnc;

    sim->set_options(options);

    // we let it know what's its identifier
    sim->set_type(type_name);

    // set the name
    // we use the type name as found in the input file as default name
    //string defaultname = Utils::extract_typename(typeid(*sim));
    string defaultname(type);
    sim->_name = sim->get_options().get_option("name", defaultname);
    sim->_options.delete_option("name");


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
SimulationInterface::do_print_info(void)
{
  cout << "  no information available";
}



bool
SimulationInterface::includes_region(ID region_id) const
{
  return get_environment().contains_region(region_id);
}


void
SimulationInterface::init(void) throw (InitFailedException)
{
  if (!_is_initialized)
  {
#ifdef DEBUG
    cerr << "Initialize " << get_name() << "... " << endl;
#endif

    // build name for equation systems
    create_equation_system_name();

    if (_environment != NULL)
    {
      _environment->prepare_for_solve();
      _scaling.set_calc_mesh_units((_environment->get_device()).get_mesh_units());
    }

    
    _verbosity = get_options().get_option("verbose", _verbosity);
    do_init();
    
  }

  _is_initialized = true;

#ifdef DEBUG
  cerr << "init of " << get_name() << " done" << endl;
#endif

  if (verbose() > 0)
  {
    cout << endl << 
      ">>================================================================<<"
      << endl << "Simulation options for " << get_name() << " (" <<
      get_default_name() << ")"
#ifdef DEBUG
      << " ptr = " << this
#endif
      << endl << endl;
    do_print_info();
    cout << endl;
    
    set<string> names;
    get_environment().get_region_names(names);
    for (set<string>::const_iterator it(names.begin()); it != names.end(); ++it)
    {
      Device& dev = get_environment().get_device();
      vector<ID> ids;
      dev.get_region_ids(*it, ids);
      assert(ids.size() != 0);
      PhysicalModel* mod =
        dev.get_material(ids[0])->get_model(get_id());
      if (mod != NULL)
      {
        cout << "  ** Model details (region " << *it << ", " << 
          dev.get_material(ids[0])->get_name() << "):" << endl;
        mod->print_info();
        cout << endl;
      }
    }
    
    cout << endl <<
      ">>================================================================<<"
      << endl << endl;
  }
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
  //catch (PetscRuntimeError& e)
  //{
  //  ostringstream s;
  //  s << "Control: Solve of " << get_name() << " failed." << endl <<
  //    "    Cause: " << e.what() << " : " << e.get_reason();
  //  throw SolveFailedException(s.str());
  //}
  catch (runtime_error& e)
  {
    ostringstream s;
    s << get_name() << ": " << e.what();
    throw SolveFailedException(s.str());
  }
  catch (...)
  {
    ostringstream s;
    s << get_name() << ": unknown error occurred";
    throw SolveFailedException(s.str());
  }

  _is_solved = true;
  
  perflog.stop_event("solve");
}



NumericVector<double>&
SimulationInterface::do_get_solution_vector(void)
{
  const EquationSystems& eq = get_equation_systems();
  const System& sys = eq.get_system(get_equation_system_name());

  return *sys.solution;
}



void
SimulationInterface::do_set_solution_vector(
    const NumericVector<double>& new_solution)
{
  get_solution_vector() = new_solution;
}





BoundaryProperties*
SimulationInterface::create_boundary_model(const ModelOptions& options) const
  throw (ModelErrorException)
{
  ignore_unused_variable(options);

  return NULL;
}


      
PhysicalModel*
SimulationInterface::create_physical_model(const ModelOptions& options,
    const Material* mat) const
  throw (ModelErrorException)
{
  ignore_unused_variable(options);
  ignore_unused_variable(mat);
  
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

  vector<double> results;
  vector<string> names;

  DataOutput data_output(dev.get_mesh(), get_control().get_output_format());


  
  // 
  // nodal values
  //
  get_nodal_results(get_control().get_plotvariables(), results, names);
  if (names.size() > 0)
  {
    string filename(outdir + "/" + get_name() + "_nodal" + suffix);
    data_output.write_nodal_data(filename, results, names);
  }



  //
  // elemental values
  //
  get_elemental_results(get_control().get_plotvariables(), results, names);
  if (names.size() > 0)
  {
    string filename(outdir + "/" + get_name() + "_elemental" + suffix);
    data_output.write_cell_data(filename, results, names);
  }




  // materials  output: for each simulation an output file with IDs of
  // physical regions activated for that simulation
  // IDs are taken from the meshdata object associated to the  device;
  // IDs are data associated to elements

  std::vector<Number> translated_data;
  std::vector<std::string> data_names;

  (dev.get_meshdata())->activate(); 
  (dev.get_meshdata())->translate_elem_data(dev.get_mesh(),
                                            translated_data,
                                            data_names);

  if (data_names.size() > 0)
  {
    string filename(outdir + "/" + get_name() + "_materials" + suffix);
    data_output.write_cell_data(filename, translated_data, data_names);
  }




  // 
  // integrated properties
  //
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

      // if nn != nr, we print data in columns, otherwise on a row
      if (nn != nr)
      {
        // TODO is completely without logic
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
  if (!has_solution_vector()) return INVALID_ID;

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



double
SimulationInterface::do_l2_norm_of_difference(ID id)
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
      double d = current(i) - old(i);
      norm += d * d;
    }
  }

  return sqrt(norm);
}




void
SimulationInterface::do_scale_solution(double factor)
{
  if (has_solution_vector())
    get_solution_vector().scale(factor);
}
  



void
SimulationInterface::do_add_scaled_remembered_solution(ID id, double factor)
{
  map<ID, NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, NumericVector<double>*>::iterator it(_remembered_solutions.find(id));
  if (it != end)
  {
    get_solution_vector().add(factor, *(it->second));
  }
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
SimulationInterface::get_nodal_results(const std::set<std::string>& variables,
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






bool
SimulationInterface::get_solution(const Elem* elem, const set<ID>& ids,
                                  vector<map<ID, double> >& values)
{

  if ((ids.size() == 0) || (elem == NULL)) return false;
  if (!is_solved()) return false;
  
  SimulationEnvironment& env = get_environment();

  const Elem* el = elem;

  bool flag = true;

  if (!env.contains_element(elem))
  {
    // perhaps the parent is?
    const Elem* parent = elem->parent();

    while ((parent != NULL) && (!env.contains_element(parent)))
      parent = parent->parent();

    el = parent; // is NULL if no parent

    if (el != NULL)
    {
      int nn = elem->n_nodes();

      // the nodes are now inner points of the parent element
      vector<Point> p(nn);
      for (int i = 0; i < nn; i++)
        p[i] = elem->point(i);
      
      values.resize(nn);
      get_solution_secure(elem, p, ids, values);
    }
    else
    {
      // no active parent, so look for children
      vector<const Elem*> tree;
      elem->family_tree(tree, false);

      set<const Elem*> elem_list;
      unsigned int len = tree.size();
      for (unsigned int i = 0; i < len; i++)
      {
        const Elem* elem_i = tree[i];
        if (env.contains_element(elem_i))
          elem_list.insert(elem_i);
      }

      unsigned int np = elem->n_nodes();
      for (unsigned int i = 0; i < np; i++)
      {
        set<const Elem*>::iterator el_it = elem_list.begin();
        set<const Elem*>::iterator el_end = elem_list.end();
        for ( ; el_it != el_end; ++el_it)
        {
          el = *el_it;
          if (el->contains_point(elem->point(i)))
          {
            get_solution_secure(el, elem->point(i), ids, values[i]);
            // we have found it, so get out of the for loop
            break;
          }
        }
        if (el_it == el_end)
          flag = false;
      }
    }
  }
  else
  {
    values.resize(elem->n_nodes());
    get_solution_secure(elem, ids, values);
  }
}



bool
SimulationInterface::get_solution(const Elem* elem, const Point& p,
                                  const set<ID>& ids, map<ID, double>& values)
{
  vector<Point> points(1, p);
  vector<map<ID, double> > vals(1);

  bool flag = get_solution(elem, points, ids, vals);
 
  if (flag)
    values = vals[0];

  return flag;
}



bool
SimulationInterface::get_solution(const Elem* elem, const Point& p,
                                  ID id, double& value)
{
  vector<Point> points(1, p);
  set<ID> ids;
  ids.insert(id);
  vector<map<ID, double> > vals(1);

  bool flag = get_solution(elem, points, ids, vals);

  if (flag)
    value = vals[0][id];

  return flag;
}



bool
SimulationInterface::get_solution(const Elem* elem, const vector<Point>& p,
                                  const set<ID>& ids,
                                  vector<map<ID, double> >& values)
{
  bool flag = true;

  unsigned int np = p.size();
  if ((np == 0) || (ids.size() == 0) || (elem == NULL)) return false;
  if (!is_solved()) return false;

  values.resize(np);

  // this will contain the element in which p lie
  const Elem* el = elem;

  SimulationEnvironment& env = get_environment();
  
  // check if elem is an active element of the simulation
  if (!env.contains_element(elem))
  {
    // do we have a parent element in the list?
    const Elem* parent = elem->parent();

    while ((parent != NULL) && (!env.contains_element(parent)))
      parent = parent->parent();

    el = parent; // is NULL if no parent
  }

  if (el != NULL) // we found it!
    get_solution_secure(el, p, ids, values);
  else
  {
    // no parent, so check for children
    vector<const Elem*> tree;
    elem->family_tree(tree, false);

    set<const Elem*> elem_list;
    unsigned int len = tree.size();
    for (unsigned int i = 0; i < len; i++)
    {
      const Elem* elem_i = tree[i];
      if (env.contains_element(elem_i))
        elem_list.insert(elem_i);
    }
    for (unsigned int i = 0; i < np; i++)
    {
      set<const Elem*>::iterator el_it = elem_list.begin();
      set<const Elem*>::iterator el_end = elem_list.end();
      for ( ; el_it != el_end; ++el_it)
      {
        el = *el_it;
        if (el->contains_point(p[i]))
        {
          get_solution_secure(el, p[i], ids, values[i]);
          // we have found it, so get out of the for loop
          break;
        }
      }
      if (el_it == el_end)
        flag = false;
    }
  }

  return flag;
}



bool
SimulationInterface::get_solution(const Elem* elem, ID id, vector<double>& values)
{
  set<ID> ids;
  ids.insert(id);
  values.resize(elem->n_nodes());

  vector<map<ID, double> > vals;
  
  bool flag = get_solution(elem, ids, vals);

  if (flag && (vals[0].size() != 0))
  {
    int n = vals.size();
    values.resize(n);
    for (int i = 0; i < n; i++)
    {
      values[i] = vals[i][id];
    }
  }

  return flag;
}



bool
SimulationInterface::get_solution(const Elem* elem, const vector<Point>& p,
                                  ID id, vector<double>& values)
{
  set<ID> ids;
  ids.insert(id);
  values.resize(p.size());

  vector<map<ID, double> > vals;
  
  bool flag = get_solution(elem, p, ids, vals);

  if (flag && (vals[0].size() != 0))
  {
    int n = vals.size();
    values.resize(n);
    for (int i = 0; i < n; i++)
    {
      values[i] = vals[i][id];
    }
  }

  return flag;
}


ModelOptions&
SimulationInterface::get_solver_options(void)
{
  if (!_options.has_submodel("$Solver"))
    _options.add_submodel("$Solver", ModelOptions());

  ModelOptions::submodel_iterator it(_options.submodels_begin("$Solver"));

  return it->second;
}
                   


Embracing*
SimulationInterface::create_embracing_region(
    SimulationInterface* other_simulation, const ModelOptions& options)
{
  Embracing* emb = NULL;
  if (other_simulation != NULL)
    if (_embracings.find(other_simulation) != _embracings.end())
      emb = _embracings[other_simulation];
    else
    {
      emb = new Embracing(this, other_simulation);
      _embracings[other_simulation] = emb;
      emb->init(options);
    }

  return emb;
}

