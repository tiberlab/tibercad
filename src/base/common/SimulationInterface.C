// $Id$

#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "MaterialBoundary.h"
#include "EdgeObject.h"
#include "NodeObject.h"
#include "Alloy.h"
#include "Embracing.h"

#include "Macrostrain.h"
#include "EnvelopFunctionApprox.h"
#include "QuantumDensity.h"
#include "OpticsKP.h"
#include "QuantumDispersion.h"
#include "Dftb.h"
#include "EmpiricalTightBinding.h"
#include "OpticsTB.h"
#include "OptRecombinSpectrum.h"
#include "Poisson.h"
#include "TunnelingCurrent.h"
#include "MaxwellEquations.h"
#include "PhononDispersion.h"
#include "CrackStrain.h"
//#include "MechanicalStructure.h"
#include "Sweep.h"
#include "RelaxationMethod.h"
#include "ModifiedBroyden.h"
#include "Utils.h"
#include "DataOutput.h"
#include "Messages.h"


// LibMesh includes
#include "system.h"

#include <sstream>



using namespace std;

SimulationInterface::SimulationMap
SimulationInterface::_simulation_map;

SolutionDescriptor
SimulationInterface::_invalid_descr = SolutionDescriptor();


SimulationInterface::SimulationInterface(const ModelOptions& options)
  : TiberModelObject(options),
    _environment(0),
    _is_initialized(false),
    _is_solved(false),
    _equilibrium_is_solved(false),
    _has_solution_vector(true),
    _verbosity(1)
{
  ID new_id = _simulation_map.size() + 1;
  _id = new_id;

  _relaxation = 1.0;

  // register in the map of simulations
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

  SimulationEnvironment::destroy(_environment);

  _simulation_map.erase(get_id());

  ostringstream os;
  os << "Deleted simulator (ID = " << get_id() <<
    " name = " << get_name() << " / type = " <<
    get_type() << ")" << " address = " << this;
  Messages::debug(os.str());

}


void
SimulationInterface::destroy(SimulationInterface* sim)
{
  SimulationMap::iterator it(_simulation_map.begin());
  for ( ; (it != _simulation_map.end()) && (it->second != sim); ++it);

  if (it != _simulation_map.end())
    delete it->second;

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

  if (type_name == "macrostrain")
    sim = Macrostrain::create(options);
  else if (type_name == "crackstrain")
    sim = CrackStrain::create(options);
  else if (type_name == "efaschroedinger")
    sim = EnvelopFunctionApprox::create(options);
  else if (type_name == "sweep")
    sim = Sweep::create(options);
 
  else if (type_name == "selfconsistent")
    sim = RelaxationMethod::create(options);
  else if (type_name == "selfconsistent_relaxation")
    sim = RelaxationMethod::create(options);
  else if (type_name == "selfconsistent_broyden")
    sim = ModifiedBroyden::create(options);
  else if (type_name == "quantumdensity")
    sim = QuantumDensity::create(options);
  else if (type_name == "opticskp")
    sim = OpticsKP::create(options);
  else if (type_name == "quantumdispersion")
    sim = QuantumDispersion::create(options);
  else if (type_name == "tunnelingcurrent")
    sim = TunnelingCurrent::create(options);
#ifdef ENABLE_DFTB
  else if (type_name == "densityfunctional_tb")
    sim = Dftb::create(options);
#endif
#ifdef ENABLE_UPTIGHT
  else if (type_name == "empirical_tb")
    sim = ETB::create(options);
  else if (type_name == "opticstb")
    sim = OpticsTB::create(options);
#endif
  else if (type_name == "opticalspectrum")
    sim = OptRecombinSpectrum::create(options);
  else if (type_name == "poisson")
    sim = Poisson::create(options);
  else if (type_name == "maxwell")
    sim = MaxwellEquations::create(options);
  else if (type_name == "phonondispersion")
    sim = PhononDispersion::create(options);
 // else if (type_name == "mechanical_structure")
  //  sim = MechanicalStructure::create(options);


  if (sim == NULL)
  {
    // try first in a module directory
    if ((type.size() > 0 ) && ((sim = create_from_library<SimulationInterface>(
           type + "/" + type_name, options)) == 0))
    {
        std::cout<<"ENTER"<<std::endl;
      sim = create_from_library<SimulationInterface>(type_name, options);
        std::cout<<sim<<std::endl;
    }
  }

  if (sim != NULL)
  {

    // we let it know what's its identifier
    sim->set_type(type);

    // set the name
    // we use the type name as found in the input file as default name
    string defaultname(type);
    sim->set_name(sim->get_options().get_option("name", defaultname));
    sim->get_options().delete_option("name");


    ostringstream os;
    os << "Added simulator" << Messages::endl;
    os << "        ID   = " << sim->get_id() << Messages::endl;
    os << "        type = " << sim->get_type() << Messages::endl;
    os << "        name = " << sim->get_name() <<
      " / default name = " << type << Messages::endl;
    os << "        address = " << sim << Messages::endl;
    Messages::debug(os.str());
  }

  return sim;
}



void
SimulationInterface::setup_solution_variables(void)
{
  // setup the set of solution variables to be plotted
  vector<string> plotvars;
  get_option("plot", plotvars);
  get_options().delete_option("plot");
  for (unsigned int i = 0; i < plotvars.size(); i++)
    _plotvariables.insert(plotvars[i]);

  // add common plot variables
  declare_solution(RegionIDs, REAL, CELL);
  add_alias("materials", RegionIDs);

  do_setup_solution_variables();
}


void
SimulationInterface::do_print_info(void)
{
}



bool
SimulationInterface::includes_region(ID region_id) const
{
  return get_environment().contains_region(region_id);
}



PhysicalModel*
SimulationInterface::new_bulk_model(const ModelOptions& options,
    const Material* material)
{
  PhysicalModel* pm = create_bulk_model(options, material);

  if (pm != NULL)
    _physical_models.insert(pm);

  return pm;
}


BoundaryProperties*
SimulationInterface::new_boundary_model(const ModelOptions& options)
{
  BoundaryProperties* bp = create_boundary_model(options);

  return bp;
}


PhysicalModel*
SimulationInterface::new_boundary_model(const ModelOptions& options,
    const Material* material_A, const Material* material_B)
{
  PhysicalModel* pm = create_boundary_model(options, material_A, material_B);

  if (pm != NULL)
    _boundary_models.insert(pm);

  return pm;
}


PhysicalModel*
SimulationInterface::new_edge_model(const ModelOptions& options)
{
  PhysicalModel* pm = create_edge_model(options);

  //if (pm != NULL)
  //  _edge_models.insert(pm);

  return pm;
}


PhysicalModel*
SimulationInterface::new_node_model(const ModelOptions& options)
{
  PhysicalModel* pm = create_node_model(options);

  //if (pm != NULL)
  //  _node_models.insert(pm);

  return pm;
}



PhysicalModel*
SimulationInterface::get_physical_model(ID region_id) const
{
  PhysicalModel* mod = NULL;
  Material* mat = get_environment().get_device().get_material(region_id);
  if (mat != NULL)
    mod = mat->get_model(get_id());

  return mod;
}



PhysicalModel*
SimulationInterface::_get_bulk_model(const Elem* elem) const
{
  PhysicalModel* mod = NULL;
  Material* mat = get_environment().get_device().get_material(elem);
  if (mat != NULL)
    mod = mat->get_model(get_id());

  return mod;
}



PhysicalModel*
SimulationInterface::_get_interface_model(const Elem* elem, int side) const
{
  PhysicalModel* mod = NULL;
  MaterialBoundary* mb =
      get_environment().get_device().get_boundary_object(elem, side);
  if (mb != NULL)
    mod = mb->get_model(get_id());

  return mod;
}



PhysicalModel*
SimulationInterface::_get_edge_model(const Elem* elem, int edge) const
{
  PhysicalModel* mod = NULL;
  EdgeObject* eo =
      get_environment().get_device().get_edge_object(elem, edge);
  if (eo != NULL)
    mod = eo->get_model(get_id());

  return mod;
}



PhysicalModel*
SimulationInterface::_get_node_model(const Elem* elem, int node) const
{
  PhysicalModel* mod = NULL;
  NodeObject* no =
      get_environment().get_device().get_node_object(elem, node);
  if (no != NULL)
    mod = no->get_model(get_id());

  return mod;
}



MeshBase&
SimulationInterface::get_mesh(void) const
{
  return get_environment().get_device().get_mesh();
}


double
SimulationInterface::get_mesh_units(void) const
{
  return get_environment().get_device().get_mesh_units();
}


void
SimulationInterface::init(void) throw (InitFailedException)
{
  if (!_is_initialized)
  {
    Messages::debug("Initialize " + get_name() + "... ");

    // build name for equation systems
    create_equation_system_name();

    if (_environment != NULL)
    {
      _environment->prepare_for_solve();
      _scaling.set_calc_mesh_units((_environment->get_device()).get_mesh_units());

      // plot the regions
      //plot_regions();
    }



    _verbosity = get_option("verbose", _verbosity);
    do_init();

  }

  _is_initialized = true;

  Messages::debug("init of " + get_name() + " done");

  if (verbose() > 0)
  {
    Messages m;
    m.newline();
    {
      ostringstream os;
      os << ">>>>";
      os.width(Messages::available_width() - 4);
      os.fill('>');
      os << "" << Messages::endl
          << "Simulation options for " << get_name() << " (" <<
          get_type() << ")";

      m.info(os.str());
    }
    m.indent();

    print_info();

    m.newline();

    set<string> names;
    get_environment().get_region_names(names);
    for (set<string>::const_iterator it(names.begin()); it != names.end(); ++it)
    {
      Device& dev = get_environment().get_device();
      vector<ID> ids;
      dev.get_active_region_ids(*it, ids);
      assert(ids.size() != 0);
      const Material* mat = dev.get_material(ids[0]);
      PhysicalModel* mod =
        mat->get_model(get_id());
      if (mod != NULL)
      {
        ostringstream os;
        os << "Region " << *it << ", " << mat->get_name();
        if (mat->is_alloy())
          os << " (x = " <<
            static_cast<const Alloy*>(mat)->get_molar_fraction() << ")";
        m.info(os.str());
        m.indent();
        mod->print_info();
        m.unindent();
      }
    }

    m.unindent();
    m.newline();
    {
      ostringstream os;
      os << "<<<<";
      os.width(Messages::available_width() - 4);
      os.fill('<');
      os << "" << Messages::endl;
      m.info(os.str());
    }
    m.newline();
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
      for ( ; (it != end) && ((it->second)->get_name() != name); ++it) {}

      if (it != end)
        sim = it->second;

      if (it == end)
      {
        // name could be model identifier
        it = _simulation_map.begin();
        for ( ; (it != end) && ((it->second)->get_type() != name); ++it) {}

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
SimulationInterface::add_plot_variable(const std::string& name)
{
  _plotvariables.insert(name);
}


void
SimulationInterface::solve(void) throw (SolveFailedException)
{


  PerfLog perflog(get_name() + ": solve", false);
  perflog.start_event("solve");

  assert(is_initialized());

  Utils::Timer tt;

  Messages m;
  m.newline();
  {
    string s(">>>>");
    s.reserve(Messages::available_width());
    for (int i = 0; i < Messages::available_width() - 4; i++)
      s += "-";
    m.info(s);
    ostringstream os;
    os << get_type() << " (name: " << get_name() << ")" << endl;
    m.info(os.str());
    m.newline();
  }
  m.indent();

  if (_environment != NULL)
    _environment->prepare_for_solve();



  try
  {
    do_solve();
  }
  catch (runtime_error& e)
  {
    ostringstream os;
    os << "Solve time: " << tt.elapsed_string();
    Messages::newline();
    Messages::info(os.str());

    ostringstream s;
    s << get_name() << ": " << e.what();
    throw SolveFailedException(s.str());
  }
  catch (...)
  {
    ostringstream os;
    os << "Solve time: " << tt.elapsed_string();
    Messages::newline();
    Messages::info(os.str());

    ostringstream s;
    s << get_name() << ": unknown error occurred";
    throw SolveFailedException(s.str());
  }

  _is_solved = true;

  m.unindent();

  ostringstream os;
  os << "Solve time: " << tt.elapsed_string();
  Messages::newline();
  Messages::info(os.str());

  int w = Messages::available_width();
  string s("<<<<");
  s.reserve(w);
  w -= 8 + get_name().size();
  for (int i = 0; i < w / 2; i++)
    s += "-";
  s += " (";
  s += get_name();
  s += ") ";
  for (int i = 0; i < (w / 2 + w % 2); i++)
    s += "-";
  m.info(s);

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
SimulationInterface::create_boundary_model(const ModelOptions&) const
  throw (ModelErrorException)
{
  return NULL;
}



PhysicalModel*
SimulationInterface::create_physical_model(const ModelOptions&,
    const Material*) const
  throw (ModelErrorException)
{
  return NULL;
}


PhysicalModel*
SimulationInterface::create_bulk_model(const ModelOptions& options,
    const Material* material) const
{
  // call the old implementation for back compatibility
  return create_physical_model(options, material);
}



PhysicalModel*
SimulationInterface::create_boundary_model(const ModelOptions&,
    const Material*, const Material*) const
{
  return NULL;
}



PhysicalModel*
SimulationInterface::create_edge_model(const ModelOptions&) const
{
  return NULL;
}



PhysicalModel*
SimulationInterface::create_node_model(const ModelOptions&) const
{
  return NULL;
}



void
SimulationInterface::get_integrated_quantities(std::vector<double>& values)
{

  if (_environment != NULL)
    get_environment().prepare_for_solve();

  values.resize(0);
  build_integrated_quantities(values);
}



void
SimulationInterface::plot(void)
{
  if (_environment != NULL)
    get_environment().prepare_for_solve();

  do_plot();
}









void
SimulationInterface::get_output_format(vector<string>& formats) const
{
  // there is at least one format in the string
  get_option("output_format", formats);
  assert(formats.size() > 0);
}


string
SimulationInterface::get_output_directory(void) const
{
  return get_option("resultpath", ".");
}


string
SimulationInterface::get_output_filename_prefix(void) const
{
  if (has_option("output_prefix"))
    return get_option("output_prefix", get_name());
  else
    return get_name();
}


string
SimulationInterface::get_output_filename(void) const
{
  return get_output_filename_prefix() + TiberCad::get_filename_suffix();
}


bool
SimulationInterface::binary_output(void) const
{
  return get_option("binary_output", true);
}




void
SimulationInterface::do_plot(void)
{
  plot_meshdata();
  plot_atomisticdata();
  plot_globaldata();
}




void
SimulationInterface::plot_meshdata(void)
{

  const MeshBase& mesh = get_mesh();

  //MeshBase::const_element_iterator it = mesh.active_local_elements_begin();
  //const MeshBase::const_element_iterator end = mesh.active_local_elements_end();
  MeshBase::const_element_iterator it = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();

  // first gather subdomain infos
  // number of elements
  map<ID, size_t> n_elem;
  // node translation table and connectivity
  map<ID, map<unsigned int, pair<unsigned int, unsigned short> > > node_conn;
  {
    for ( ; it != end; ++it)
    {
      const Elem* elem = *it;

      ID subdomain = elem->subdomain_id();

      if (!n_elem.count(subdomain))
        n_elem[subdomain] = 1;
      else
        n_elem[subdomain]++;

      for (unsigned int n = 0; n < elem->n_nodes(); ++n)
      {
        if (node_conn[subdomain].count(elem->node(n)) == 0)
        {
          unsigned int nodeid = node_conn[subdomain].size();
          node_conn[subdomain][elem->node(n)] = make_pair(nodeid, 1);
        }
        else
          (node_conn[subdomain][elem->node(n)].second)++;
      }
    }
  }


  // the container for the solution values
  map<ID, vector<double> > solutions;

  // the container for the data to be handed over to the writer
  map<ID, map<SolutionDescriptor, vector<double> > > data;

  // we prepare them to contain all mesh localized solution variables
  IDSet::const_iterator varit(_plotvariable_ids.begin());
  for ( ; varit != _plotvariable_ids.end(); ++varit)
  {
    const SolutionDescriptor& descr = get_solution_descriptor(*varit);
    if (descr.on_mesh())
    {
      solutions.insert(make_pair(descr.id(), vector<double>(0)));
      unsigned int ncomp = descr.n_components();
      map<ID, size_t>::iterator it(n_elem.begin());
      for ( ; it != n_elem.end(); ++it)
      {
        size_t len = 0;
        switch (descr.location())
        {
          case SolutionDescriptor::NODES:
            len = ncomp * node_conn[it->first].size();
            data[it->first].insert(make_pair(descr, vector<double>(len)));
            break;
          case SolutionDescriptor::CELL:
            len = ncomp * it->second;
            data[it->first].insert(make_pair(descr, vector<double>(0)));
            data[it->first][descr].reserve(len);
            break;
          default:
            break;
        }
      }
    }
  }


  for (it = mesh.active_elements_begin(); it != end; ++it)
  {
    const Elem* elem = *it;

    ID subdomain = elem->subdomain_id();

    get_solution(elem, solutions);

    // put them into the right vectors
    map<SolutionDescriptor, vector<double> >::iterator
                                    dit(data[subdomain].begin());
    const map<SolutionDescriptor, vector<double> >::iterator
                                    dend(data[subdomain].end());
    for ( ; dit != dend; ++dit)
    {
      const SolutionDescriptor& descr = dit->first;
      vector<double>& vec = dit->second;
      vector<double>& sol = solutions[descr.id()];
      
      if (sol.size() == 0)
        continue;

      if (descr.id() == RegionIDs)
      {
        vec.push_back(subdomain);
        continue;
      }

      switch (descr.location())
      {
        case SolutionDescriptor::CELL:
          vec.insert(vec.end(), sol.begin(), sol.end());
          break;

        case SolutionDescriptor::NODES:
          for (unsigned int n = 0; n < elem->n_nodes(); n++)
          {
            unsigned int ncomp = descr.n_components();
            unsigned int index =
              ncomp * node_conn[subdomain][elem->node(n)].first;
            unsigned short w = node_conn[subdomain][elem->node(n)].second;
            for (unsigned int i = 0; i < ncomp; i++)
              vec[index + i] += sol[ncomp * n + i] / w;

          }
          break;

        default:
          break;
      }
    }
  }



  vector<string> formats;
  get_output_format(formats);
  for (unsigned int i = 0; i < formats.size(); i++)
  {
    DataOutput* writer = DataOutput::create(formats[i]);
    if ((writer != NULL) && (_solution_descriptors.size() > 0))
    {
      writer->set_output_directory(get_output_directory());
      writer->set_filename(get_output_filename());
      if (binary_output())
        writer->set_binary();
      else
        writer->set_ascii();

      writer->set_mesh(get_mesh());

      map<ID, map<SolutionDescriptor, vector<double> > >::iterator dit(data.begin());
      const map<ID, map<SolutionDescriptor, vector<double> > >::iterator dend(data.end());
      for ( ; dit != dend; ++dit)
        writer->set_data(dit->second, dit->first);

      writer->write();
    }
  }

  do_plot_old();

}



void
SimulationInterface::plot_atomisticdata(void)
{
}



void
SimulationInterface::plot_globaldata(void)
{
  map<ID, vector<double> > values;

  IDSet::iterator it(_plotvariable_ids.begin());
  const IDSet::iterator end(_plotvariable_ids.end());
  for ( ; it != end; ++it)
  {
    ID id = *it;
    if (get_solution_descriptor(id).location() == SolutionDescriptor::GLOBAL)
      values.insert(make_pair(id, vector<double>(0)));
  }

  get_solution_secure(values);
  map<ID, vector<double> >::iterator ii(values.begin());
  for ( ; ii != values.end(); ++ii)
    if (ii->second.size() == 0)
      values.erase(ii);

  if (values.size() > 0)
  {
    string outdir = get_output_directory();

    string filename(outdir + "/" + get_output_filename() + ".dat");
    ofstream file;
    file.open(filename.c_str());
    if (file.good())
    {
      // header
      int width[5] = {22, 15, 10};
      file << "# Global data for simulation: " << get_name() << endl;
      //file << "# Data:" << endl;
      file << "#" << endl;
      {
        ostringstream os;
        os << "# name";
        int w = width[0];
        os.width(w - os.tellp());
        os << "" << "units";
        w += width[1];
        os.width(w - os.tellp());
        os << "" << "type" << endl;
        file << os.str();
      }

      map<ID, vector<double> >::iterator it(values.begin());
      for ( ; it != values.end(); ++it)
      {
        const SolutionDescriptor& descr = get_solution_descriptor(it->first);
        ostringstream os;
        os << "# " << descr.name();
        int w = width[0];
        os.width(w - os.tellp());
        os << "" << descr.units();
        w += width[1];
        os.width(w - os.tellp());
        os << "";
        if (descr.type() == SolutionDescriptor::NTUPLE)
          os << descr.n_components() << "-tuple";
        else
          os << descr.type();

        file << os.str() << endl;
      }

      file << "#" << endl;
      file << "# ";
      for (it = values.begin(); it != values.end(); ++it)
      {
        const SolutionDescriptor& descr = get_solution_descriptor(it->first);
        file << descr.name() << "  ";
      }
      file << endl;

      for (it = values.begin(); it != values.end(); ++it)
      {
        const vector<double>& vals = it->second;
        for (unsigned int i = 0; i < vals.size(); i++)
          file << vals[i] << "  ";
      }

      file << endl;

      file.close();
    }
  }
  else // temporary only !
  {
    string outdir = get_output_directory();

    vector<double> results;
    vector<string> names;

    vector<string> description;
    get_integrated_quantities_description(names, description);
    if (names.size() > 0)
    {
      string filename(outdir + "/" + get_output_filename() + ".dat");
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

        build_integrated_quantities(results);

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
}



void
SimulationInterface::do_plot_old(void)
{

  const Device& dev = get_environment().get_device();

  string suffix = TiberCad::get_filename_suffix();
  string outdir = get_output_directory();

  vector<double> results;
  vector<string> names;

  vector<string> formats;
  get_output_format(formats);
  string formatstr("(");
  for (unsigned int i = 0; i < formats.size(); i++)
    formatstr += formats[i] + ",";
  formatstr += ")";

  DataOutput data_output(dev.get_mesh(), formatstr);
  data_output.set_output_directory(outdir);



  //
  // nodal values
  //
  get_nodal_results(results, names);
  if (names.size() > 0)
  {
    string filename(get_name() + "_nodal" + suffix);
    data_output.write_nodal_data(filename, results, names);
  }



  //
  // elemental values
  //
  get_elemental_results(results, names);
  if (names.size() > 0)
  {
    string filename(get_name() + "_elemental" + suffix);
    data_output.write_cell_data(filename, results, names);
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

    // the solution vector might be unclosed which produces an error
    // in debug mode
    get_solution_vector().close();
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


#define __create_finite_elem(d, f) \
  { \
    FiniteElement<d, libMeshEnums::f>* fem = \
          new FiniteElement<d, libMeshEnums::f>(type); \
    fem->set_symmetry(sym); \
    fem->set_scaling(x0, mu); \
    fe = fem; \
  }

#define __switch_family(d) \
  switch (type.family) \
  { \
    case libMeshEnums::LAGRANGE: \
      __create_finite_elem(d, LAGRANGE); \
      break; \
    case libMeshEnums::MONOMIAL: \
      __create_finite_elem(d, MONOMIAL); \
      break; \
  }

AutoPtr<FEBase>
SimulationInterface::build_finite_element(unsigned int dim, FEType type,
                                          bool scale)
{

  double x0 = scale ? get_scaling().get_length_scaling() : 1.0;
  double mu = get_scaling().get_calc_mesh_units();
  TiberCad::Symmetry sym = get_environment().get_device().get_symmetry();

  FEBase* fe;

  switch (dim)
  {
    case 1:
      __switch_family(1);
      break;

    case 2:
      __switch_family(2);
      break;

    case 3:
      __switch_family(3);
      break;

    default:
      fe = NULL;
  }
  assert(fe != NULL);

  return AutoPtr<FEBase>(fe);
}



void
SimulationInterface::get_elemental_results(std::vector<double>& results,
    std::vector<std::string>& legend)
{
  results.resize(0);
  legend.resize(0);
  if (_plotvariables.size() > 0)
  {
    build_elemental_results(_plotvariables, results, legend);

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
SimulationInterface::get_nodal_results(std::vector<double>& results,
    std::vector<std::string>& legend)
{
  results.resize(0);
  legend.resize(0);
  if (_plotvariables.size() > 0)
  {
    build_nodal_results(_plotvariables, results, legend);

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
    std::vector<std::string>& legend,
    std::vector<std::string>& description)
{
  legend.resize(0);
  description.resize(0);
  build_integrated_quantities_description(legend, description);
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

  return flag;
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





bool
SimulationInterface::get_solution(const Elem* elem,
    map<ID, vector<double> >& values,
    const vector<Point>& p, bool local_coord)
{
  if (!is_solved()) return false;

  vector<Point> points(p);
  unsigned int nn = points.size();
  if (nn == 0)
  {
    nn = elem->n_nodes();
    points.resize(nn);
    for (unsigned int i = 0; i < nn; i++)
    {
      points[i] = elem->local_node(elem->type(), i);
    }
  }
  else if (!local_coord)
  {
    FEInterface::inverse_map(get_mesh().mesh_dimension(), FEType(), elem, p, points);
  }

  SimulationEnvironment& env = get_environment();

  const Elem* el = elem;

  bool flag = true;

  // first resize all data vectors to the right size
  map<ID, vector<double> >::iterator it(values.begin());
  map<ID, vector<double> >::iterator end(values.end());
  for ( ; it != end; ++it)
  {
    const SolutionDescriptor& sd = get_solution_descriptor(it->first);
    assert(sd.id() != INVALID_ID);
    unsigned int n_comp = sd.n_components();
    switch (sd.location())
    {
      case SolutionDescriptor::CELL:
        values[it->first].resize(n_comp);
        break;

      case SolutionDescriptor::NODES:
        values[it->first].resize(nn * n_comp);
        break;

      default:
      {
        ostringstream os;
        os << "In get_solution(): solution variable \'" << sd.name()
                        << "\' seems not to be associated to the mesh.";
        throw ModelErrorException(os.str());
      }
    }
  }

  if (!env.contains_element(elem))
  {
    // perhaps the parent is?
    const Elem* parent = elem->parent();

    while ((parent != NULL) && (!env.contains_element(parent)))
      parent = parent->parent();

    el = parent; // is NULL if no parent

    if (el != NULL)
    {
      get_solution_secure(elem, values, points);
    }
    else
    {
      // no active parent, so look for children
      // TODO This part is not tested at all !!!
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

      if (elem_list.size() == 0)
        flag = false;
      else
      {

        // we need a copy of values
        map<ID, vector<double> > valcopy(values);

        // area of elem
        double elem_vol = elem->volume();

        set<unsigned int> nodenr;
        for (unsigned int i = 0; i < nn; i++)
          nodenr.insert(i);

        set<const Elem*>::iterator el_it = elem_list.begin();
        set<const Elem*>::iterator el_end = elem_list.end();
        for ( ; el_it != el_end; ++el_it)
        {
          el = *el_it;

          vector<Point> my_p;
          vector<unsigned int> pid;
          my_p.reserve(1);
          pid.reserve(1);

          set<unsigned int>::iterator nit(nodenr.begin());
          const set<unsigned int>::iterator nend(nodenr.end());
          for ( ; nit != nend; ++nit)
          {
            unsigned int node = *nit;
            if (FEInterface::on_reference_element(points[node], el->type()))
            {
              nodenr.erase(nit);
              my_p.push_back(points[node]);
              pid.push_back(node);
            }
            else
              my_p.push_back(el->centroid());
          }

          get_solution_secure(el, valcopy, my_p);

          double weight = el->volume() / elem_vol;

          // now copy nodal values to their right position
          // make mean value of cell values
          for (it = values.begin(); it != end; ++it)
          {
            const SolutionDescriptor& sd = get_solution_descriptor(it->first);
            unsigned int n_comp = sd.n_components();
            vector<double>& vals = values[it->first];
            vector<double>& newvals = valcopy[it->first];
            switch (sd.location())
            {
              case SolutionDescriptor::CELL:
                for (unsigned int i = 0; i < n_comp; i++)
                  vals[i] += newvals[i] * weight;
                break;

              case SolutionDescriptor::NODES:
                for (unsigned int i = 0; i < p.size(); i++)
                  for (unsigned int j = 0; j < n_comp; j++)
                    vals[pid[i] + j] = newvals[i + j];
                break;

              default:
                break;
            }
          }
        }
      }
    }
  }
  else
    get_solution_secure(elem, values, points);

  return flag;
}



bool
SimulationInterface::get_solution(const Elem* elem,
    ID id, vector<double>& values,
    const vector<Point>& p, bool local_coord)
{
  map<ID, vector<double> > tmp;
  tmp[id].resize(0);
  bool success = get_solution(elem, tmp, p, local_coord);
  if (success) values = tmp[id];

  return success;
}


bool
SimulationInterface::get_solution(const Atom*, map<ID, vector<double> >&)
{
  if (!is_solved()) return false;

  // TODO needs to be implemented
  return false;
}


bool
SimulationInterface::get_solution(map<ID, vector<double> >& values)
{
  if (!is_solved()) return false;

  if (values.size() == 0)
  {
    IDSet::iterator it(_plotvariable_ids.begin());
    const IDSet::iterator end(_plotvariable_ids.end());
    for ( ; it != end; ++it)
      if (get_solution_descriptor(*it).location() == SolutionDescriptor::GLOBAL)
        values[*it] = vector<double>(0);
  }

  get_solution_secure(values);
  return true;
}



void
SimulationInterface::get_solution_secure(const Elem*,
    std::map<ID, std::vector<double> >&,
    const std::vector<Point>&)
{
}

void
SimulationInterface::get_solution_secure(const Atom*,
    std::map<ID, std::vector<double> >&)
{
}

void
SimulationInterface::get_solution_secure(std::map<ID, std::vector<double> >&)
{
}



void
SimulationInterface::get_solution_secure(const Elem*,
        const std::set<ID>&, std::vector<std::map<ID, double> >&)
{
}


void
SimulationInterface::get_solution_secure(const Elem*,
        const std::vector<Point>&, const std::set<ID>&,
        std::vector<std::map<ID, double> >&)
{
}


void
SimulationInterface::build_integrated_quantities(std::vector<double>&)
{
}


void
SimulationInterface::build_integrated_quantities_description(
        std::vector<std::string>&,
        std::vector<std::string>&)
{
}



void
SimulationInterface::declare_solution_ext(const std::string& name, ID id,
    SolutionDescriptor::Type type, SolutionDescriptor::Location location,
    const std::string& units, unsigned int n_comp)
{
  _solution_descriptors[id] = SolutionDescriptor(name, id, type, location,
          units, n_comp);
  _solution_ids[name] = id;

  bool plot_all = _plotvariables.count("all");

  // check if it should be plotted
  if (_plotvariables.count(name) || plot_all)
    _plotvariable_ids.insert(id);
}


void
SimulationInterface::add_alias(const std::string& alias, ID id)
{
  _solution_ids[alias] = id;

  // check if it should be plotted
  if (_plotvariables.count(alias))
    _plotvariable_ids.insert(id);
}


ID
SimulationInterface::get_solution_id(const std::string& variable_name) const
{
  return convert_variable_name_to_id(variable_name);
}


const SolutionDescriptor&
SimulationInterface::get_solution_descriptor(ID id) const
{
  SolutionDescrMap::const_iterator it(_solution_descriptors.find(id));
  if (it != _solution_descriptors.end())
    return it->second;

  // If it does not exist, return an invalid solution descriptor
  return _invalid_descr;
}



const SolutionDescriptor&
SimulationInterface::get_solution_descriptor(const std::string& solution_name) const
{
  ID id = INVALID_ID;
  map<const string, ID>::const_iterator it(_solution_ids.find(solution_name));
  if (it != _solution_ids.end())
    id = it->second;

  return get_solution_descriptor(id);
}


ID
SimulationInterface::convert_variable_name_to_id(
    const std::string& variable_name) const
{
  return get_solution_descriptor(variable_name).id();
}



ModelOptions&
SimulationInterface::get_solver_options(void)
{
  if (!get_options().has_submodel("$Solver"))
    get_options().add_submodel("$Solver", ModelOptions());

  ModelOptions::submodel_iterator it(get_options().submodels_begin("$Solver"));

  return it->second;
}



Embracing*
SimulationInterface::create_embracing_region(
    SimulationInterface* other_simulation,
    const ModelOptions& options, bool need_mixing_coeff)
{
  Embracing* emb = NULL;
  if (other_simulation != NULL)
  {
    if (_embracings.find(other_simulation) != _embracings.end())
      emb = _embracings[other_simulation];
    else
    {
      emb = new Embracing(this, other_simulation);
      _embracings[other_simulation] = emb;
      emb->need_mixing_coeff(need_mixing_coeff);
      emb->init(options);
    }
  }

  return emb;
}



void
SimulationInterface::print_info(void)
{
  if (verbose() > 0)
  {
    int width[5] = {25, 15, 10, 12, 4};
    int tot_width = width[0] + width[1] + width[2] + width[3] + width[4];

    Messages::newline();
    Messages::info("Available solution variables:");
    ostringstream line;
    line.width(tot_width);
    line.fill('-');
    line << "";
    Messages::info(line.str());

    {
      ostringstream os;
      os << "Name";
      int w = width[0];
      os.width(w - os.tellp());
      os << "" << "Units";
      w += width[1];
      os.width(w - os.tellp());
      os << "" << "Type";
      w += width[2];
      os.width(w - os.tellp());
      os << "" << "Association";
      w += width[3];
      os.width(w - os.tellp());
      os << "" << "Plot" << Messages::endl;
      os.width(w + width[4]);
      os.fill('-');
      os << "";
      Messages::info(os.str());
    }

    SolutionDescrMap::const_iterator it(_solution_descriptors.begin());
    for ( ; it != _solution_descriptors.end(); ++it)
    {
      if (it->first == INVALID_ID) continue;

      int w = width[0];
      ostringstream os;
      os << it->second.name();
      os.width(w - os.tellp());
      os << "" << it->second.units();
      w += width[1];
      os.width(w - os.tellp());
      os << "" << it->second.type();
      w += width[2];
      os.width(w - os.tellp());
      os << "" << it->second.location();
      w += width[3];
      os.width(w - os.tellp());
      os << "" << (plot_solution(it->first) ? "y" : "n");
      Messages::info(os.str());
    }
    Messages::info(line.str());


    Messages::newline();
    ostringstream os;
    os << "Output directory      : " << get_output_directory() << endl;
    os << "Output file format    : " << get_option("output_format", "-") << endl;
    os << "Output file basename  : " << get_output_filename_prefix() << endl;
    Messages::info(os.str());
    Messages::newline();
  }

  if (verbose() > 0) do_print_info();
}
