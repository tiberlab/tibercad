// $Id$


#include <boost/filesystem/operations.hpp>

#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/solver/TiberEqSystem.h"
#include "tibercad/solver/TiberLinearSystem.h"
#include "tibercad/physics/Material.h"
#include "tibercad/atomistic/Atom.h"
#include "tibercad/geom/Boundary.h"
#include "tibercad/physics/MaterialBoundary.h"
#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/geom/EdgeObject.h"
#include "tibercad/geom/NodeObject.h"
#include "tibercad/physics/Alloy.h"
#include "tibercad/embracing/Embracing.h"
#include "tibercad/base/RuntimeException.h"

#include "tibercad/base/Variable.h"

#include "tibercad/physics/tightbinding/OpticsTB.h"

#include "solver/Sweep.h"
#include "solver/RelaxationMethod.h"
#include "tibercad/utils/Utils.h"
#include "tibercad/io/DataOutput.h"
#include "tibercad/io/Messages.h"
#include "tibercad/atomistic/AtomisticStructure.h"
#include "tibercad/atomistic/GridCells.h"
#include "tibercad/geom/MeshUtils.h"

#include "tibercad/base/libMeshDefs.h"


// LibMesh includes
#include "libmesh/system.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/mesh.h"
#include "libmesh/fe_interface.h"
#include "libmesh/periodic_boundary.h"

#include <sstream>
#include <algorithm>



using namespace std;


SimulationInterface::SimulationMap
SimulationInterface::_simulation_map;

map<string, list<boost::function<void(void)>>>
SimulationInterface::_callback_functions;

SolutionDescriptor
SimulationInterface::_invalid_descr = SolutionDescriptor();


SimulationInterface::SimulationInterface(const ModelOptions& options)
  : TiberModelObject(options),
    _environment(0),
    _is_initialized(false),
    _solve_sequence_nr(0),
    _is_task(false),
    _equilibrium_is_solved(false),
    _has_solution_vector(true),
    _symmetry(TiberCad::NONE),
    _use_cache(false),
    _verbosity(1),
    _mesh(0),
    _atomistic_structure(0)
{
  ID new_id = _simulation_map.size() + 1;
  _id = new_id;

  // register in the map of simulations
  _simulation_map[new_id] = this;

  // dummy read
  get_options().get_option("regions", "");
}


SimulationInterface::~SimulationInterface(void)
{
  map<ID, libMesh::NumericVector<double>*>::iterator it(_remembered_solutions.begin());

  map<ID, libMesh::NumericVector<double>*>::iterator end(_remembered_solutions.end());
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
  SimulationInterface* sim = nullptr;

  string type_name(type);
  string flavour = options.get_option("module_subtype", "");
  if (flavour.size() != 0)
    type_name += "_" + flavour;

  if (type_name == "sweep")
    sim = Sweep::create(options);
  else if (type_name == "selfconsistent")
    sim = RelaxationMethod::create(options);
  else if (type_name == "selfconsistent_relaxation")
    sim = RelaxationMethod::create(options);
  else if (type_name == "opticstb")
    sim = OpticsTB::create(options);

  if (sim == nullptr)
  {
    // try first without a module directory
    if ((sim = create_from_library<SimulationInterface>(type_name, options)) == 0)
    {
      sim = create_from_library<SimulationInterface>(type + "/" + type_name, options);
    }
  }

  if (sim != nullptr)
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
SimulationInterface::create_boundary(const ModelOptions& opts)
{
  if (is_task()) return;

  SimulationEnvironment* env = &(get_environment());

  if (env == nullptr)
    throw InitFailedException("Simulation which is not a task "
        "must have environment associated.");

  Device& device = env->get_device();

  // first get region names
  string id_str = opts.get_option("regions", "");
  vector<string> ids_strings;
  Utils::extract_vector(id_str, ids_strings);

  string boundary_name = opts.get_name();

  // for the numeric IDs
  vector<ID> ids;

  unsigned int n_ids = ids_strings.size();
  // if no numbers are specified we try to get them from the region name
  if (n_ids == 0)
    device.get_boundary_region_ids(boundary_name, ids);
  else
  {
    vector<ID> tmp_id;
    for (unsigned int i = 0; i < n_ids; i++)
    {
      // either it is a name or a number
      // try first name
      device.get_boundary_region_ids(ids_strings[i], tmp_id);
      if (tmp_id.size() == 0)
      {
        ostringstream s;
        s << "Physical region \'" << ids_strings[i]
          << "\' (in boundary  \'" << boundary_name
          << "\') does not exist in mesh.";
        throw InitFailedException(s.str());
      }
      ids.insert(ids.end(), tmp_id.begin(), tmp_id.end());
    }
  }

  if (ids.size() == 0)
  {
    ostringstream s;
    s << "Boundary region \'" << boundary_name <<
        "\' is not consistent with mesh.";
    throw InitFailedException(s.str());
  }

  {
    ostringstream os;
    os << "Adding boundary \'" << boundary_name << "\'";
    Messages::info(os.str());
  }


  Boundary* bnd = env->get_boundary(boundary_name);
  if (bnd == NULL)
  {
    bnd = new Boundary(boundary_name, opts);
    bnd->set_region_ids(ids);
    env->add_boundary(bnd);
  }


  for (unsigned int i = 0; i < ids.size(); i++)
  {

    bool found = false;

    MaterialBoundary* bd;
    if ((bd = device.get_boundary_object(ids[i])) != NULL)
    {
      PhysicalModel* pm = bd->get_model(get_id());
      if (pm != NULL)
      {
        ostringstream os;
        os << "Trying to add already existing boundary \'"
            << boundary_name << "\' for module "
            << get_name();
        throw InitFailedException(os.str());
      }

      ID ida = bd->get_id_A();
      ID idb = bd->get_id_B();
      if (!includes_region(ida) && !includes_region(idb))
      {
        continue;

        Messages::warning("At least part of boundary \'" + boundary_name +
            "\' does not touch any region of simulation \'" +
            get_name() + "\'");
      }

      pm = new_boundary_model(opts, bd);
      bd->add_model(pm, get_id());
      bnd->add_model(ids[i], pm);
      found = true;
    }

    EdgeObject* eo;
    if ((eo = device.get_edge_object(ids[i])) != NULL)
    {
      PhysicalModel* pm = eo->get_model(get_id());
      if (pm != NULL)
      {
        ostringstream os;
        os << "Trying to add already existing boundary \'"
            << boundary_name << "\' for module "
            << get_name();
        throw InitFailedException(os.str());
      }

      pm = new_edge_model(opts, eo);
      eo->add_model(pm, get_id());
      bnd->add_model(ids[i], pm);
      found = true;
    }

    NodeObject* no;
    if ((no = device.get_node_object(ids[i])) != NULL)
    {
      PhysicalModel* pm = no->get_model(get_id());
      if (pm != NULL)
      {
        ostringstream os;
        os << "Trying to add already existing boundary \'"
            << boundary_name << "\' for module "
            << get_name();
        throw InitFailedException(os.str());
      }

      pm = new_node_model(opts, no);
      no->add_model(pm, get_id());
      bnd->add_model(ids[i], pm);
      found = true;
    }

    if (!found)
    {
      ostringstream os;
      os << "Boundary \'" << boundary_name
                  << "\' does not exist.";
      throw InitFailedException(os.str());
    }

  }
}



void
SimulationInterface::register_callback(string& name,
    boost::function<void(void)> callback)
{
  if (!name.empty())
  {
    vector<string> tokens;

    Utils::tokenize(name, tokens, ".");

    _callback_functions[tokens[0]].push_back(callback);
  }
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
SimulationInterface::set_symmetry(TiberCad::Symmetry symmetry)
{
  _symmetry = symmetry;
}


TiberCad::Symmetry
SimulationInterface::get_symmetry(void) const
{
  return(_symmetry);
}



void
SimulationInterface::do_print_info(void)
{
}


void
SimulationInterface::restrict_solve_to_subdomains(const set<ID>& ids,
    const vector<string>& variables)
{
  // find_excluded_dofs wants the set of ids where variables have to
  // be excluded
  set<ID> all_ids, excluded_ids;
  get_region_ids(all_ids);
  set_difference(all_ids.begin(), all_ids.end(), ids.begin(), ids.end(),
      std::inserter(excluded_ids, excluded_ids.end()));

  exclude_solve_from_subdomains(excluded_ids, variables);
}



void
SimulationInterface::exclude_solve_from_subdomains(const set<ID>& ids,
    const vector<string>& variables)
{
  // if we have no environment then we have probably no mesh
  // and we can go out immediately
  if (!has_environment()) return;

  // we have to activate our elements
  get_environment().prepare_for_solve();

  find_excluded_dofs(ids, variables);

  for (set<ID>::iterator it(ids.begin()); it != ids.end(); ++it)
    for (int i = 0; i < variables.size(); ++i)
      _excluded_domains[*it].insert(variables[i]);
}



void
SimulationInterface::find_excluded_dofs(const std::set<ID>& ids,
    const vector<string>& variables)
{
  // if we have no environment then we have probably no mesh
  // and we can go out immediately
  if (!has_environment())
    return;

  // a set of DoFs for each system
  vector<IDHashSet> dofsets(1);
  vector<IDHashSet> bd_dofs(1);

  TiberEqSystem& tiber_sys = get_equation_system<TiberEqSystem>();
  libMesh::System* system = tiber_sys.get_libmesh_system();

  // In the remote case that system is nullptr we return immediately
  if (system == nullptr)
    return;


  const libMesh::DofMap& dof_map = system->get_dof_map();
  vector<unsigned int> dof_indices;

  // contains the 'active' DoFs, false means inactive
  vector<bool> used_dofs(dof_map.n_dofs(), false);
  // contains the 'inactive' DoFs, true means inactive
  vector<bool> unused_dofs(dof_map.n_dofs(), false);

  // for each variable, tell if it is used in the
  // excluded domains:
  // true means include, false means exclude
  vector<bool> var(dof_map.n_variables(), true);


  for (int i = 0; i < variables.size(); ++i)
    var[system->variable_number(variables[i])] = false;

  if (variables.empty())
  {
    // exclude all variables
    var = vector<bool>(var.size(), false);
  }


  const MeshBase& mesh = get_mesh();

  MeshBase::const_element_iterator el =
                                  this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  this->active_local_elements_end();

  // loop over all active elements
  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    ID subdomain = elem->subdomain_id();

    // if the subdomain ID is not found in the given set
    // it means all its DoFs need to be included
    int restricted = ids.count(subdomain);

    for (int i = 0; i < var.size(); ++i)
    {
      dof_map.dof_indices(elem, dof_indices, i);
      if (!restricted || var[i])
      {
        for (int j = 0; j < dof_indices.size(); ++j)
          used_dofs[dof_indices[j]] = true;
      }
      else
      {
        for (int j = 0; j < dof_indices.size(); ++j)
          unused_dofs[dof_indices[j]] = true;
      }
    }
  }

  // now we have all used DoFs in used_dofs marked with true
  for (size_t i = 0; i < used_dofs.size(); ++i)
    if (!used_dofs[i])
      dofsets[0].insert(i);
    else
    {
      // it may be on the boundary, in that case we have to treat it
      // in a special way
      if (unused_dofs[i])
        bd_dofs[0].insert(i);
    }

  tiber_sys.set_excluded_dofs(dofsets[0], bd_dofs[0], ids);
}


bool
SimulationInterface::includes_region(ID region_id) const
{
  return get_environment().contains_region(region_id);
}


bool
SimulationInterface::includes_regions(std::set<ID> region_ids) const
{
  return includes(get_environment().get_region_ids().begin(),
      get_environment().get_region_ids().end(),
      region_ids.begin(), region_ids.end());
}


PhysicalModel*
SimulationInterface::new_bulk_model(const ModelOptions& options,
    const Material* material)
{
  PhysicalModel* pm = create_bulk_model(options, material);

  if (pm != nullptr)
    _physical_models.insert(pm);

  return pm;
}


//BoundaryProperties*
//SimulationInterface::new_boundary_model(const ModelOptions& options)
//{
//  BoundaryProperties* bp = create_boundary_model(options);

//  return bp;
//}


PhysicalModel*
SimulationInterface::new_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary)
{
  PhysicalModel* pm = create_boundary_model(options, boundary);

  if (pm != nullptr)
    _boundary_models.insert(pm);

  return pm;
}


PhysicalModel*
SimulationInterface::new_edge_model(const ModelOptions& options,
    const EdgeObject* edge)
{
  PhysicalModel* pm = create_edge_model(options, edge);

  //if (pm != nullptr)
  //  _edge_models.insert(pm);

  return pm;
}


PhysicalModel*
SimulationInterface::new_node_model(const ModelOptions& options,
    const NodeObject* node)
{
  PhysicalModel* pm = create_node_model(options, node);

  //if (pm != nullptr)
  //  _node_models.insert(pm);

  return pm;
}





Material*
SimulationInterface::get_material(const Elem* elem) const
{
  return get_environment().get_device().get_material(elem);
}



PhysicalModel*
SimulationInterface::_get_bulk_model(const Elem* elem) const
{
  PhysicalModel* mod = nullptr;
  const Material* mat = get_environment().get_device().get_material(elem);
  if (mat != nullptr)
    mod = mat->get_model(get_id());

  return mod;
}


/*
PhysicalModel*
SimulationInterface::_get_bulk_model(const Atom& atom, bool parent) const
{
  PhysicalModel* mod = nullptr;
  const Material* mat = get_atomistic_structure()->get_material(atom, parent);
  if (mat != nullptr)
    mod = mat->get_model(get_id());
 
  return mod;
}
*/

PhysicalModel*
SimulationInterface::_get_bulk_model(const Atom& atom1, const Atom& atom2,
    bool parent) const
{
  PhysicalModel* mod = nullptr;
  const Material* mat = get_atomistic_structure()->get_material(atom1, atom2, parent);

  if (mat != nullptr) mod = mat->get_model(get_id());

  return mod;
}


PhysicalModel*
SimulationInterface::_get_interface_model(const Elem* elem, int side) const
{
  PhysicalModel* mod = nullptr;
  MaterialBoundary* mb =
      get_environment().get_device().get_boundary_object(elem, side);
  if (mb != nullptr)
    mod = mb->get_model(get_id());

  if (mod != nullptr)
    mod->set_material(get_material(elem));

  return mod;
}



PhysicalModel*
SimulationInterface::_get_edge_model(const Elem* elem, int edge) const
{
  PhysicalModel* mod = nullptr;
  EdgeObject* eo =
      get_environment().get_device().get_edge_object(elem, edge);
  if (eo != nullptr)
    mod = eo->get_model(get_id());

  if (mod != nullptr)
    mod->set_material(get_material(elem));

  return mod;
}



PhysicalModel*
SimulationInterface::_get_node_model(const Elem* elem, int node) const
{
  PhysicalModel* mod = nullptr;
  NodeObject* no =
      get_environment().get_device().get_node_object(elem, node);
  if (no != nullptr)
    mod = no->get_model(get_id());

  if (mod != nullptr)
    mod->set_material(get_material(elem));

  return mod;
}


libMesh::MeshBase::const_element_iterator
SimulationInterface::active_local_elements_begin(void) const
{
  return(this->get_mesh().active_local_subdomain_set_elements_begin(
      this->get_environment().get_region_ids()));
}

libMesh::MeshBase::const_element_iterator
SimulationInterface::active_local_elements_end(void) const
{
  return(this->get_mesh().active_local_subdomain_set_elements_end(
      this->get_environment().get_region_ids()));
}



libMesh::MeshBase::element_iterator
SimulationInterface::active_local_elements_begin(void)
{
  return(this->get_mesh().active_local_subdomain_set_elements_begin(
      this->get_environment().get_region_ids()));
}

libMesh::MeshBase::element_iterator
SimulationInterface::active_local_elements_end(void)
{
  return(this->get_mesh().active_local_subdomain_set_elements_end(
      this->get_environment().get_region_ids()));
}

double
SimulationInterface::get_mesh_units(void) const
{
  return get_environment().get_device().get_mesh_units();
}



void
SimulationInterface::prepare(void)
{

  Messages m;

  //m.newline();
  //m.info("Setting up simulation of type \'"
  //    + get_name() + "\' ...");
  m.indent();


  // the following stuff does not need to be done if the simulation is a task
  if (!is_task())
  {

    SimulationEnvironment* env = &get_environment();
    Device& device = env->get_device();

    IDSet phys_regions;
    string physreg = get_option("regions", "all");
    device.extract_physical_regions(physreg, phys_regions);


    //
    // and now... the boundary conditions
    //

    // we accept the following keywords for boundaries:
    //   Contact, Boundary, Interface
    vector<string> keys(4);
    keys[0] = "Contact";
    keys[1] = "Boundary";
    keys[2] = "Interface";
    keys[3] = "BoundaryCondition";

    for (size_t i = 0; i < keys.size(); ++i)
    {
      ModelOptions::const_submodel_iterator it =
          get_options().submodels_begin(keys[i]);
      const ModelOptions::const_submodel_iterator end =
          get_options().submodels_end(keys[i]);

      for ( ; it != end; ++it)
        create_boundary(it->second);

      // we remove them, so in the following we have only models
      //physopts.delete_submodels(keys[i]);
    }




    // the bulk physical models, all together
    ModelOptions bulk_opts;


    //
    // Next, we create all lower dimensional submodels
    // At the same time, we put all bulk Physics blocks together
    //
    // NOTE: only definitions with correct space dimensions will be
    //       added to the different PhysicalObject instances

    ModelOptions::const_submodel_iterator ph_it =
        get_options().submodels_begin("Physics");
    const ModelOptions::const_submodel_iterator ph_end =
        get_options().submodels_end("Physics");

    for ( ; ph_it != ph_end; ++ph_it)
    {
      ModelOptions physopts = ph_it->second;

      ModelOptions::submodel_iterator it = physopts.submodels_begin();
      const ModelOptions::submodel_iterator end = physopts.submodels_end();
      while (it != end)
      {
        ModelOptions::submodel_iterator tmpit(it);
        ++it;

        const ModelOptions& bdopts = tmpit->second;

        string physreg = bdopts.get_option("regions", "all");

        // if "all", it cannot be a lower dimensional region
        if (physreg == "all")
          continue;

        bool add = false;

        // It might be some lower dim model

        vector<string> ids_strings;
        Utils::extract_vector(physreg, ids_strings);

        for (unsigned int i = 0; i < ids_strings.size(); i++)
        {
          vector<ID> region_ids;
          device.get_boundary_region_ids(ids_strings[i], region_ids);

          // if it is no boundary region, we continue to the next region name
          if (region_ids.size() == 0)
          {
            // it could be the name of the Boundary object
            Boundary* bnd = env->get_boundary(ids_strings[i]);
            if (bnd != NULL)
              bnd->get_region_ids(region_ids);
          }

          if (region_ids.size() == 0)
            continue;

          // now it must be a lower dim model
          add = true;

          // NOTE: the model will only be added if a corresponding boundary ID
          // has been found
          for (unsigned int reg = 0; reg < region_ids.size(); reg++)
          {
            ID id = region_ids[reg];

            MaterialBoundary* bd;
            if ((bd = device.get_boundary_object(id)) != NULL)
            {
              PhysicalModel* pm = bd->get_model(get_id());
              if (pm == NULL)
              {
                // create the default model on the fly
                pm = new_boundary_model(ModelOptions(), bd);
                bd->add_model(pm, get_id());
              }

              // now it's there
              pm->get_options().add_submodel(tmpit->first, bdopts);
            }

            EdgeObject* eo;
            if ((eo = device.get_edge_object(id)) != NULL)
            {
              PhysicalModel* pm = eo->get_model(get_id());
              if (pm == NULL)
              {
                // create the default model on the fly
                pm = new_edge_model(ModelOptions(), eo);
                eo->add_model(pm, get_id());
              }

              // now it's there
              pm->get_options().add_submodel(tmpit->first, bdopts);
            }

            NodeObject* no;
            if ((no = device.get_node_object(id)) != NULL)
            {
              PhysicalModel* pm = no->get_model(get_id());
              if (pm == NULL)
              {
                // create the default model on the fly
                pm = new_node_model(ModelOptions(), no);
                no->add_model(pm, get_id());
              }

              // now it's there
              pm->get_options().add_submodel(tmpit->first, bdopts);
            }
          }
        }

        // remove it from the map
        // NOTE: this means, we cannot specify the same model for bulk
        //       boundaries in a single block
        // TODO: maybe this can be relaxed?
        if (add)
          physopts.delete_submodel(tmpit);
      }

      // now add the remaining content to the common Physics options
      bulk_opts += physopts;
      bulk_opts.set_key((ph_it->second).get_key());
    }
    m.unindent();


    //
    // now we have to create the bulk models
    //
    //m.newline();
    m.info("Creating bulk models... ");
    m.indent();

    // we have to do this for each material!
    IDSet::iterator reg_it(phys_regions.begin());
    const IDSet::iterator reg_end(phys_regions.end());
    for ( ; reg_it != reg_end; ++reg_it)
    {
      ID reg_id = *reg_it;
      Material* mat = device.get_material(reg_id);

      m.info("Region: "+device.get_region_name(reg_id));
      if (mat == NULL)
      {
        ostringstream s;
        s << "Physical region " << device.get_region_name(reg_id) <<
            " has no material associated!";
        Messages::warning(s.str());
        continue;
        //throw InitFailedException(s.str());
      }

      // we only continue if the model has not already been added
      // this is important as a material can be assigned to different
      // regions
      // TODO (Is this true ??)
      if (mat->get_model(get_id()) == NULL)
      {

        // the crystal structure
        //string crystal_structure(mat->get_structure());

        // we make a copy so we can safely delete submodels
        ModelOptions opts(bulk_opts);

        // we add the crystal structure for bulk materials as this could
        // lead to different model implementations
        //opts["crystal_structure"] = crystal_structure;

        // recursively eliminate all models not defined on the current region
        _eliminate_unneeded_submodels(opts, reg_id);


        // here we actually create the model
        PhysicalModel* model = new_bulk_model(opts, mat);

        // NOTE: model could be NULL, but we don't care about. Who tells us that
        // every simulation necessarily needs a model?
        mat->add_model(model, get_id());
      }
    }
    m.unindent();

    // prepare some of the environments internals (lists of elements etc.)
    _environment->prepare();
  }


  // setup the solution variables
  setup_solution_variables();
}


void
SimulationInterface::_eliminate_unneeded_submodels(ModelOptions& opts, ID reg_id) const
{
  const Device& device = get_environment().get_device();

  auto it = opts.submodels_begin();
  const auto end = opts.submodels_end();
  for ( ; it != end; )
  {
    auto current = it++;
    // shall we delete the block?
    bool del = false;

    const ModelOptions& modopts = current->second;

    // we have to check if it should be built for the current region
    // TODO to not allow for errors

    IDSet regs;
    string physreg = modopts.get_option("regions", "all");
    device.extract_physical_regions(physreg, regs);

    if (regs.count(reg_id) == 0) del = true;


    if (del)
    {
      // we add the crystal structure for bulk materials as this could
      // lead to different model implementations
      //modopts.set_option("crystal_structure", crystal_structure);

      // we set the name to the model type if not explicitly
      // given by user
      //if (!(mapit->second).find_option("name"))
      //  (mapit->second)["name"] = mapit->first;
      opts.delete_submodel(current);
    }
    else
      _eliminate_unneeded_submodels(current->second, reg_id);
  }
}





void
SimulationInterface::setup_environment(Device& device, const set<ID>& region_numbers)
{
  if (!is_task())
  {
    if (_environment != nullptr) delete _environment;
    _environment = new SimulationEnvironment(device, region_numbers);

    this->set_communicator(device.get_communicator());

    // get atomistic structure
    setup_atomistic_structure();

    // get the mesh pointer
    setup_mesh();
    if (_mesh == nullptr)
      throw InitFailedException("No simulation mesh provided for \'" + get_name() + "\'");

  }
}


void
SimulationInterface::set_mesh(MeshBase* mesh)
{
  get_environment().set_mesh(mesh);
  _mesh = mesh;
}



void
SimulationInterface::setup_mesh(void)
{
  if (get_option("atomistic_mesh", false))
  {
    libMesh::UnstructuredMesh* mesh = new Mesh(get_solver_communicator(), 3);
    if (get_atomistic_structure() == nullptr)
      throw InitFailedException(get_name() + ": could not find atomistic structure");

    get_atomistic_structure()->create_conformal_grid(*mesh);
    get_environment().set_mesh(mesh);
  }

  // TODO allow reading individual mesh
  if (has_option("meshfile"))
  {
    //NOTE: we need to also fix from where to get mesh_units!
    throw InitFailedException("Reading mesh in modules not yet implemented.");
  }

  _mesh = &get_environment().get_mesh();
}


void
SimulationInterface::setup_atomistic_structure(void)
{
  string name(get_option("atomistic_structure", ""));
  if (!name.empty())
  {
    _atomistic_structure = get_environment().get_device().get_atomistic_structure(name);
    if (_atomistic_structure == nullptr)
      throw InitFailedException("No atomistic structure \'" + name + "\' found "
          "for simulation \'" + get_name());

    if (!includes_regions(_atomistic_structure->get_IDset()))  
      Messages::error("Module will restrict the atomistic structure");      
    
  }
}


void
SimulationInterface::reinit(void)
{
  // reinitialize all models
  set<PhysicalModel*>::iterator it(get_physical_models().begin());
  set<PhysicalModel*>::iterator end(get_physical_models().end());
  for ( ; it != end; ++it)
    (*it)->reinit();

  it = get_interface_models().begin();
  end = get_interface_models().end();
  for ( ; it != end; ++it)
    (*it)->reinit();

  do_reinit();
}


void
SimulationInterface::setup_mpi_comm(void)
{
  if (this->has_environment())
  {

    // this is just a guess
    if (_atomistic_structure != nullptr)
      this->set_solver_communicator(this->get_communicator());
    else
      this->set_solver_communicator(this->get_mesh().comm());
  }
  else
    this->set_communicator(TiberCad::get_mpi_comm());

}



void
SimulationInterface::init(void)
{
  if (verbose() > 0)
  {
    Messages m;
    m.newline();
    m.frameline(">>>>",'>');

    ostringstream os;
    os << "Simulation options for " << get_name() << " (" <<
        get_type() << ")";

    m.info(os.str());
    m.indent();
  }

  if (this->has_environment())
  {
    if (!get_option("ignore_symmetry", false))
      set_symmetry(get_environment().get_device().get_symmetry());
    else
    {
      Messages::warning("'ignore_symmetry' is specified, so any spatial"
          " symmetry is ignored in this module.");
    }
  }


  bool plot_all = _plotvariables.count("all");

  // check which variables should be plotted
  for (auto&& var : _solution_descriptors)
  {
    if (_plotvariables.count(var.second.name()) || plot_all)
      _plotvariable_ids.insert(var.first);

    // a string of type "-pippo" prevents plotting of "pippo"
    if (_plotvariables.count("-" + var.second.name()))
      _plotvariable_ids.erase(var.first);
  }

  if (!_is_initialized)
  {
    Messages::debug("Initialize " + get_name() + "... ");

    // build name for equation systems
    create_equation_system_name();

    if (_environment != nullptr)
    {
      _environment->prepare_for_solve();
      _scaling.set_calc_mesh_units(get_mesh_units());
    }

    this->setup_mpi_comm();

    if (verbose() > 0)
    {
      ostringstream os;
      os << "MPI: rank " << get_communicator().rank() <<
          " of communicator with size " <<
          get_communicator().size() << "\n";
      os << "MPI: solver parallelized on " <<
          get_solver_communicator().size() << " processes";
      Messages::info(os.str());
    }

    using namespace boost::filesystem;

    // check output and scratch path
    path outpath(get_output_directory());
    if (!exists(outpath))
    {
      // we catch any error here without doing anything yet
      try {
        create_directories(outpath);
      }
      catch (...) {}
    }

    if (!(exists(outpath) && is_directory(outpath)))
    {
      string msg("Cannot create or use '");
      msg += outpath.string() + "' as output directory.";
      throw InitFailedException(msg);
    }

    if (get_output_directory() != get_scratch_directory())
    {
      path outpath(get_scratch_directory());
      if (!exists(outpath))
      {
        // we catch any error here without doing anything yet
        try {
          create_directories(outpath);
        }
        catch (...) {}
      }

      if (!(exists(outpath) && is_directory(outpath)))
      {
        string msg("Cannot create or use '");
        msg += outpath.string() + "' as output directory.";
        throw InitFailedException(msg);
      }
    }

    _use_cache = get_option("use_data_cache", _use_cache);

    _verbosity = get_option("verbose", _verbosity);
    _verbosity = get_option("verbosity", _verbosity);
    do_init();

  }

  _is_initialized = true;

  Messages::debug("init of " + get_name() + " done");


  if (verbose() > 0)
  {
    print_info();

    Messages m;
    m.newline();

    set<string> names;
    if (has_environment())
    {
      get_environment().get_region_names(names);
      for (set<string>::const_iterator it(names.begin()); it != names.end(); ++it)
      {
        Device& dev = get_environment().get_device();
        set<ID> ids;
        dev.get_active_region_ids(*it, ids);
        assert(ids.size() != 0);
        const Material* mat = dev.get_material(*ids.begin());
        PhysicalModel* mod =
            mat->get_model(get_id());
        if (mod != nullptr)
        {
          ostringstream os;
          os << "# Region " << *it << ", " << mat->get_name();
          if (mat->is_alloy())
            os << " (x = " <<
            static_cast<const Alloy*>(mat)->get_molar_fraction() << ")";
          m.info(os.str());
          m.indent();
          mod->print_info();
          m.unindent();
          m.newline();
        }
      }
    }

    // some dummy reads to shut up the automatic check
    get_option("save_state", "");
    if (has_option("load_state", ""))
      get_option("solve_after_load", "");

    get_options().check_unused(1);

    m.unindent();
    m.newline();
    m.frameline("<<<<",'<');

  }
}


unsigned int
SimulationInterface::get_equation_system_id(ID i) const
{
  const libMesh::System& sys = get_equation_system<libMesh::System>(i);
  return(sys.number());
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
  SimulationInterface* sim = nullptr;

  SimulationMap::iterator it(_simulation_map.begin());
  SimulationMap::iterator end(_simulation_map.end());

  if (it != end)
  {
    // this is actually a bad idea
    //if (name == "") // we just take the first we can find ...
    //  sim = it->second;
    //else
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




SolutionProvider
SimulationInterface::find_solution_provider(const string& simulation,
    const string& solution)
{
  SolutionProvider result(nullptr, INVALID_ID);

  if (!simulation.empty())
  {
    vector<string> tokens;

    Utils::tokenize(simulation, tokens, ".");

    result.first = find_simulation(tokens[0]);
    if (result.first != nullptr)
    {
      if (tokens.size() > 1)
        result.second = result.first->get_solution_id(tokens[1]);
      else
        result.second = result.first->get_solution_id(solution);
    }
  }

  return result;
}

void
SimulationInterface::check_nonlinear_step(libMesh::NumericVector<Number>& dx)
{
  do_check_nonlinear_step(dx);
}


void
SimulationInterface::do_check_nonlinear_step(libMesh::NumericVector<Number>& dx)
{
  static_cast<void>(dx);
}

void
SimulationInterface::get_region_ids(std::set<ID>& region_ids) const
{
  region_ids.clear();
  if (has_environment())
    region_ids = get_environment().get_region_ids();
}


const std::set<ID>&
SimulationInterface::get_region_ids(void) const
{
  // this is a bit of a quirk, but as long as region ids are in
  // the environment, I don't know a better way without having to
  // explicitly call all the time get_environment()
  static const std::set<ID> empty_set;
  if (has_environment())
    return(get_environment().get_region_ids());

  return(empty_set);
}




libMesh::EquationSystems&
SimulationInterface::get_equation_systems(void) const
{
  return _environment->get_device().get_equation_systems(_mesh);
}


void
SimulationInterface::clear_systems(void)
{
  for (size_t i = 0; i < _systems.size(); ++i)
    delete _systems[i];

  _systems.clear();
}


ID
SimulationInterface::create_equation_system(const std::string& type,
    const std::string& block,
    const ModelOptions& options)
{
  ID newid = _systems.size();

  ostringstream os;
  os << get_equation_system_name();
  
  //if (newid > 0)
  os << "_" << newid;

  ModelOptions opts(options);
  opts += get_solver_options(block);

  opts.set_option("simulation", get_name());
  opts.find_option("simulation");

  TiberEqSystem* sys = TiberEqSystem::create(get_equation_systems(),
      os.str(), type, opts);

  ModelOptions::submodel_iterator pbit(
      get_options().submodels_begin("PeriodicBoundary"));

  for ( ; pbit != get_options().submodels_end("PeriodicBoundary"); ++pbit)
  {
    ModelOptions& pb_opts = pbit->second;

    libMesh::RealVectorValue period;
    pb_opts.get_option("periodicity", period);

    vector<string> bdpair(2, "");
    pb_opts.get_option("boundary_pair", bdpair);

    vector<ID> bdids;
    const Device& dev = get_environment().get_device();
    dev.get_boundary_region_ids(bdpair[0], bdids);
    ID tmp = bdids[0];
    dev.get_boundary_region_ids(bdpair[1], bdids);
    bdids.resize(2);
    bdids[1] = bdids[0];
    bdids[0] = tmp;

    libMesh::PeriodicBoundary pb(period);
    pb.myboundary = bdids[0];
    pb.pairedboundary = bdids[1];
    libMesh::DofMap& dof_map = sys->get_libmesh_system()->get_dof_map();
    dof_map.add_periodic_boundary(pb);
    ostringstream os;
    os << "set up periodic boundary: " << endl
        << "    periodicity: (" << period(0) << ", " << period(1) << ", "
        << period(2) << ")" << endl
        << "    pair       : " << bdpair[0] << " (" << bdids[0] << "), "
                               << bdpair[1] << " (" << bdids[1] << ")";
    Messages::info(os.str());
  }



  _systems.push_back(sys);

  return newid;
}


void
SimulationInterface::solve_equilibrium(void)
{

  if (!_equilibrium_is_solved)
  {

    assert(is_initialized());

    if (_environment != nullptr)
      _environment->prepare_for_solve();


    do_equilibrium();

    increment_solve_sequence_number();


    _equilibrium_is_solved = true;

  }
}




void
SimulationInterface::increment_solve_sequence_number(void)
{
  ++_solve_sequence_nr;

  // here we also call all callbacks
  map<string, list<boost::function<void(void)>>>::iterator mit =
      _callback_functions.find(get_name());

  if (mit != _callback_functions.end())
  {
    list<boost::function<void(void)>>::iterator it((mit->second).begin());
    for ( ; it != (mit->second).end(); ++it)
    {
      (*it)();
    }
  }
}




void
SimulationInterface::add_plot_variable(const std::string& name)
{
  _plotvariables.insert(name);
}



void
SimulationInterface::add_plot_variable(ID id)
{
  _plotvariable_ids.insert(id);
}

void
SimulationInterface::remove_plot_variable(ID id)
{
  _plotvariable_ids.erase(id);
}

void
SimulationInterface::solve(void)
{

  Messages m;
  m.newline();
  m.frameline(">>>>",'-',get_name());

  m.indent();

  if (_environment != nullptr)
    _environment->prepare_for_solve();

  // call reinitialization
  reinit();

  assert(is_initialized());

  Utils::Timer tt;


  try
  {
    if (!load_state() || get_option("solve_after_load", false))
    {
      do_solve();
      save_state();
    }
  }
  catch (SolveFailedException& e)
  {
    ostringstream os;
    os << "Solve time: " << tt.elapsed_string();
    Messages::newline();
    Messages::info(os.str());

    ostringstream s;
    s << get_name() << ": " << e.what();
    throw SolveFailedException(s.str());
  }
  catch (runtime_error& e)
  {
    ostringstream os;
    os << "Solve time: " << tt.elapsed_string();
    Messages::newline();
    Messages::info(os.str());

    ostringstream s;
    s << get_name() << ": " << e.what();
    throw RuntimeException(s.str());
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

  string unit("kB");
  size_t mem = _data_cache.get_memory_size() / 1024;
  if (mem >= 1024)
  {
    mem /= 1024;
    unit = "MB";
    if (mem >= 1024)
    {
      mem /= 1024;
      unit = "GB";
    }
    ostringstream os;
    os << "Current result cache size: " << mem << " " << unit;
    Messages::info(os.str());
  }

  // plot at this point
  //plot();

  // invalidate cache
  _data_cache.flush();

  increment_solve_sequence_number();

  m.unindent();

  ostringstream os;
  os << "Solve time: " << tt.elapsed_string();
  Messages::newline();
  Messages::info(os.str());

  this->analyze_errors();

  Messages::frameline("<<<<",'-');

}


void
SimulationInterface::analyze_errors(void)
{
  if (get_options().has_submodel("ErrorAnalysis"))
  {
    Messages m;
    m.info("Performing analysis of numerical errors");
    m.indent();
    this->do_analyze_errors((get_options().submodels_begin("ErrorAnalysis"))->second);
    m.unindent();
  }
}

void
SimulationInterface::do_analyze_errors(const ModelOptions& options)
{
  static_cast<void>(options);
  Messages::info("Not implemented yet");
}

libMesh::NumericVector<double>&
SimulationInterface::do_get_solution_vector(void)
{
  assert(_systems.size() > 0);

  get_equation_system<TiberEqSystem>(0).get_local_solution_vector().close();
  return get_equation_system<TiberEqSystem>(0).get_local_solution_vector();
}



void
SimulationInterface::do_set_solution_vector(
    const libMesh::NumericVector<double>& new_solution)
{
  get_solution_vector() = new_solution;
  get_solution_vector().close();
  if (_systems.size() > 0)
    get_equation_system<libMesh::System>(0).update();

}






PhysicalModel*
SimulationInterface::create_bulk_model(const ModelOptions&,
    const Material*) const
{
  return nullptr;
}



PhysicalModel*
SimulationInterface::create_boundary_model(const ModelOptions&,
    const MaterialBoundary*) const
{
  return nullptr;
}



PhysicalModel*
SimulationInterface::create_edge_model(const ModelOptions&,
    const EdgeObject*) const
{
  return nullptr;
}



PhysicalModel*
SimulationInterface::create_node_model(const ModelOptions&,
const NodeObject*) const
{
  return nullptr;
}


/*
void
SimulationInterface::get_integrated_quantities(std::vector<double>& values)
{

  if (_environment != nullptr)
    get_environment().prepare_for_solve();

  values.resize(0);
  build_integrated_quantities(values);
}
*/


void
SimulationInterface::plot(void)
{
  if (_environment != nullptr)
    get_environment().prepare_for_solve();

  do_plot();
  project_on_tensor_grid();


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
SimulationInterface::get_scratch_directory(void) const
{
  return get_option("scratchpath", ".");
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
  //if (this->get_communicator().rank() != 0)
  //  return;

  plot_meshdata();
  plot_atomisticdata();

  if (get_communicator().rank() == 0)
    plot_globaldata();
}




void
SimulationInterface::plot_meshdata(void)
{

  const MeshBase& mesh = get_mesh();

  // The device communicator might be larger than the mesh communicator.
  // In that case, we do the work only on one group of processes

  bool do_write = true;

  if (this->get_communicator().size() > this->get_mesh().comm().size())
  {
    unsigned int dev_rank = this->get_communicator().rank();
    this->get_mesh().comm().broadcast(dev_rank);
    unsigned int min_rank = dev_rank;
    this->get_communicator().min(min_rank);

    if (dev_rank != min_rank)
      return;
      //do_write = false;
  }

  // we write only on mesh part associated to this simulation

  MeshBase::const_element_iterator it = this->active_local_elements_begin();
  MeshBase::const_element_iterator end = this->active_local_elements_end();
  //libMesh::MeshBase::const_element_iterator it = mesh.active_elements_begin();
  //const libMesh::MeshBase::const_element_iterator end = mesh.active_elements_end();

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
        if (node_conn[subdomain].count(elem->node_id(n)) == 0)
        {
          unsigned int nodeid = node_conn[subdomain].size();
          node_conn[subdomain][elem->node_id(n)] = make_pair(nodeid, 1);
        }
        else
          (node_conn[subdomain][elem->node_id(n)].second)++;
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


  end = this->active_local_elements_end();
  for (it = this->active_local_elements_begin(); it != end; ++it)
  {
    const Elem* elem = *it;

    ID subdomain = elem->subdomain_id();

    // if subdomain is not in data, we go to the next element
    if (data.find(subdomain) == data.end())
      continue;

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
          // Attention: we get the same data on every node, but we need to put only
          // a single data set (for one node)
          vec.insert(vec.end(), sol.data(), sol.data() + descr.n_components());
          break;

        case SolutionDescriptor::NODES:
          for (unsigned int n = 0; n < elem->n_nodes(); n++)
          {
            unsigned int ncomp = descr.n_components();
            unsigned int index =
              ncomp * node_conn[subdomain][elem->node_id(n)].first;
            unsigned short w = node_conn[subdomain][elem->node_id(n)].second;
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
    unique_ptr<DataOutput> writer(DataOutput::create(formats[i]));
    if ((writer.get() != nullptr) && (data.size() > 0))
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
    const SolutionDescriptor& sd = get_solution_descriptor(id);
    if (sd.location() == SolutionDescriptor::GLOBAL)
      values.insert(make_pair(id, vector<double>(sd.n_components())));
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
      file << "# Global data for simulation: " << get_name() << endl;
      //file << "# Data:" << endl;
      file << "#\n";

      map<ID, vector<double> >::iterator it(values.begin());


      for (it = values.begin(); it != values.end(); ++it)
      {
        const SolutionDescriptor& descr = get_solution_descriptor(it->first);
        file << "\n# " << descr.name() << "  ";
        if (descr.units() != "")
          file << "(" << descr.units() << ")  ";
        if (descr.type() == SolutionDescriptor::NTUPLE)
          file << descr.n_components() << "-tuple";
        else
          file << descr.type();
        file << "\n";

        const vector<double>& vals = it->second;
        for (unsigned int i = 0; i < vals.size(); i++)
          file << vals[i] << "  ";

        file << "\n";
      }

      file.close();
    }
  }
  /*else // temporary only !
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
  }*/
}



void
SimulationInterface::save_state(const string& file)
{
  if (get_option("save_state", false) || !file.empty())
  {
    string f(file);
    if (f.empty())
      f = get_output_directory() + "/" +
        get_output_filename() + ".tsv";

    ofstream of(f.c_str(), ios_base::binary);
    if (!of.good()) throw InitFailedException("Cannot use " + f + " for writing.");

    Messages::newline();
    Messages::info("Writing state to " + f);
    do_save_data(of);
  }
}



bool
SimulationInterface::load_state(const string& file)
{
  bool loaded = false;

  string f(file);
  if (f.empty())
    f = get_option("load_state", "");

  if (!f.empty())
  {

    ifstream in(f.c_str(), ios_base::binary);
    if (!in.good()) throw InitFailedException("Cannot use " + f + " for reading.");

    Messages::info("Reading state from " + f);
    do_load_data(in);
    loaded = true;

    get_options().set_option("load_state", "");

  }

  return loaded;
}



void
SimulationInterface::do_save_data(ostream& os)
{

  // NOTE we always use windows line endings to
  //      have better portability of the files
  string eol("\r\n");


  // first write all variables
  os << "<variables>" << eol;
  VariableValue::iterator vit(VariableValue::begin());
  const VariableValue::iterator vend(VariableValue::end());
  for ( ; vit != vend; ++vit)
  {
    os << (*vit)->get_name() << " " << (*vit)->get_value_string() << eol;
  }
  os << "</variables>" << eol;


  for (size_t i = 0; i < _systems.size(); ++i)
  {
    // then the data
    os << "<data>" << eol;

    const libMesh::NumericVector<Number>& solution =
        get_equation_system<TiberEqSystem>(i).get_solution_vector();

    for (size_t i = 0; i < solution.size(); ++i)
    {
      double val = solution(i);
      os.write(reinterpret_cast<char*>(&val), sizeof(double));
    }

    os << eol << "</data>" << eol << flush;
  }
}




void
SimulationInterface::do_load_data(istream& is)
{

  const streamsize bufsize = 256;
  char buf[bufsize];

  map<string, double> values;

  string keyword("<variables>");

  // NOTE we compare with an explicit number of characters
  //      to not get confused if there is a \r

  is.getline(buf, bufsize);
  while (is.good() && (keyword.compare(0, keyword.size(), buf, keyword.size()) != 0))
  {
    is.getline(buf, bufsize);
  }

  if (!is.good()) throw InitFailedException("Bad datafile (missing variables block)");

  keyword = "</variables>";
  is.getline(buf, bufsize);
  while (is.good() && (keyword.compare(0, keyword.size(), buf, keyword.size()) != 0))
  {
    istringstream ss(buf);
    string name;
    double value;
    ss >> name >> value;
    values[name] = value;
    is.getline(buf, bufsize);
  }

  map<string, double>::iterator vit(values.begin());
  const map<string, double>::iterator vend(values.end());
  for ( ; vit != vend; ++vit)
  {
    VariableValue::set_variable_value(vit->first, vit->second);
  }

  if (!is.good()) throw InitFailedException("Bad datafile (missing data block?)");
  values.clear();

  bool has_read = false;

  for (size_t i = 0; i < _systems.size(); ++i)
  {
    keyword = "<data>";
    is.getline(buf, bufsize);
    while (is.good() && (keyword.compare(0, keyword.size(), buf, keyword.size()) != 0))
    {
      is.getline(buf, bufsize);
    }

    libMesh::NumericVector<Number>& solution = get_solution_vector();

    for (size_t i = 0; i < solution.size(); ++i)
    {
      if (!is.good()) throw InitFailedException("Bad datafile (corrupted data?)");

      double val;
      is.read(buf, sizeof(double));
      val = *(reinterpret_cast<double*>(buf));

      solution.set(i, val);
    }
    keyword = "</data>";
    is.getline(buf, bufsize);
    while (is.good() && (keyword.compare(0, keyword.size(), buf, keyword.size()) != 0))
    {
      is.getline(buf, bufsize);
    }

    solution.close();

    if (_systems.size() > 0)
      get_equation_system<libMesh::System>(0).update();

    has_read = true;
  }

  if (has_read)
  {
    equilibrium_done(true);
    increment_solve_sequence_number();
  }
}





void
SimulationInterface::do_plot_old(void)
{

  if (!has_environment()) return;
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

  DataOutput data_output(get_mesh(), formatstr);
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

  map<ID, libMesh::NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, libMesh::NumericVector<double>*>::iterator it(_remembered_solutions.find(id));

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

  map<ID, libMesh::NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, libMesh::NumericVector<double>*>::iterator it(_remembered_solutions.find(id));

  if (it != end)
  {
    get_solution_vector() = *(it->second);
    get_solution_vector().close();

    if (_systems.size() > 0)
      get_equation_system<libMesh::System>(0).update();
  }
}



void
SimulationInterface::do_delete_remembered_solution(ID id)
{

  map<ID, libMesh::NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, libMesh::NumericVector<double>*>::iterator it(_remembered_solutions.find(id));
  if (it != end)
  {
    delete it->second;
    _remembered_solutions.erase(it);
  }
}



libMesh::NumericVector<double>*
SimulationInterface::get_remembered_solution(ID id)
{
  libMesh::NumericVector<double>* vec = nullptr;

  map<ID, libMesh::NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, libMesh::NumericVector<double>*>::iterator it(_remembered_solutions.find(id));
  if (it != end)
    vec = it->second;

  return vec;
}


double
SimulationInterface::do_maximum_norm_of_difference(ID id)
{
  double norm = 0.0;

  map<ID, libMesh::NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, libMesh::NumericVector<double>*>::iterator it(_remembered_solutions.find(id));
  if (it != end)
  {
    libMesh::NumericVector<double>& old = *(it->second);
    libMesh::NumericVector<double>& current = get_solution_vector();

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

  map<ID, libMesh::NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, libMesh::NumericVector<double>*>::iterator it(_remembered_solutions.find(id));
  if (it != end)
  {
    libMesh::NumericVector<double>& old = *(it->second);
    libMesh::NumericVector<double>& current = get_solution_vector();

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
  {
    get_solution_vector().scale(factor);

    if (_systems.size() > 0)
      get_equation_system<libMesh::System>(0).update();
  }
}




void
SimulationInterface::do_add_scaled_remembered_solution(ID id, double factor)
{
  map<ID, libMesh::NumericVector<double>*>::iterator end(_remembered_solutions.end());
  map<ID, libMesh::NumericVector<double>*>::iterator it(_remembered_solutions.find(id));
  if (it != end)
  {
    get_solution_vector().add(factor, *(it->second));

    if (_systems.size() > 0)
      get_equation_system<libMesh::System>(0).update();
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

std::unique_ptr<libMesh::FEBase>
SimulationInterface::build_finite_element(unsigned int dim, libMesh::FEType type,
                                          bool scale) const
{

  double x0 = scale ? get_scaling().get_length_scaling() : 1.0;
  double mu = get_scaling().get_calc_mesh_units();
  TiberCad::Symmetry sym = get_symmetry();

  libMesh:: FEBase* fe;

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
      fe = nullptr;
  }
  assert(fe != nullptr);

  return std::unique_ptr<libMesh::FEBase>(fe);
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

    unsigned int n = get_mesh().n_active_elem();



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

    unsigned int n = get_mesh().n_nodes();
    if (results.size() != n * legend.size())
    {
      ostringstream s;
      s << "SimulationInterface::get_nodal_results: simulation "
        << this->get_name() << " returned broken nodal results.";
      throw runtime_error(s.str());
    }
  }
}


bool
SimulationInterface::get_solution(map<ID, vector<double> >& values,
  const vector<libMesh::Point>& p)
{
  return get_solution(nullptr, values, p, false);
}


bool
SimulationInterface::get_solution(const libMesh::Elem* elem,
    map<ID, vector<double> >& values,
    const vector<libMesh::Point>& p, bool local_coord)
{
  if (!is_solved()) return false;

  if (elem == nullptr)
  {
    if ((p.size() == 0) || local_coord)
      throw RuntimeException("Need at least one point in global coordinates"
         " in get_solution(Elem* elem, ...) with elem = nullptr.");

    elem = MeshUtils::search_element(&get_mesh(), p[0]);
  }

  if (elem == nullptr)
    throw RuntimeException("No element provided in get_solution(Elem*, ...).");
  
  vector<libMesh::Point> points(p);


  // if no points are given, we use the vertices
  unsigned int nn = points.size();
  if (nn == 0)
  {

    nn = elem->n_nodes();
    points.resize(nn);
    for (unsigned int i = 0; i < nn; i++)
    {
      points[i] = elem->master_point(i);
    }
    local_coord = true;
  }

  // first resize all data vectors to the right size
  map<ID, vector<double> >::iterator it(values.begin());
  const map<ID, vector<double> >::iterator end(values.end());
  for ( ; it != end; ++it)
  {
    const SolutionDescriptor& sd = get_solution_descriptor(it->first);
    //assert(sd.id() != INVALID_ID);
    unsigned int n_comp = sd.n_components();
    switch (sd.location())
    {
      case SolutionDescriptor::CELL:
      case SolutionDescriptor::NODES:
        values[it->first].resize(nn * n_comp, 0.0);
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


  // perhaps is in a quantum_contact ?
  QuantumContact* qc = get_environment().get_device().get_quantum_contact(elem);
  if (qc != nullptr)
  {
    if (local_coord)
    {
      for (unsigned int qp = 0; qp < points.size(); qp++)
      {
        // old style. In new libMESH this is in FEInterface
        switch (get_mesh().mesh_dimension())
        {
          case 1:
            points[qp] = libMesh::FE<1, libMesh::LAGRANGE>::map(elem, points[qp]);
            break;
          case 2:
            points[qp] = libMesh::FE<2, libMesh::LAGRANGE>::map(elem, points[qp]);
            break;
          case 3:
            points[qp] = libMesh::FE<3, libMesh::LAGRANGE>::map(elem, points[qp]);
        }
      }
    }

    pair<const Elem*, vector<Point>> pair = qc->project_on_boundary(elem, points);

    elem = pair.first;

    libMesh::FEInterface::inverse_map(get_mesh().mesh_dimension(), libMesh::FEType(), elem, pair.second, points);
   
    local_coord = true;
  }


  if (!local_coord)
  {
    libMesh::FEInterface::inverse_map(get_mesh().mesh_dimension(), libMesh::FEType(), elem, p, points);
  }

  SimulationEnvironment& env = get_environment();

  const libMesh::Elem* el = elem;

  bool flag = true;


  map<ID, vector<double> > new_values(values);
  set<unsigned int> req_points_id;
  vector<Point> req_points;
  req_points.reserve(points.size());

  for (it = values.begin(); it != end; ++it)
  {
    vector<double> cached_data;
    set<unsigned int> found = _data_cache.get_data(el, points, it->first, cached_data);

    unsigned int nd = found.size();
    // size of single dataset
    nd = nd > 0 ? cached_data.size() / nd : 0;

    unsigned int ctr = 0;
    for (unsigned int i = 0; i < points.size(); ++i)
    {
      if (found.count(i))
      {
        for (unsigned int j = 0; j < nd; ++j)
          values[it->first][nd*i+j] = cached_data[nd*ctr+j];

        ++ctr;
      }
      else // add missing points to the required points array
        req_points_id.insert(i);
    }
  }

  // prepare vector with required coordinates, ordered as the set with IDs
  for (auto mit(req_points_id.begin()); mit != req_points_id.end(); ++mit)
  {
    req_points.push_back(points[*mit]);
  }

  // resize data container accordingly
  for (auto it(new_values.begin()); it != new_values.end(); ++it)
  {
    it->second.resize(req_points.size() * it->second.size() / nn);
  }

  if (req_points.size() > 0)
  {
    if (!env.contains_element(elem))
    {

      // perhaps the parent is?
      const Elem* parent = elem->parent();

      while ((parent != nullptr) && (!env.contains_element(parent)))
        parent = parent->parent();

      el = parent; // is nullptr if no parent

      if (el != nullptr)
      {
        get_solution_secure(el, new_values, req_points);
      }
      else
      {
        MeshUtils::GridMapper& mapper =
            MeshUtils::GridMapper::get_mapper(get_mesh(), get_region_ids());

        vector<Point> req_points_r;
        vector<Point> req_points_r_loc;
        req_points_r.reserve(points.size());
        for (auto mit(req_points_id.begin()); mit != req_points_id.end(); ++mit)
        {
          req_points_r.push_back(p[*mit]);
          req_points_r_loc.push_back(points[*mit]);
        }


        vector<const Elem*> elems;
        for (unsigned int np = 0; np < req_points_r.size(); ++np)
        {
          //const Elem* my_el = MeshUtils::search_element(&get_mesh(), req_points_r[np]);
          const Elem* my_el = mapper.get_element(req_points_r[np]);
          elems.push_back(my_el);
        }

        for (unsigned int nel = 0; nel < elems.size(); ++nel)
        {
          map<ID, vector<double>> elem_values(new_values);

          for (auto nit(elem_values.begin()); nit != elem_values.end(); ++nit)
          {
            int ncomp = nit->second.size() / req_points.size();
            nit->second.resize(0);
            nit->second.resize(ncomp, 0);
          }

          if (elems[nel] != nullptr)
          {
            vector<Point> my_p(1, libMesh::FEInterface::inverse_map(
                get_mesh().mesh_dimension(), libMesh::FEType(), elems[nel], req_points_r[nel]));

            get_solution_secure(elems[nel], elem_values, my_p);
          }
          else
          {
            // TODO we have foreseen only one flag for the whole element
            // to indicate whether data is available. For now, we invalidate
            // all data on the element, if not all points can be obtained
            flag = false;
          }

          for (auto it(new_values.begin()); it != new_values.end(); ++it)
          {
            int ncomp = it->second.size() / req_points.size();

            for (unsigned int nc = 0; nc < ncomp; ++nc)
            {
              it->second[nel * ncomp + nc] = elem_values[it->first][nc];
            }

          }
        }
      }
    }
    else
      get_solution_secure(elem, new_values, req_points);


    if (_use_cache)
    {
      _data_cache.put_data(elem, req_points, new_values);
    }

    for (it = new_values.begin(); it != new_values.end(); ++it)
    {
      unsigned int nd = (it->second).size() / req_points.size();

      unsigned int ctr = 0;

      // add new data
      for (auto&& p_id : req_points_id)
      {
        for (unsigned int j = 0; j < nd; ++j)
          values[it->first][nd*p_id+j] = new_values[it->first][nd*ctr+j];

        ++ctr;
      }
    }
  }


  // cell based solutions are copied to all requested points
  it = values.begin();
  for ( ; it != end; ++it)
  {
    const SolutionDescriptor& sd = get_solution_descriptor(it->first);
    unsigned int n_comp = sd.n_components();
    if (sd.location() == SolutionDescriptor::CELL)
    {
      vector<double>& vec = values[it->first];
      for (size_t i = 1; i < nn; ++i)
      {
        size_t shift = i * n_comp;
        for (size_t c = 0; c < n_comp; ++c)
          vec[shift + c] = vec[c];
      }
    }
  }

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
SimulationInterface::get_solution(const Atom* atom, map<ID, vector<double> >& values)
{
  if (!is_solved()) return false;

  // Two possibilities:
  //  1) we have an atomistic structure, therefore try to get results on atoms
  //  2) we have only a mesh, therefore get the result on the atomic position

  bool ret = false;

  if (_atomistic_structure != nullptr)
  {
    get_solution_secure(atom, values);
    ret = true;
  }
  else
  {
    vector<Point> p(1, atom->get_position());
    ret = get_solution(atom->get_elem(), values, p);
  }

  return ret;
}


bool
SimulationInterface::get_solution(map<ID, vector<double> >& values)
{
  if (!is_solved())
    return(false);

  if (values.size() == 0)
  {
    IDSet::iterator it(_plotvariable_ids.begin());
    const IDSet::iterator end(_plotvariable_ids.end());
    for ( ; it != end; ++it)
      if (get_solution_descriptor(*it).location() == SolutionDescriptor::GLOBAL)
        values[*it] = vector<double>(0);
  }

  // first resize all data vectors to the right size
  map<ID, vector<double> >::iterator it(values.begin());
  map<ID, vector<double> >::iterator end(values.end());
  for ( ; it != end; ++it)
  {
    const SolutionDescriptor& sd = get_solution_descriptor(it->first);
    assert(sd.id() != INVALID_ID);
    unsigned int n_comp = sd.n_components();
    if (sd.location() != SolutionDescriptor::GLOBAL)
    {
      ostringstream os;
      os << "In simulation \'" << get_name() << "\' :"
          << "cannot access solution variable \'" << sd.name()
          << "\' as global (mesh-independent) quantity.";
      throw RuntimeException(os.str());
    }
    values[it->first].resize(n_comp);
  }

  get_solution_secure(values);

  return(true);
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


/*
void
SimulationInterface::get_solution_secure(const Elem*,
        const std::set<ID>&, std::vector<std::map<ID, double> >&)
{
}
*/
/*
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
*/



void
SimulationInterface::declare_solution_ext(const std::string& name, ID id,
    SolutionDescriptor::Type type, SolutionDescriptor::Location location,
    const std::string& units, unsigned int n_comp)
{
  _solution_descriptors[id] = SolutionDescriptor(name, id, type, location,
          units, n_comp);
  _solution_ids[name] = id;



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


SolutionDescriptor&
SimulationInterface::get_solution_descriptor(ID id)
{
  return(const_cast<SolutionDescriptor&>(
      static_cast<const SimulationInterface&>(*this).get_solution_descriptor(id)));
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


SolutionDescriptor&
SimulationInterface::get_solution_descriptor(const std::string& solution_name)
{
  return(const_cast<SolutionDescriptor&>(
      static_cast<const SimulationInterface&>(*this).get_solution_descriptor(solution_name)));
}

const SolutionDescriptor&
SimulationInterface::get_solution_descriptor(const std::string& solution_name) const
{
  ID id = INVALID_ID;

  // first split any strain:xx:yy:zz
  vector<string> tokens;
  Utils::tokenize(solution_name, tokens, ":");

  if (tokens.size() > 0)
  {
    // the first token is the solution name
    map<const string, ID>::const_iterator it(_solution_ids.find(tokens[0]));
    if (it != _solution_ids.end())
      id = it->second;
  }

  return get_solution_descriptor(id);
}


ID
SimulationInterface::convert_variable_name_to_id(
    const std::string& variable_name) const
{
  return get_solution_descriptor(variable_name).id();
}



ModelOptions&
SimulationInterface::get_solver_options(const std::string& block)
{
  if (!get_options().has_submodel(block))
    get_options().add_submodel(block, ModelOptions());

  ModelOptions::submodel_iterator it(get_options().submodels_begin(block));

  return it->second;
}



Embracing*
SimulationInterface::create_embracing_region(
    SimulationInterface* other_simulation,
    const ModelOptions& options, bool need_mixing_coeff)
{
  Embracing* emb = nullptr;
  if (other_simulation != nullptr)
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
    if (_solution_descriptors.size() > 1)
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
    }


    Messages::newline();
    ostringstream os;
    os << "Output directory      : " << get_output_directory() << endl;
    os << "Output file format    : " << get_option("output_format", "-") << endl;
    os << "Output file basename  : " << get_output_filename_prefix() << endl;
    Messages::info(os.str());
    Messages::newline();

    get_option("binary_output", true);
  }

  if (verbose() > 0) do_print_info();
}




double
SimulationInterface::build_map_elem_atoms(double sigma, double cutoff)
{
  
  // Get total number of elements
  // (the map is oversized, but faster since the elem ID is used as key)
  _elem_to_atoms.resize(get_mesh().n_elem());

  double scale = get_atomistic_structure()->get_scale();

  // Maximum cutoff distance
  //const double tau = 1.0 / projection_length; // projection in Angstroms
  //const double deltar_max = (5.0*log(10.0) - log( 8.0*3.141593/(tau*tau*tau) )  ) / tau;   

  const double sigma2 = 2.0*sigma*sigma;
  const double deltar2_max = -sigma2*(log(cutoff) + 1.5*log(2.0*3.141593*sigma)); 
  const double deltar_max = sqrt(deltar2_max);
  
  unsigned int N_wo_H = get_atomistic_structure()->get_N_without_H();
  // Estimate number of atoms in a sphere
  unsigned int Nat = round( sqrt(3.0)*3.1416/2.0 * pow(deltar_max/1.90, 3.0) );

  //std::cout<<"Scale: "<<scale<<std::endl;
  //std::cout<<"Sigma: "<<sigma<<std::endl;
  //std::cout<<"Rmax: "<<deltar_max<<std::endl;
  //std::cout<<"Atoms in sphere Rmax: "<<Nat<<std::endl;

  const std::vector<Atom>& structure = get_atomistic_structure()->get_structure_atoms();

  // Partition the structure into cells for O(N) scheme
  Tensor2 period = get_atomistic_structure()->get_ttype_lattice_vectors(); 
  GridCells cells(structure, period,
      get_atomistic_structure()->get_origin(),
      deltar_max, get_mesh().mesh_dimension());
  //cells.print_statistics();

  unsigned int notassociated = 0;

  MeshBase::const_element_iterator elit = this->active_local_elements_begin();
  const MeshBase::const_element_iterator elend = this->active_local_elements_end();
  
  for ( ; elit != elend; ++elit)
  { 
    const Elem* elem = *elit;
 
    Point pc = elem->vertex_average() * scale;
    //double x=pc(0), y=pc(1), z=pc(2);

    std::vector<unsigned int> temp;
    temp.reserve(Nat);

    // Find the cell containing the element centroid
    unsigned int l, m, n;
    cells.get_cell(pc, l, m, n); 

    // Loop on all 27 neighboring cells (periodicity is taken care by the iterator)    
    GridCells::NeighborIterator it = cells.begin(l,m,n);
    GridCells::NeighborIterator end = cells.end(l,m,n);
    unsigned int count = 0;

    for ( ; it != end; ++it)
    {
      unsigned int c1 = (*it).first;
      const Tensor1& shift = *((*it).second);
   
      //unsigned int u,v,w; cells.index(c1,u,v,w);
      //std::cout<<"cell: "<<u<<" "<<v<<" "<<w<<" natoms: "<<cells[c1].size()<<std::endl;
      // Loop over all atoms in each cell
      for (unsigned int i = 0; i < cells[c1].size(); i++)
      {     
        unsigned int iatm = cells[c1][i];

        if (structure[iatm].get_elem() == nullptr ) continue;

        Point delta_r(structure[iatm].get_position() - pc);
        
        // the following is quite a hack, just to make it take all atoms in the
        // directions orthogonal to the simulation domain
        switch (get_mesh().mesh_dimension())
        {
          case 1:
            delta_r(1) = 0.0;

          case 2:
            delta_r(2) = 0.0;

          default:
            break;
        }
        
        // do we need the checks in single directions?

        if (delta_r * delta_r > deltar2_max) continue;


        temp.push_back(iatm);
        count++;
      }    
    }


    unsigned int id = elem->id();

    _elem_to_atoms[id].resize(count);

    for (unsigned int iatm = 0; iatm  <  count; iatm++)
    	_elem_to_atoms[id][iatm] = temp[iatm];

    temp.clear();
  }
  
  //std::cerr<<"Map done"<<std::endl;

  return deltar_max;

}




void
SimulationInterface::project_on_tensor_grid(void)
{
  if (!get_options().has_submodel("Projection"))
    return;

  if (this->get_communicator().rank() != 0)
    return;

  Messages msg;
  msg.info("Project solutions on 2D tensor grid");
  msg.indent();

  const ModelOptions& opts = (get_options().submodels_begin("Projection"))->second;

  if (get_mesh().mesh_dimension() < 2)
  {
    //Messages::warning("Projection is implemented only for 2D meshes.");
    Messages::warning("Projection is implemented only for 2D and 3D meshes.");
    return;
  }

  string format = opts.get_option("format", "ascii");
  if (format != "ascii")
  {
    Messages::warning("Projection is implemented only for ascii format.");
    return;
  }

  IDSet sol_ids;
  vector<string> solutionnames;
  opts.get_option("solutions", solutionnames);
  for (int i = 0; i < solutionnames.size(); ++i)
  {
    sol_ids.insert(get_solution_id(solutionnames[i]));
  }

  if (sol_ids.empty())
    sol_ids = get_plotvariable_ids();



  // get the bounding box

  auto bbox(get_environment().get_bounding_box(true));
  Point& pmax = bbox.second;
  Point& pmin = bbox.first;

  ostringstream os;
  os << "Bounding box: (" << pmin(0) << ", " << pmin(1) << ", " << pmin(2)
      << ") - (" << pmax(0) << ", " << pmax(1) << ", " << pmax(2) << ")\n";
  msg.info(os.str());
  pmin(0) += 1e-11;
  pmin(1) += 1e-11;
  pmin(2) += 1e-11;
  pmax(0) -= 1e-11;
  pmax(1) -= 1e-11;
  pmax(2) -= 1e-11;

  // get the projection plane
  string pplane_s = opts.get_option("projection_plane", "xy");
  int pplane = 0; // xy

  if (pplane_s == "xz")
  {
    pplane = 1;
  }
  else if (pplane_s == "yz")
  {
    pplane = 2;
  }

  double dx = opts.get_option("spacing", get_mesh_units());
  dx /= get_mesh_units();

  int Nx = ceil((pmax(0) - pmin(0)) / dx);
  int Ny = ceil((pmax(1) - pmin(1)) / dx);

  switch (pplane)
  {
    case 1:
      Nx = ceil((pmax(0) - pmin(0)) / dx);
      Ny = ceil((pmax(2) - pmin(2)) / dx);
      pmin(1) = 0.5 * (pmax(1) + pmin(1));
      break;

    case 2:
      Nx = ceil((pmax(1) - pmin(1)) / dx);
      Ny = ceil((pmax(2) - pmin(2)) / dx);
      pmin(0) = 0.5 * (pmax(0) + pmin(0));
      break;

    default:
      pmin(2) = 0.5 * (pmax(2) + pmin(2));
      break;
  }

  // we allow two formats:
  // 0:
  //    x y data1 data2 ...
  //
  // 1: (one data per file)
  //   data1(x1y1) data1(x2y1) ...
  //   data1(x1y2) data1(x2y2) ...
  int ascii_form = 1;
  if (opts.get_option("ascii_format", "list") != string("list"))
    ascii_form = 0;



  string basename = get_output_directory() + "/" + get_output_filename() + "_projected";
  string extension = ".dat";

  // to get solutions
  map<ID, vector<double>> solutions;

  // the file streams
  map<ID, vector<ofstream*>> fstreams;

  if (ascii_form == 0)
  {
    // in this case we use one single file
    fstreams[0] = vector<ofstream*>(1);
    fstreams[0][0] = new ofstream(basename + extension);
    (*fstreams[0][0]) << "% pmin = (" << pmin(0) << ", " << pmin(1) << ", " << pmin(2) <<
            "), dx = " << dx*get_mesh_units() << " Nx = " << Nx << " Ny = " << Ny << endl;
    (*fstreams[0][0]) << "% ";
    switch (pplane)
    {
      case 1:
        (*fstreams[0][0]) << "x z ";
        break;

      case 2:
        (*fstreams[0][0]) << "y z";
        break;

      default:
        (*fstreams[0][0]) << "x y";
        break;
    }
  }
  else
  {
    basename += "_";
  }

  IDSet::const_iterator idit(sol_ids.begin());
  IDSet::const_iterator idend(sol_ids.end());

  for ( ; idit != idend; ++idit)
  {
    const SolutionDescriptor& descr = get_solution_descriptor(*idit);
    if (descr.on_mesh())
    {
      solutions[*idit].clear();
      int ncomp = descr.n_components();
      fstreams[*idit].resize(ncomp, nullptr);

      vector<string> comp(ncomp, "");
      if (ncomp == 3)
      {
        comp[0] = "_x";
        comp[1] = "_y";
        comp[2] = "_z";
      }
      else if (ncomp == 6)
      {
        comp[0] = "_xx";
        comp[1] = "_yy";
        comp[2] = "_zz";
        comp[3] = "_xy";
        comp[4] = "_yz";
        comp[5] = "_xz";
      }
      else
      {
        for (int i = 0; i < ncomp; ++i)
        {
          ostringstream os;
          os << "_" << i;
          comp[i] = os.str();
        }
      }

      for (int i = 0; i < ncomp; ++i)
      {
        if (ascii_form == 1)
        {
          fstreams[*idit][i] = new ofstream(basename + descr.name() + comp[i] + extension);
          (*fstreams[*idit][i]) << "% ";
          (*fstreams[*idit][i]) << get_name() << " " << descr.name() << comp[i] <<
            "\n% origin = (" << pmin(0) << ", " << pmin(1) << ", " << pmin(2) <<
            "), dx = " << dx*get_mesh_units() << " Nx = " << Nx << " Ny = " << Ny << endl;
        }
        else
        {
          (*fstreams[0][0]) << " " << descr.name() << comp[i];
        }
      }
    }
  }

  if (ascii_form == 0)
    (*fstreams[0][0]) << "\n";


  map<ID, vector<double>>::iterator mit(solutions.begin());
  map<ID, vector<double>>::iterator mend(solutions.end());

  for (unsigned int i = 0; i < Nx; ++i)
  {
    for (unsigned int j = 0; j < Ny; ++j)
    {
      Point p(pmin(0) + i * dx, pmin(1) + j * dx, pmin(2));
      if (pplane == 1)
        p = Point(pmin(0) + i * dx, pmin(1), pmin(2) + j * dx);
      if (pplane == 2)
        p = Point(pmin(0), pmin(1) + i * dx, pmin(2) + j * dx);

      if (ascii_form == 0)
      {
        switch (pplane)
        {
          case 1:
            (*fstreams[0][0]) << p(0) << " " << p(2) << " ";
            break;

          case 2:
            (*fstreams[0][0]) << p(1) << " " << p(2) << " ";
            break;

          default:
            (*fstreams[0][0]) << p(0) << " " << p(1) << " ";
            break;
        }
      }

      const Elem* elem = MeshUtils::search_element(&get_mesh(), p);

      if ((elem == nullptr) || !get_solution(elem, solutions, vector<Point>(1, p)))
      {
        for (mit = solutions.begin(); mit != mend; ++mit)
        {
          int ncomp = get_solution_descriptor(mit->first).n_components();
          (mit->second).clear();
          (mit->second).resize(ncomp, 0.0);
        }
      }

      for (mit = solutions.begin(); mit != mend; ++mit)
      {
        int ncomp = get_solution_descriptor(mit->first).n_components();
        for (unsigned int c = 0; c < ncomp; ++c)
        {
          if (ascii_form == 1)
          {
            (*fstreams[mit->first][c]) << (mit->second)[c] << " ";
          }
          else
          {
            (*fstreams[0][0]) << (mit->second)[c] << " ";
          }
        }
      }

      if (ascii_form == 0)
        (*fstreams[0][0]) << "\n";
    }

    if (ascii_form == 1)
    {
      for (mit = solutions.begin(); mit != mend; ++mit)
      {
        int ncomp = get_solution_descriptor(mit->first).n_components();
        for (unsigned int c = 0; c < ncomp; ++c)
        {
          (*fstreams[mit->first][c]) << endl;
        }
      }
    }
  }


  if (ascii_form == 1)
  {
    for ( ; idit != idend; ++idit)
    {
      const SolutionDescriptor& descr = get_solution_descriptor(*idit);
      if (descr.on_mesh())
      {
        for (int i = 0; i < descr.n_components(); ++i)
        {
          fstreams[*idit][i]->flush();
          fstreams[*idit][i]->close();
          delete fstreams[*idit][i];
        }
      }
    }
  }
  else
  {
    fstreams[0][0]->flush();
    fstreams[0][0]->close();
    delete fstreams[0][0];
  }
}
