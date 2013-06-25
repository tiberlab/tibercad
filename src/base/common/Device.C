// $Id$

#include "Device.h"
#include "Material.h"
#include "MeshUtils.h"
#include "MeshReader.h"
#include "DataOutput.h"
#include "ReadGMSH.h"
#include "ReadISEGrid.h"
#include "BoundaryRegions.h"
#include "MaterialBoundary.h"
#include "EdgeObject.h"
#include "NodeObject.h"
#include "MeshRegionInfo.h"
#include "SimulationOptions.h"
#include "AtomisticStructure.h"
#include "QuantumContact.h"

#include "gmsh_io.h"
#include "equation_systems.h"
#include "mesh.h"
#include "unstructured_mesh.h"
#include "elem.h"
#include "Messages.h"

#include <iostream>
#include <memory>


using namespace std;


namespace
{
  class IDPair
  {
    public:
      IDPair(ID a, ID b) : _a(a), _b(b)
      { if (b < a) {_b = a; _a = b;} }

      bool operator==(const IDPair& rhs) const
          { return ((_a == rhs._a) && (_b == rhs._b)); }
      bool operator!=(const IDPair& rhs) const
          { return !(*this == rhs); }
      bool operator<(const IDPair& rhs) const
          { return (_a < rhs._a) || (!(rhs._a < _a) && (_b < rhs._b)); }

    private:
      ID _a, _b;
  };
}

Device::Device(void)
  : _mesh(NULL),
    _mesh_units(1e-6),
    _eq_system(NULL),
    _boundary_nodes(NULL),
    _symmetry(TiberCad::NONE)
{
  _material_map.clear();
}


Device::~Device(void)
{
  // we put them first into a set because a material can be associated
  // to several IDs
  MaterialMap::iterator it(_material_map.begin());
  const MaterialMap::iterator end(_material_map.end());
  set<Material*> mats;
  for ( ; it != end; ++it)
    mats.insert(it->second);

  set<Material*>::iterator matit(mats.begin());
  const set<Material*>::iterator matend(mats.end());
  for ( ; matit != matend; ++matit)
    delete *matit;

  // every ID has its own object
  BoundaryMap::iterator bdit(_boundary_map.begin());
  const BoundaryMap::iterator bdend(_boundary_map.end());
  for ( ; bdit != bdend; ++bdit)
    delete bdit->second;

  // every ID has its own object
  EdgeObjMap::iterator eit(_edge_map.begin());
  const EdgeObjMap::iterator eend(_edge_map.end());
  for ( ; eit != eend; ++eit)
    delete eit->second;

  // every ID has its own object
  NodeObjMap::iterator nit(_node_map.begin());
  const NodeObjMap::iterator nend(_node_map.end());
  for ( ; nit != nend; ++nit)
    delete nit->second;

  _material_map.clear();

  delete _eq_system;
  delete _boundary_nodes;
  delete _mesh;
  delete _mesh_region_info;
  delete _bd_regions;
}


Device*
Device::create(const ModelOptions& options)
{

  Device* device = new Device();

  device->set_options(options);

  return device;
}


void
Device::prepare(void)
{
  setup_mesh();
  setup_regions();
  setup_clusters();
  setup_quantum_contacts();
  setup_atomistic_structures();

  // for now, this works only for D > 1
  if (_options.get_option("write_boundary_mesh", false) &&
      (get_mesh().mesh_dimension() > 1))
  {
    auto_ptr<DataOutput> writer(DataOutput::create("vtk"));
    if (writer.get() != NULL)
    {
      writer->set_output_directory(_options.get_option("output_path", "./"));
      writer->set_filename(Utils::basename(_options["meshfile"]) + "_bnd");

      AutoPtr<MeshBase> bdmesh = MeshUtils::create_boundary_mesh(get_mesh());

      writer->set_mesh(*bdmesh);
      writer->write(true);
    }
  }
}


void
Device::setup_mesh(void)
{

  delete _mesh;

  Messages m;
  m.info("Setup mesh");
  m.indent();

  _mesh_units = _options.get_option("mesh_units", _mesh_units);

  // this is backup solution if dim cannot be guessed from the mesh file
  int dim = _options.get_option("dimension", 1);

  const string& meshfile = _options["meshfile"];

  {
    ostringstream os;
    os << "Reading mesh file " << meshfile << " ...";
    m.info(os.str(), false);
  }

  _mesh = new Mesh(dim);

  _mesh_region_info = new MeshRegionInfo;
  _bd_regions = new BoundaryRegions;

  MeshReader::read_mesh(meshfile, *_mesh, *_mesh_region_info, *_bd_regions);


  // TODO only for now (back compatibility)
  delete _boundary_nodes;
  _boundary_nodes = new map<unsigned int, vector<ID> >();
  _bd_regions->get_bc_node_map(*_boundary_nodes);


  // update mesh dimension
  dim = _mesh->mesh_dimension();

  m.info(" done.");
  {
    ostringstream os;

    os << "mesh units          : " << _mesh_units << " m" << endl
       << "mesh dimension      : " << setw(7) << setfill(' ') << dim << endl
       << "number of nodes     : " << setw(7) << setfill(' ') << _mesh->n_nodes() << endl
       << "number of elements  : " << setw(7) << setfill(' ') << _mesh->n_elem() << endl
       << "number of subdomains: " << setw(7) << setfill(' ') << _mesh_region_info->n_subdomains();
    m.info(os.str());
  }
  m.newline();

  //_mesh_region_info->print_info();

  const string& sym = _options.get_option("symmetry", "");
  if (sym == "cylindrical")
  {
    _symmetry = TiberCad::CYLINDRICAL;
    m.info("Using cylinder symmetry (=> 3D simulation)");
  }


  /*
   * NOTE:
   * In parallel case, the local mesh does not contain all elements
   * and therefore only a part of the region IDs are present
   */


  _eq_system = new EquationSystems(*_mesh);
}





void
Device::setup_regions(void)
{

  Messages::debug("Control::create_materials() begin");

  Messages m;
  m.info("Create materials ...");
  m.indent();

  //
  // first we process the physical regions
  //

  // iterate the regions and create the materials
  ModelOptions::const_submodel_iterator rit(_options.submodels_begin("Region"));
  const ModelOptions::const_submodel_iterator rend(_options.submodels_end("Region"));
  for ( ; rit != rend; ++rit)
  {
    ModelOptions data(rit->second);

    // we read the region numbers as strings as they could be region names
    vector<string> region_ids_str;
    string reg("");
    reg = data.get_option("mesh_regions", reg);
    Utils::extract_vector(reg, region_ids_str);

    // for the numeric region IDs
    vector<ID> region_ids;

    unsigned int n_ids = region_ids_str.size();
    // if no numbers are specified we try to get them from the region name
    if (n_ids == 0)
      get_mesh_region_ids(data.get_name(), region_ids);
    else
    {
      vector<ID> tmp_id;
      for (unsigned int i = 0; i < n_ids; i++)
      {
        get_mesh_region_ids(region_ids_str[i], tmp_id);
        if (tmp_id.size() == 0)
        {
          ostringstream s;
          s << "Physical region \'" << region_ids_str[i]
            << "\' (in Region \'" <<  data.get_name()
            << "\') does not exist in mesh.";
          throw InitFailedException(s.str());
        }
        region_ids.insert(region_ids.end(), tmp_id.begin(), tmp_id.end());
      }
    }


    if (region_ids.size() == 0)
    {
      ostringstream s;
      s << "Physical region \'" << data.get_name() <<
        "\' is not consistent with mesh.";
      throw InitFailedException(s.str());
    }

    // some options can be provided globally for all regions

    // The default material is Si
    string material = _options.get_option("material", "Si");
    string x_frac = _options.get_option("x", "0.5");
    string xdir = _options.get_option("x-growth-direction", "");
    string ydir = _options.get_option("y-growth-direction", "");
    string zdir = _options.get_option("z-growth-direction", "");


    material = data.get_option("material", material);
    data["material"] = material;
    x_frac = data.get_option("x", x_frac);
    data["x"] = x_frac;
    xdir = data.get_option("x-growth-direction", xdir);
    if (!xdir.empty())
      data["x-growth-direction"] = xdir;
    ydir = data.get_option("y-growth-direction", ydir);
    if (!ydir.empty())
      data["y-growth-direction"] = ydir;
    zdir = data.get_option("z-growth-direction", zdir);
    if (!zdir.empty())
      data["z-growth-direction"] = zdir;

    data.set_option("dimension", get_mesh().mesh_dimension());

    Material* mat = Material::create(material, data);
    set_material(mat, region_ids, data.get_name());
  }

  m.unindent();
  m.info("Creation of materials done.");

  Messages::debug("Control::create_materials() end");
}




void
Device::setup_clusters(void)
{

  ModelOptions::const_submodel_iterator it(_options.submodels_begin("Cluster"));
  const ModelOptions::const_submodel_iterator end(_options.submodels_end("Cluster"));
  for ( ; it != end; ++it)
  {
    const ModelOptions& data = it->second;

    // we read the region numbers as strings as they could be region names
    string region_ids_str("");
    region_ids_str = data.get_option("regions", region_ids_str);

    // for the numeric region IDs
    set<ID> ids;
    extract_physical_regions(region_ids_str, ids);
    vector<ID> region_ids;
    region_ids.reserve(ids.size());
    region_ids.insert(region_ids.begin(), ids.begin(), ids.end());


    if (region_ids.size() > 0)
    {
      ostringstream os;
      os << "Setting up Cluster \'" << data.get_name()
        << "\' containing regions " << region_ids[0];
      for (size_t i = 1; i < region_ids.size(); i++)
        os << ", " << region_ids[i];
      Messages::newline();
      Messages::info(os.str());

      set_cluster(data.get_name(), region_ids);
    }
    else
      Messages::warning("Cluster \'" + data.get_name() + "\' is empty.");
  }
}



void
Device::setup_atomistic_structures(void)
{

  Messages::debug("Control::create_atomistic_structures() begin");

  ModelOptions::const_submodel_iterator it(_options.submodels_begin("Atomistic"));
  const ModelOptions::const_submodel_iterator end(_options.submodels_end("Atomistic"));
  for ( ; it != end; ++it)
  {
    const ModelOptions& data = it->second;

    const string& st_name = data.get_name();

    AtomisticStructure* st = AtomisticStructure::create();

    //WARNING: For debugging purposes, initialization of
    //atomistic structures is here, but it's not the right place! (maybe it is...)
    st->init(st_name, this, data);

    // Defined atomistic structure is put in the atomistic_structure_map
    _atomistic_structure_map[st_name] = st;

    //UnstructuredMesh* mesh = new Mesh(3);
    //st->create_conformal_grid(*mesh);
    //DataOutput* dto = DataOutput::create("vtk");
    //dto->set_filename("pippo");
    //dto->set_output_directory("./");
    //dto->set_mesh(*mesh);
    //dto->write(true);
  }

  Messages::debug("Control::create_atomistic_structures() end");
}

void
Device::setup_quantum_contacts(void)
{

  Messages::debug("Control::create_quantum_contacts() begin");

  ModelOptions::const_submodel_iterator it(_options.submodels_begin("Quantum_contact"));
  const ModelOptions::const_submodel_iterator end(_options.submodels_end("Quantum_contact"));

  for ( ; it != end; ++it)
  {

    const ModelOptions& data = it->second;

    const string& name = data.get_name();

    string regions("");
    regions = data.get_option("regions",regions);

    double length = data.get_option("length", 0.5);

    Messages::info(" ");
    Messages::info("Creating quantum contact "+name);

    std::vector<ID> rg_ids, bd_ids;
    std::set<ID> ids;

    // We get mesh regions (not active regions).

    extract_physical_regions(regions, ids);

    rg_ids = QuantumContact::set2vec(ids);

    //for (int n=0; n<rg_ids.size();n++)
    //std::cout << "rg_ids: "<< rg_ids[n]<<std::endl;

    if (rg_ids.size()==0)
      throw InitFailedException("In Quantum_contact user must define a valid mesh region");

    get_boundary_region_ids(name, bd_ids);
    //for (int n=0; n<bd_ids.size();n++)
    //std::cout << "bd_ids: "<< bd_ids[n]<<std::endl;

    ID  newid = _mesh_region_info->next_id();

    QuantumContact* st = QuantumContact::create();

    st->init(newid, name, this, _bd_regions, rg_ids, bd_ids, length);

    // Defined quantum contact is put in the quantum_contact_map
    _quantum_contact_map[name] = st;

    _mesh_region_info->add_id(newid);
    _mesh_region_info->set_name(newid, name);

    std::vector<ID> vid(1, newid);

    // Only one material/region can be created !!
    set_material(get_material(rg_ids[0]), vid, name);

    // We have to erase the region id from the list of active regions, otherwise
    // we mess up the other modules (quantum contact regions should be invisible
    // when asking for real device regions)
    _active_region_ids.erase(newid);

  }

  if (_options.has_submodel("Quantum_contact") )
  {

    const string& meshfile = _options["meshfile"];

    std::cout<<"write meshfile new_"+meshfile<<std::endl;

    GmshIO(*_mesh).write("new_" + meshfile);

   }

  Messages::debug("Control::create_quantum_contacts() end");
}


void
Device::extract_physical_regions(const std::string& str, IDSet& ids) const
{
  ids.clear();

  // the IDs that have to be excluded ("-pippo" syntax)
  IDSet exclude;

  // we have to get it as vector (for the moment at least)
  vector<string> preg;
  Utils::extract_vector(str, preg);

  set<ID> preg_ids;

  unsigned int n = preg.size();
  for (unsigned int i = 0; i < n; i++)
  {
    if (preg[i].at(0) == '-')
      get_active_region_ids(preg[i].substr(1), preg_ids);
    else
      get_active_region_ids(preg[i], preg_ids);

    if (preg_ids.size() == 0)
    {
      ostringstream s;
      s << "Physical region " << preg[i] <<
      " does not exist in mesh file.";
      throw InitFailedException(s.str());
    }

    if (preg[i].at(0) == '-')
      exclude.insert(preg_ids.begin(), preg_ids.end());
    else
      ids.insert(preg_ids.begin(), preg_ids.end());
  }

  if (ids.empty() && !exclude.empty())
    ids = get_active_region_ids();

  for (IDSet::iterator it(exclude.begin()); it != exclude.end(); ++it)
    ids.erase(*it);
}


/*
void
Device::prepare_boundaries(void)
{
  BCNodeMap::const_iterator bd_it;
  const BCNodeMap::const_iterator bd_end(_boundary_nodes->end());

  // we only look on level zero
  MeshBase::const_element_iterator el(get_mesh().level_elements_begin(0));
  const MeshBase::const_element_iterator el_end(get_mesh().level_elements_end(0));

  for ( ; el != el_end; ++el)
  {
    Elem* elem = *el;
    const ID id = elem->subdomain_id();

    // loop over the sides
    int n_sides = elem->n_sides();
    for (int s = 0; s < n_sides; s++)
    {

      // check if the neighbouring element is inexistent (outer boundary)
      // or in another simulation region (inner boundary)
      //
      // we allow inner 'boundaries', i.e. we don't really consider
      // boundaries but n-1 dimensional domains
      const Elem* neighbour = elem->neighbor(s);

      if ((neighbour == NULL) ||
          (neighbour->subdomain_id() != id))
      {

      }
    }
  }

}
*/


void
Device::init(void)
{

  Messages m;

  // init all materials
  m.info("Setting up bulk models ...");
  m.indent();
  MaterialMap::iterator it(_material_map.begin());
  const MaterialMap::iterator end(_material_map.end());
  for ( ; it != end; ++it)
    (it->second)->init();
  m.unindent();
  m.newline();

  // init all lower dimensional regions
  m.info("Setting up lower dimensional (boundary) models ...");
  m.indent();
  BoundaryMap::iterator bdit(_boundary_map.begin());
  const BoundaryMap::iterator bdend(_boundary_map.end());
  for ( ; bdit != bdend; ++bdit)
    (bdit->second)->init();

  EdgeObjMap::iterator eit(_edge_map.begin());
  const EdgeObjMap::iterator eend(_edge_map.end());
  for ( ; eit != eend; ++eit)
    (eit->second)->init();

  NodeObjMap::iterator nit(_node_map.begin());
  const NodeObjMap::iterator nend(_node_map.end());
  for ( ; nit != nend; ++nit)
    (nit->second)->init();
  m.unindent();

}



void
Device::get_atomistic_structures(const string& names,
    vector<AtomisticStructure*>& structures)
{
  if (names == "all")
  {
    structures.resize(_atomistic_structure_map.size());
    atomistic_structure_iterator it(_atomistic_structure_map.begin());
    for (size_t i = 0; it != _atomistic_structure_map.end(); ++it, ++i)
      structures[i] = it->second;
  }
  else
  {
    vector<string> ns;
    Utils::extract_vector(names, ns);
    structures.resize(ns.size());
    for (size_t i = 0; i < ns.size(); ++i)
      structures[i] = get_atomistic_structure(ns[i]);
  }
}


EquationSystems&
Device::get_equation_systems(MeshBase* mesh)
{
  EquationSystems* eqsys;
  if (mesh == _mesh)
    eqsys = _eq_system;
  else
  {
    map<const MeshBase*, EquationSystems*>::iterator it =
        _eq_sys_map.find(mesh);

    if (it == _eq_sys_map.end())
    {
      eqsys = new EquationSystems(*mesh);
      it = (_eq_sys_map.insert(make_pair(mesh, eqsys))).first;
    }

    eqsys = it->second;
  }

  return(*eqsys);
}



void
Device::set_material(Material* material, ID region_id)
{
  assert(material != NULL);

  if (!_mesh_region_info->has_id(region_id))
  {
    /*
     * In single processor case this has to be considered an error,
     * in parallel we should do the check in a different manner
     */
    if (libMesh::n_processors() == 1)
    {
      ostringstream s;
      s << "Device: region " << region_id <<
        " does not exist in mesh file.";
      throw InitFailedException(s.str());
    }
  }
  if (_material_map.find(region_id) != _material_map.end())
  {
    ostringstream s;
    s << "Device: trying to redefine mesh region " << region_id << ".";
    throw InitFailedException(s.str());
  }

  _material_map[region_id] = material;
  _active_region_ids.insert(region_id);
}



void
Device::set_material(Material* material, const vector<ID>& region_ids,
                     const string& region_name)
{
  assert(material != NULL);

  for (unsigned int i = 0; i < region_ids.size(); ++i)
    set_material(material, region_ids[i]);

  set_region_name(region_name, region_ids);

  ostringstream os;
  os << "Added material " << material->get_name()
     << " for region \'" << region_name
     << "\' (mesh regions " << region_ids[0];
  for (unsigned int i = 1; i < region_ids.size(); ++i)
    os << ", " << region_ids[i];
  os << ")" << endl
     << "  (using datafile " << material->get_database().get_data_file() << ")";
  Messages::info(os.str());
  material->info();
}





Material*
Device::get_material(const std::string& name)
{
  Material* mat = NULL;

  map<ID, string>::const_iterator it(_region_names.begin());
  const map<ID, string>::const_iterator end(_region_names.end());
  for ( ; it != end; ++it)
    if (it->second == name)
    {
      mat = get_material(it->first);
      break;
    }

  return mat;
}


const Material*
Device::get_material(const std::string& name) const
{
  const Material* mat = NULL;

  map<ID, string>::const_iterator it(_region_names.begin());
  const map<ID, string>::const_iterator end(_region_names.end());
  for ( ; it != end; ++it)
    if (it->second == name)
    {
      mat = get_material(it->first);
      break;
    }

  return mat;
}


MaterialBoundary*
Device::get_boundary_object(ID id)
{
  MaterialBoundary* mb = NULL;

  if (id != INVALID_ID)
  {
    BoundaryMap::iterator it(_boundary_map.find(id));
    if (it != _boundary_map.end())
      mb = it->second;
    else
    {
      if (_bd_regions->is_side(id))
      {
        // first get ids of the regions on both sides
        const IDSet& ids = _bd_regions->get_contiguous_regions_for_side(id);
        assert(ids.size() == 2);
        ID idA = *(ids.begin());
        ID idB = *(++(ids.begin()));
        Material* matA = get_material(idA);
        Material* matB = get_material(idB);
        ModelOptions opts;
        mb = MaterialBoundary::create(idA, matA, idB, matB, opts);
        _boundary_map[id] = mb;
      }
    }
  }

  return mb;
}



MaterialBoundary*
Device::get_boundary_object(const Elem* elem, int side)
{
  MaterialBoundary* mb = NULL;

  ID id = _bd_regions->get_side_id(elem, side);
  if (id != INVALID_ID)
  {
    BoundaryMap::iterator it(_boundary_map.find(id));
    if (it != _boundary_map.end())
      mb = it->second;
  }

  return mb;
}


EdgeObject*
Device::get_edge_object(ID id)
{
  EdgeObject* mb = NULL;

  if (id != INVALID_ID)
  {
    EdgeObjMap::iterator it(_edge_map.find(id));
    if (it != _edge_map.end())
      mb = it->second;
    else
    {
      if (_bd_regions->is_edge(id))
      {
        ModelOptions opts;
        mb = EdgeObject::create(opts);
        _edge_map[id] = mb;
      }
    }
  }

  return mb;
}


EdgeObject*
Device::get_edge_object(const Elem* elem, int edge)
{
  EdgeObject* mb = NULL;

  ID id = _bd_regions->get_edge_id(elem, edge);
  if (id != INVALID_ID)
  {
    EdgeObjMap::iterator it(_edge_map.find(id));
    if (it != _edge_map.end())
      mb = it->second;
  }

  return mb;
}



NodeObject*
Device::get_node_object(ID id)
{
  NodeObject* mb = NULL;

  if (id != INVALID_ID)
  {
    NodeObjMap::iterator it(_node_map.find(id));
    if (it != _node_map.end())
      mb = it->second;
    else
    {
      if (_bd_regions->is_node(id))
      {
        ModelOptions opts;
        mb = NodeObject::create(opts);
        _node_map[id] = mb;
      }
    }
  }

  return mb;
}


NodeObject*
Device::get_node_object(const Elem* elem, int node)
{
  NodeObject* mb = NULL;

  ID id = _bd_regions->get_node_id(elem->get_node(node));
  if (id != INVALID_ID)
  {
    NodeObjMap::iterator it(_node_map.find(id));
    if (it != _node_map.end())
      mb = it->second;
  }

  return mb;
}


NodeObject*
Device::get_node_object(const Node* node)
{
  NodeObject* mb = NULL;

  ID id = _bd_regions->get_node_id(node);
  if (id != INVALID_ID)
  {
    NodeObjMap::iterator it(_node_map.find(id));
    if (it != _node_map.end())
      mb = it->second;
  }

  return mb;
}



void
Device::get_active_region_ids(const string& name, set<ID>& ids) const
{
  ids.clear();

  if (name == "all")
    ids = get_active_region_ids();
  else
  {
    ClusterMap::const_iterator clit(_cluster_map.find(name));
    if (clit != _cluster_map.end())
    {
      for (size_t i = 0; i < clit->second.size(); ++i)
        ids.insert((clit->second)[i]);
    }
    else
    {
      map<ID, string>::const_iterator it(_region_names.begin());
      const map<ID, string>::const_iterator end(_region_names.end());

      for ( ; it != end; ++it)
        if (it->second == name)
          ids.insert(it->first);

      // next we look in the mesh region list
      if (ids.size() == 0)
      {
        const IDSet& idset = _mesh_region_info->get_ids(name);
        IDSet::iterator it(idset.begin());
        const IDSet::iterator end(idset.end());

        for ( ; it != end; ++it)
          if (_region_names.count(*it))
            ids.insert(*it);
      }

      // as last resort we look for material names
      if (ids.size() == 0)
      {
        MaterialMap::const_iterator it(_material_map.begin());
        const MaterialMap::const_iterator end(_material_map.end());

        for ( ; it != end; ++it)
          if (it->second->get_name() == name)
            ids.insert(it->first);
      }
    }
  }
}



void
Device::get_mesh_region_ids(const string& name, vector<ID>& ids) const
{
  const IDSet& idset = _mesh_region_info->get_ids(name);
  ids.clear();
  ids.reserve(idset.size());

  IDSet::iterator it(idset.begin());
  const IDSet::iterator end(idset.end());
  for ( ; it != end; ++it)
    ids.push_back(*it);

  ids.reserve(ids.size());
}




void
Device::get_boundary_region_ids(const string& name, vector<ID>& ids) const
{
  IDSet idset(_bd_regions->get_ids(name));

  // if the set is empty, we try to interpret the string as a specification
  // of the boundary between two materials
  if (idset.empty())
  {
    vector<string> comp;
    Utils::tokenize(name, comp, "/");

    // there have to be two _different_ components
    if ((comp.size() == 2) && (comp[0] != comp[1]))
    {
      // it is an interface specification

      // a map to keep track of existing IDs
      map<IDPair, ID> known_ids;

      // we have to loop over all element sides!
      // we only look on level zero
      MeshBase::const_element_iterator el(get_mesh().level_elements_begin(0));
      const MeshBase::const_element_iterator el_end(get_mesh().level_elements_end(0));

      for ( ; el != el_end; ++el)
      {
        Elem* elem = *el;
        const ID id = elem->subdomain_id();

        // the material of this element
        const Material* mat = get_material(id);

        // the name of this elements region
        const string& reg_name = get_region_name(id);

        // the other component
        size_t other;
        if ((comp[0] == mat->get_name()) || (comp[0] == reg_name))
          other = 1;
        else if ((comp[1] == mat->get_name()) || (comp[1] == reg_name))
          other = 0;
        else // go to the next element
          continue;


        // loop over the sides
        int n_sides = elem->n_sides();
        for (int s = 0; s < n_sides; s++)
        {

          // check if the neighbouring element is inexistent (outer boundary)
          // or in another simulation region (inner boundary)
          //
          // we allow inner 'boundaries', i.e. we don't really consider
          // boundaries but n-1 dimensional domains
          const Elem* neighbour = elem->neighbor(s);
          ID neighbour_id = INVALID_ID;

          // subdomain IDs have to be different
          if ((neighbour == NULL) ||
              ((neighbour_id = neighbour->subdomain_id()) != id))
          {
            map<IDPair, ID>::iterator it(known_ids.find(IDPair(id, neighbour_id)));
            // if the ID pair already exists, we can just add the elem side
            if (it != known_ids.end())
            {
              _bd_regions->add_side(elem, s, it->second);
            }
            else
            {

              // check if it has to be added
              bool add = false;
              if (neighbour_id == INVALID_ID)
                add = !comp[other].compare("-");
              else
                add = ((comp[other] == get_material(neighbour_id)->get_name()) ||
                    (comp[other] == get_region_name(neighbour_id)));

              if (add)
              {
                // first check if the side has already a boundary ID
                ID newid = _bd_regions->get_side_id(elem, s);
                if (newid == INVALID_ID)
                {
                  // we have to create a new one
                  newid = _bd_regions->next_id();
                  _bd_regions->add_side(elem, s, newid);
                }
                _bd_regions->set_name(newid, name);
                known_ids[IDPair(id, neighbour_id)] = newid;
                idset.insert(newid);
              }
            }
          }
        }
      }
    }
  }

  ids.clear();
  ids.reserve(idset.size());

  IDSet::iterator it(idset.begin());
  const IDSet::iterator end(idset.end());
  for ( ; it != end; ++it)
    ids.push_back(*it);
}



void
Device::set_region_name(const string& name, const vector<ID>& ids)
{
  for (unsigned int i = 0 ; i < ids.size(); ++i)
    _region_names[ids[i]] = name;
}


/*
void
Device::set_boundary_region_name(const string& name, const vector<ID>& ids)
{
  for (unsigned int i = 0 ; i < ids.size(); ++i)
    _boundary_region_names[ids[i]] = name;
}
*/


void
Device::set_cluster(const string& name, const vector<ID>& ids)
{
  if (_cluster_map.find(name) != _cluster_map.end()) {
    string msg("Cluster ");
    msg += name;
    msg += " already defined.";
    throw InitFailedException(msg);
  }

  _cluster_map[name] = ids;
}






