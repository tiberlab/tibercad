// $Id$

#include "Device.h"
#include "InputParser.h"
#include "Material.h"
#include "Alloy.h"
#include "MeshUtils.h"
#include "MeshReader.h"
#include "DataOutput.h"
#include "BoundaryRegions.h"
#include "MaterialBoundary.h"
#include "EdgeObject.h"
#include "NodeObject.h"
#include "MeshRegionInfo.h"
#include "SimulationOptions.h"
#include "AtomisticStructure.h"
#include "QuantumContact.h"
#include "ExternalProfile.h"

#include "libmesh/gmsh_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/elem.h"

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

Device::Device(const ModelOptions& options)
  : _mesh(NULL),
    _mesh_units(1e-6),
    _eq_system(NULL),
    _boundary_nodes(NULL),
    _mesh_region_info(NULL),
    _bd_regions(NULL),
    _options(options),
    _symmetry(TiberCad::NONE)
{
  _material_map.clear();

  libMesh::Parallel::Communicator& comm = TiberCad::get_mpi_comm();
  unsigned int nodes_per_device = comm.size();
  unsigned int nodes_per_mesh = comm.size();

  int color = 0;

  if (_options.has_submodel("Parallel"))
  {
    const ModelOptions& mpi_opts = _options.submodels_begin("Parallel")->second;

    nodes_per_device = mpi_opts.get_option("mpi_processes_per_device", nodes_per_device);
    nodes_per_mesh = mpi_opts.get_option("mpi_processes_per_mesh", nodes_per_device);

    if (nodes_per_mesh > nodes_per_device)
      throw InitFailedException("Cannot distribute mesh on more nodes than "
          "requested for device.");
  }

  if (nodes_per_device < comm.size())
  {
    int proc_id = comm.rank();
    color = proc_id / nodes_per_device;

    //MPI_Comm local_comm;
    //MPI_Comm_split(comm.get(), color, 0, &local_comm);
    //_mpi_comm = local_comm;

    // this did not work before 1.0.0, because there was a bug in libmesh
    // parallel_implementation.h, line 470, missing this->assign(comm)
    comm.split(color, 0, _mpi_comm);
  }
  else
  {
    // TODO not sure if we should duplicate here?
    //_mpi_comm = TiberCad::get_mpi_comm();
    _mpi_comm.duplicate(TiberCad::get_mpi_comm());

    if (nodes_per_device > comm.size())
      throw InitFailedException("Too many MPI nodes requested for device");
  }


  if (nodes_per_mesh < nodes_per_device)
  {
    int proc_id = _mpi_comm.rank();
    int mesh_color = proc_id / nodes_per_mesh;

    _mpi_comm.split(mesh_color, 0, _mesh_comm);
  }
  else
  {
    _mesh_comm.duplicate(_mpi_comm);
  }

  ostringstream os;
  os << color;
  InputParser::add_defined("MPI_DEV_KEY", os.str(), false);
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

  _mesh_comm.clear();
  _mpi_comm.clear();

  delete _eq_system;
  delete _boundary_nodes;
  delete _mesh;
  delete _mesh_region_info;
  delete _bd_regions;
}


Device*
Device::create(const ModelOptions& options)
{

  Device* device = new Device(options);

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
    libMesh::UniquePtr<DataOutput> writer(DataOutput::create("vtk"));
    if (writer.get() != NULL)
    {
      writer->set_output_directory(_options.get_option("output_path", "./"));
      writer->set_filename(Utils::basename(_options["meshfile"]) + "_bnd");

      libMesh::UniquePtr<MeshBase> bdmesh = MeshUtils::create_boundary_mesh(get_mesh());

      writer->set_mesh(*bdmesh);
      writer->write(true);
    }
  }

  if (_options.get_option("plot_alloy_composition", false))
  {
    if (get_mesh().comm().size() > 1)
      Messages::warning("plot_alloy_composition not implemented in parallel.");
    else
    {
    string format = "vtk";
    if (get_mesh().mesh_dimension() == 1)
      format = "dat";

    DataOutput writer(get_mesh(), format);
    writer.set_output_directory(_options.get_option("output_path", "./"));
    //writer->set_filename(Utils::basename(_options["meshfile"]) + "_alloy_comp");
    //writer.set_mesh(get_mesh());

    vector<double> data(get_mesh().n_active_elem(), 0.0);
    size_t ctr = 0;
    MeshBase::element_iterator el(get_mesh().active_elements_begin());
    MeshBase::element_iterator end(get_mesh().active_elements_end());
    for ( ; el != end; ++el, ++ctr)
    {
      const Elem* elem = *el;
      const Material* mat = get_material(elem);
      if (mat->is_alloy())
        data[ctr] = static_cast<const Alloy*>(mat)->get_molar_fraction();
    }

    vector<string> legend(1, "x");
    string filename = Utils::basename(_options["meshfile"]) + "_alloy_comp";
    writer.write_cell_data(filename, data, legend);
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

  if (meshfile.empty())
    throw InitFailedException("'meshfile' must be provided in Device block.");

  {
    ostringstream os;
    os << "Reading mesh file \'" << meshfile << " ...";
    m.info(os.str());
    m.indent();
  }

  _mesh = new libMesh::Mesh(_mesh_comm, dim);

  _mesh_region_info = new MeshRegionInfo(*_mesh);
  _bd_regions = new BoundaryRegions(*_mesh);

  MeshReader::read_mesh(meshfile, *_mesh, *_mesh_region_info, *_bd_regions);

  // TODO only for now (back compatibility)
  delete _boundary_nodes;
  _boundary_nodes = new map<ID, vector<unsigned int> >();
  _bd_regions->get_bc_node_map(*_boundary_nodes);


  // update mesh dimension
  dim = _mesh->mesh_dimension();

  //m.info("done.");
  m.newline();
  {
    ostringstream os;

    os << "mesh units          : " << _mesh_units << " m" << endl
       << "mesh dimension      : " << setw(7) << setfill(' ') << dim << endl
       << "number of nodes     : " << setw(7) << setfill(' ') << _mesh->n_nodes() << endl
       << "number of elements  : " << setw(7) << setfill(' ') << _mesh->n_elem() << endl
       << "number of subdomains: " << setw(7) << setfill(' ') 
       << _mesh_region_info->n_subdomains() << endl
       << "number of boundaries: " << setw(7) << setfill(' ') << _bd_regions->n_subdomains();
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


  _eq_system = new libMesh::EquationSystems(*_mesh);
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

  //
  // 2014-04-03 A way to treat variable alloy composition:
  //   subdivide the region in subregions according to intervals
  //   in the composition.
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

    // No default material
    string material = _options.get_option("material", "none");
    string x_frac = _options.get_option("x", "0.0");
    string xdir = _options.get_option("x-growth-direction", "");
    string ydir = _options.get_option("y-growth-direction", "");
    string zdir = _options.get_option("z-growth-direction", "");

    material = data.get_option("material", material);
    data["material"] = material;
    x_frac = data.get_option("x", x_frac);

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

    if (!data.has_submodel("Doping"))
    {
      ModelOptions::const_submodel_iterator dop_it(
          _options.submodels_begin("Doping"));
      if (dop_it != _options.submodels_end("Doping"))
        data.add_submodel("Doping", dop_it->second);
    }

    // create the materials
    data["x"] = x_frac;
    Material* mat = Material::create(material, data);

    set_material(mat, region_ids, data.get_name());

    // now we check if we should use variable alloy composition
    if (data.has_submodel("alloy_composition"))
    {
      const ModelOptions& alloy =
          data.submodels_begin("alloy_composition")->second;
      Messages::info("decompose region " + data.get_name() +
          " to account for variable alloy composition.");

      vector<double> comp;

      // this will change the region_ids, but it needs the original ones
      decompose_region(alloy, region_ids, comp);

      for (int i = 0; i < comp.size(); ++i)
      {
        //cerr << region_ids[i] << " -> " << comp[i] << endl;

        data.set_option("x", comp[i]);
        Material* mat = Material::create(material, data);

        ostringstream os;
        os << data.get_name() << "_" << i;
        set_material(mat, vector<ID>(1, region_ids[i]), os.str());
      }

      set_cluster(data.get_name(), region_ids);

    }


  }

  m.unindent();
  m.info("Creation of materials done.");

  Messages::debug("Control::create_materials() end");
}


const std::string&
Device::get_region_name(ID id) const
{
  std::map<ID, std::string>::const_iterator it(_region_names.find(id));

  if (it == _region_names.end())
  {
    std::ostringstream s;
    s << "Tried to access unknown region with id " << id;
    throw (DeviceException(s.str()));
  }
  return it->second;
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
    ModelOptions data(it->second);
    // use device global material as reference  
    string material = _options.get_option("material", "none");  
    string x_frac = _options.get_option("x", "0.0");
    string xdir = _options.get_option("x-growth-direction", "");
    string ydir = _options.get_option("y-growth-direction", "");
    string zdir = _options.get_option("z-growth-direction", "");
    material = data.get_option("reference_material", material);

    if (material != "")
    {
      ModelOptions refopts;
      refopts.set_name(material);
      refopts["x"] = x_frac;
      refopts.set_option("dimension", get_mesh().mesh_dimension());
 
      xdir = refopts.get_option("x-growth-direction", xdir);
      if (!xdir.empty())
        refopts["x-growth-direction"] = xdir;
      ydir = refopts.get_option("y-growth-direction", ydir);
      if (!ydir.empty())
        refopts["y-growth-direction"] = ydir;
      zdir = refopts.get_option("z-growth-direction", zdir);
      if (!zdir.empty())
        refopts["z-growth-direction"] = zdir;
     
      ModelOptions::submodel_iterator it(data.submodels_begin("reference_material"));
      if (it != data.submodels_end("reference_material"))
      {
        refopts += it->second;
        data.delete_submodels("reference_material");
      }
      data.add_submodel("reference_material", refopts);
    }

    const string& st_name = data.get_name();

    AtomisticStructure* st = AtomisticStructure::create();

    // Defined atomistic structure is put in the atomistic_structure_map
    _atomistic_structure_map[st_name] = st;

    //WARNING: For debugging purposes, initialization of
    //atomistic structures is here, but it's not the right place! (maybe it is...)
    st->init(st_name, this, data);

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
    Messages::info("Creating quantum contact " + name);

    std::set<ID> rg_ids, bd_ids;

    // We get mesh regions (not active regions).

    extract_physical_regions(regions, rg_ids);

    if (rg_ids.size()==0)
      throw InitFailedException("In Quantum_contact user must define a valid mesh region");

    std::vector<ID> bd_ids_v;
    get_boundary_region_ids(name, bd_ids_v);
    bd_ids.insert(bd_ids_v.begin(), bd_ids_v.end());

    ID  newid = _mesh_region_info->next_id();

    QuantumContact* st = QuantumContact::create();

    st->init(newid, name, this, _bd_regions, rg_ids, bd_ids, length);

    // Defined quantum contact is put in the quantum_contact_map
    _quantum_contact_map[name] = st;

    _mesh_region_info->add_id(newid);
    _mesh_region_info->set_name(newid, name);

    std::vector<ID> vid(1, newid);

    // Only one material/region can be created (for now) !!

    const Material* refmat = get_material(*rg_ids.begin());
    const ModelOptions& refopts = refmat->get_options();
    // here we can be sure that material is given
    string material = refopts.get_option("material", "");

    Material* newmat = Material::create(material, refopts);
    set_material(newmat, vid, name);

    // We have to erase the region id from the list of active regions, otherwise
    // we mess up the other modules (quantum contact regions should be invisible
    // when asking for real device regions)
    _active_region_ids.erase(newid);

  }

  if (_options.has_submodel("Quantum_contact") )
  {

    string meshfile = _options["meshfile"];
    meshfile = "new_" + meshfile;

    //ostringstream os;
    //os << "Writing mesh with quantum contacts: " << meshfile;
    //Messages::info(os.str());

    //libMesh::GmshIO(*_mesh).write(meshfile);

   }

  Messages::debug("Control::create_quantum_contacts() end");
}


void
Device::extract_physical_regions(const std::string& str, IDSet& ids) const
{
  ids.clear();

  // the IDs that have to be excluded ("-pippo" syntax)
  IDSet exclude;
  // the IDs that have to be included
  IDSet include;

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
      include.insert(preg_ids.begin(), preg_ids.end());
  }

  if (!exclude.empty())
  {
    ids = get_active_region_ids();

    for (IDSet::iterator it(exclude.begin()); it != exclude.end(); ++it)
      ids.erase(*it);
  }

  ids.insert(include.begin(), include.end());
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
  m.newline();
  m.info("Setting up bulk models ...");
  m.indent();
  MaterialMap::iterator it(_material_map.begin());
  const MaterialMap::iterator end(_material_map.end());
  for ( ; it != end; ++it) {
    (it->second)->init();
  }

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


	libMesh::EquationSystems&
Device::get_equation_systems(MeshBase* mesh)
{
	libMesh::EquationSystems* eqsys;
  if (mesh == _mesh)
    eqsys = _eq_system;
  else
  {
    map<const MeshBase*, libMesh::EquationSystems*>::iterator it =
        _eq_sys_map.find(mesh);

    if (it == _eq_sys_map.end())
    {
      eqsys = new libMesh::EquationSystems(*mesh);
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
    if (libMesh::global_n_processors() == 1)
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
        //assert(ids.size() == 2);
        ID idA = *(ids.begin());
        ID idB = idA;
        if (ids.size() == 2)
          idB = *(++(ids.begin()));

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
  IDSet idset;
  get_boundary_region_ids(name, idset);

  ids.resize(idset.size());
  std::copy(idset.begin(), idset.end(), ids.begin());
}

void
Device::get_boundary_region_ids(const string& name, IDSet& ids) const
{
  ids = _bd_regions->get_ids(name);

  // if the set is empty, we try to interpret the string as a specification
  // of the boundary between two materials
  if (ids.empty())
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

        // if id is not used, skip to next element
        if (!get_active_region_ids().count(id))
          continue;

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

          // subdomain IDs have to be different, and be part of the simulation
          if ((neighbour == nullptr) ||
              ((neighbour_id = neighbour->subdomain_id()) != id) &&
               get_active_region_ids().count(neighbour_id))
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
              {
                add = ((comp[other] == get_material(neighbour_id)->get_name()) ||
                    (comp[other] == get_region_name(neighbour_id)));
              }

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
                ids.insert(newid);
              }
            }
          }
        }
      }
    }
  }
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



void
Device::reassign_alloy_regions(const string& source,
    const vector<ID>& region_ids,
    const vector<double>& composition)
{
  Messages m;
  m.info("Re-assigning alloy subregion IDs...");
  m.indent();

  pair<SimulationInterface*, ID> provider(
      SimulationInterface::find_solution_provider(source));

  set<ID> reg_ids;
  for (int i = 0; i < region_ids.size(); ++i)
    reg_ids.insert(region_ids[i]);

  if (provider.first != NULL)
  {
    m.info("Alloy composition provided by " + provider.first->get_name());

    MeshBase::const_element_iterator el(get_mesh().level_elements_begin(0));
    const MeshBase::const_element_iterator el_end(get_mesh().level_elements_end(0));

    for ( ; el != el_end; ++el)
    {
      Elem* elem = *el;
      ID id = elem->subdomain_id();

      if (reg_ids.count(id))
      {
        double x;
        provider.first->get_solution(
            elem, provider.second, x, elem->centroid());

        int i = 0;
        while ((i < (composition.size() - 1)) &&
               (x > 0.5*(composition[i] + composition[i+1])))
          i++;

        elem->subdomain_id() = region_ids[i];
      }
    }
  }
  else
  {
    // it might be an Atomistic:tb like string
    vector<string> tokens;
    Utils::tokenize(source, tokens, ":");
    if ((tokens[0] != "Atomistic") || (tokens.size() < 2))
      throw InitFailedException("You must provide a module name or atomistic "
          "structure for external alloy profile source");

    const AtomisticStructure* str = get_atomistic_structure(tokens[1]);

    if (str == NULL)
      throw InitFailedException("Unknown atomistic structure: " + tokens[1]);

    m.info("Alloy composition provided by atomistic structure " + tokens[1]);
    m.indent();

    double cutoff = str->get_options().get_option("control_volume_radius", 1.0);
    {
      ostringstream os;
      os << "control volume radius is " << cutoff << " nm";
      m.info(os.str());
    }

    //map<Specie, vector<unsigned int>> stats;
    //str->extract_statistics(stats, reg_ids, cutoff);

    double scale = 1e-9 / this->get_mesh_units();
    double control_vol = cutoff * scale;
    switch (this->get_mesh().mesh_dimension())
    {
      case 3:
        control_vol *= 4.0/3.0 * cutoff * scale;
        [[fallthrough]];
      case 2:
        control_vol *= M_PI * cutoff * scale;
        break;
    }

    // Idea: get the local concentration at the nodes for a given control volume,
    //       then take the mean value for the element.This should work well on small
    //       elements and worse on big elements.

    MeshBase::element_iterator el(get_mesh().local_elements_begin());
    MeshBase::element_iterator end(get_mesh().local_elements_end());

    for ( ; el != end; ++el)
    {
      Elem* elem = *el;
      ID id = elem->subdomain_id();

      if (reg_ids.count(id))
      {
        Material* mat = get_material(elem);

        if (mat->is_alloy())
        {
          const Alloy* alloy = static_cast<const Alloy*>(mat);
          const Material* matA = alloy->get_component_A();
          const Material* matB = alloy->get_component_B();

          // this is the local composition
          double x = alloy->get_molar_fraction();

          //
          // if the element is small
          if (elem->volume() < control_vol)
          {
            const vector<unsigned int>& atoms = str->get_atoms_in_elem(elem);
            int atom = -1;

            //
            // it may not even contain atoms
            if (!atoms.empty())
            {
              // ok, there is at least one atom inside
              // atomic coordinates are in Angstrom
              Point center = 10 * elem->centroid();
              double min_dist = Point(center -
                  str->get_structure_atom(atoms[0]).get_position()).size();
              unsigned int nearest = 0;

              // look for the atom nearest to the center
              for (unsigned int i = 1; i < atoms.size(); ++i)
              {
                double dist = Point(center -
                    str->get_structure_atom(atoms[i]).get_position()).size();
                if (dist < min_dist)
                {
                  min_dist = dist;
                  nearest = i;
                }
              }
              atom = atoms[nearest];
            }
            else
            {
              // we try to find some nearby atom
              set<const Elem*> processed_elems;
              set<const Elem*> to_process;
              to_process.insert(elem);

              unsigned int nearest = 0;
              double min_dist = 100 * 10 * cutoff * scale;

              while (!to_process.empty())
              {
                set<const Elem*>::iterator it(to_process.begin());
                const Elem* next_el = *it;

                const vector<unsigned int>& atoms = str->get_atoms_in_elem(next_el);
                if (!atoms.empty())
                {
                  // ok, there is at least one atom inside
                  // atomic coordinates are in Angstrom
                  Point center = 10 * elem->centroid();
                  double min_dist = Point(center -
                      str->get_structure_atom(atoms[0]).get_position()).size();
                  unsigned int nearest = 0;

                  // look for the atom nearest to the center
                  for (unsigned int i = 1; i < atoms.size(); ++i)
                  {
                    double dist = Point(center -
                        str->get_structure_atom(atoms[i]).get_position()).size();
                    if (dist < min_dist)
                    {
                      min_dist = dist;
                      nearest = i;
                    }
                  }
                  atom = atoms[nearest];
                }

                for (int s = 0; s < next_el->n_sides(); s++)
                {

                  const Elem* neigh = next_el->neighbor(s);

                  if ((neigh != NULL) &&
                      reg_ids.count(neigh->subdomain_id()) &&
                      !processed_elems.count(neigh) &&
                      (Point(elem->centroid() - neigh->centroid()).size() <
                          scale * 3 * cutoff))
                  {
                    to_process.insert(neigh);
                  }
                }

                processed_elems.insert(next_el);
                to_process.erase(it);
              }
            }

            //cerr << atom << " " << str->get_structure_atom(atom).get_specie() << " : ";

            map<Specie, unsigned int> counts;
            str->extract_statistics(atom, counts, reg_ids, cutoff);


            double sum = 0;
            map<Specie, unsigned int>::iterator stat_it(counts.begin());
            map<Specie, unsigned int>::iterator stat_end(counts.end());
            for ( ; stat_it != stat_end; ++stat_it)
            {
              bool inA = matA->has_specie(stat_it->first);
              bool inB = matB->has_specie(stat_it->first);
              if (inA && !inB) 
              {
                x = stat_it->second;
                sum += x;
              }
              else if (!inA && inB)
                sum += stat_it->second;
            }
            //cerr << x << " " << sum << endl;
            x /= sum;
          }
          else
          {
            // there must be at least an atom inside, otherwise there is something strange

            const vector<unsigned int>& atoms = str->get_atoms_in_elem(elem);
            int atom = -1;


            x = 0;
            for (unsigned int i = 0; i < atoms.size(); ++i)
            {
              map<Specie, unsigned int> counts;
              str->extract_statistics(atoms[i], counts, reg_ids, cutoff);

              double conc = 0;
              double sum = 0;
              map<Specie, unsigned int>::iterator stat_it(counts.begin());
              map<Specie, unsigned int>::iterator stat_end(counts.end());
              for ( ; stat_it != stat_end; ++stat_it)
              {
                bool inA = matA->has_specie(stat_it->first);
                bool inB = matB->has_specie(stat_it->first);
                if (inA && !inB)
                {
                  conc = stat_it->second;
                  sum += conc;
                }
                else if (!inA && inB)
                  sum += stat_it->second;
              }

              x += conc / sum;
            }
            x /= atoms.size();

          }


          int i = 0;
          while ((i < (composition.size() - 1)) &&
              (x > 0.5*(composition[i] + composition[i+1])))
            i++;

          elem->subdomain_id() = region_ids[i];
        }
      }
    }
  }
}


void
Device::decompose_region(const ModelOptions& options,
    vector<ID>& region_ids, vector<double>& composition)
{

  ExternalProfile* profile = nullptr;
  if (!options.get_name().empty() && (options.get_name() != "extern"))
    profile = ExternalProfile::create(options);

  double min_x = 0.0;
  double max_x = 1.0;

  if (profile != nullptr)
  {
    min_x = profile->get_min_max().first;
    max_x = profile->get_min_max().second;
  }

  int intervals = options.get_option("intervals", 10);
  double delta_rel = 1.0 / intervals;
  min_x = options.get_option("min_content", min_x);
  max_x = options.get_option("max_content", max_x);
  double delta = (max_x - min_x) * delta_rel;

  int n_orig = region_ids.size();
  double mean_x = 0.5 * (max_x + min_x);

  vector<ID> reg_ids;
  reg_ids.reserve(intervals + 1);
  composition.reserve(intervals + 1);

  ID id = _mesh_region_info->next_id();

  for (double x = min_x; x <= max_x; x += delta, ++id)
  {
    _mesh_region_info->add_id(id);
    reg_ids.push_back(id);
    composition.push_back(x);
  }

  region_ids.insert(region_ids.begin(), reg_ids.begin(), reg_ids.end());

  if (profile != nullptr)
  {

    set<ID> used_ids;
    for (int i = 0; i < region_ids.size(); ++i)
      used_ids.insert(region_ids[i]);

    //srand(time(NULL));

    MeshBase::const_element_iterator el(get_mesh().level_elements_begin(0));
    const MeshBase::const_element_iterator el_end(get_mesh().level_elements_end(0));

    for ( ; el != el_end; ++el)
    {
      Elem* elem = *el;
      ID id = elem->subdomain_id();

      if (!used_ids.count(id))
        continue;

      double x = profile->get_data(elem);
      //cerr << "x = " << x << endl;

      int i = 0;
      while ((i < (composition.size() - 1)) &&
          (x > 0.5*(composition[i] + composition[i+1])))
        i++;

      elem->subdomain_id() = region_ids[i];

      //double rnd = static_cast<double>(rand()) / RAND_MAX;
      //int offset = floor(rnd / delta_rel);

      //id = first_id + offset;
      //elem->subdomain_id() = id;
      //used_ids.insert(offset);
    }


    /*
    set<ID>::iterator it(used_ids.begin());
    const set<ID>::iterator end(used_ids.end());

    for ( ; it != end; ++it)
    {
      ID id = *it;
      _mesh_region_info->add_id(first_id + id);
      region_ids.push_back(first_id + id);
      composition.push_back(id*delta + min_x);
    }
    */
  }
  else
  {
    string source = options.get_option("source", "");

    // it might be an Atomistic:tb like string
    vector<string> tokens;
    Utils::tokenize(source, tokens, ":");

    if ((tokens[0] == "Atomistic") && (tokens.size() == 2))
    {
      AtomisticStructure::register_callback(tokens[1],
          boost::bind(&Device::reassign_alloy_regions, this, source,
              region_ids, composition));
    }
    else
      SimulationInterface::register_callback(source,
          boost::bind(&Device::reassign_alloy_regions, this, source,
              region_ids, composition));
  }

}

