// $Id$

#include "LatticeMismatch.h"
#include "tibercad/atomistic/BulkCrystal.h"
#include "tibercad/physics/Material.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/io/Messages.h"

#include "tibercad/module/TiberModule.h"

#include "libmesh/tensor_value.h"


using namespace std;


LatticeMismatch::LatticeMismatch(const ModelOptions& options) :
    BodyForceModel(options)
{
}

LatticeMismatch::~LatticeMismatch(void)
{

}


void
LatticeMismatch::do_init(void)
{


  libMesh::RealGradient body_force(0);
  
  //Get reference lattice
  std::string name;
  ModelOptions matopts;

  if (get_options().has_submodel("reference_material"))
  {
    auto ref = get_options().submodels_begin("reference_material");
    matopts = ref->second;

    if (!matopts.find_option("material"))
    {
      matopts["material"] = matopts.get_name();
    }
    
    name = matopts["material"];
  }
  else
  {
    name = get_option("reference_material", "");

    SimulationInterface *sim = SimulationInterface::get_simulation(get_simulator_id());
    const Device& dev = sim->get_environment().get_device();
    const Material* mat = dev.get_material(name);

    if (mat != nullptr)
    {
      matopts = mat->get_options();
      name = mat->get_name();
    }
    else if (!name.empty())
    {
      matopts = get_options();

      get_option("x-growth-direction", "");
      get_option("y-growth-direction", "");
      get_option("z-growth-direction", "");
      get_option("euler_angles", "");
      get_option("x", "");
    }
  }

  if (name.empty())
    throw InitFailedException("Lattice mismatch: reference material is not defined");


  // for "official" materials (created from Device) this is always present
  unsigned int dim = get_material()->get_options().get_option("dimension", 3);
  matopts.set_option("dimension", dim);

  Material* ref_mat = Material::create(name, matopts);
  ref_mat->init();
  const BulkCrystal* ref_bulk = ref_mat->get_bulk_crystal();

  const Material* mat = get_material();
  const BulkCrystal* mat_bulk = mat->get_bulk_crystal();

  libMesh::RealTensor eps0(0);

  if ((mat_bulk != nullptr) && (ref_bulk != nullptr))
    mat_bulk->get_lattice_matching_strain(*ref_bulk, eps0);
  
  double relax = get_option("relaxation_factor", 1.0);

  eps0 *= relax;


  set_strain_source(eps0);

  delete ref_mat;
 
}




 
