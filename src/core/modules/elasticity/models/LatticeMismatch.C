// $Id$

#include "LatticeMismatch.h"
#include "BulkCrystal.h"
#include "Material.h"
#include "TiberModule.h"
#include "Messages.h"

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
  const std::string name = get_option("reference_material", "");

  // TODO: we should be able to take a material from a region, too,
  //       and to use the crystal directions from the region.
  //       Also, we need to implement reading a subblock, here.
 
  // maybe the reference material should be defined in a subblock
  // but this should be done after a change in the core
  // dummy reads
  get_option("x-growth-direction", "");
  get_option("y-growth-direction", "");
  get_option("z-growth-direction", "");
  get_option("x", "");

  //  get_options().print_all();
  if (name.empty())
    throw InitFailedException("Lattice mismatch: reference material is not defined");

  ModelOptions matopts(get_options());

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

  //RealTensor dummy_tens(0);
  //set_stress_source(dummy_tens);

  //RealGradient dummy_grad(0);
  //set_force_source(dummy_grad);

  delete ref_mat;
 
}




 
