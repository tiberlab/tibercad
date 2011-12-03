// $Id$

#include "LatticeMismatch.h"
#include "Material.h"
#include "RotatedCrystal.h"
#include "tensor_value.h"
// The first string is the class name, the second one
// is the type of the model (here it is a bulk model),
// the third one is the specific model implementation.
// The library name will then be bulk_default.so

TIBER_MODULE(LatticeMismatch, body_force, lattice_mismatch)

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


  RealGradient body_force(0);
  
 //Get reference lattice
  const std::string name = get_option("reference_material", "");
 
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


  const Material* mat = get_material();
  const RotatedCrystal* ref_crystal =  &(ref_mat->get_rotated_crystal());
  
  double ref_lat_const[3];
   
  ref_crystal->get_lat_const(ref_lat_const);
 
  //Get local lattice

  
  const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());
       
  Tensor2Sym _eps0 = crystal_el->get_eps0(ref_lat_const);

  //double loc_lat_const[3];
  //crystal_el->get_lat_const(loc_lat_const);

  //Rotation
  Tensor2Gen RotMatrix = crystal_el->RotMatrix;
  _eps0 = sym(RotMatrix * ( _eps0 * (RotMatrix.transpose())));


  RealTensor eps0(0);
  for (ID i = 0; i<3; i ++)
    for (ID j = i; j<=i; j ++)
    {
      eps0(i,j) = _eps0(i+1,j+1);
      eps0(j,i) =  eps0(i,j);
    }


  double relax = get_option("relaxation_factor", 1.0);

  eps0 *= relax;


  set_strain_source(eps0);

  //RealTensor dummy_tens(0);
  //set_stress_source(dummy_tens);

  //RealGradient dummy_grad(0);
  //set_force_source(dummy_grad);

  delete ref_mat;
 
}




 
