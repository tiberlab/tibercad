// $Id$

#include "LatticeMismatch.h"
#include "Material.h"
#include "RotatedCrystal.h"

// The first string is the class name, the second one
// is the type of the model (here it is a bulk model),
// the third one is the specific model implementation.
// The library name will then be bulk_default.so

TIBER_MODULE(LatticeMismatch, body_force, lattice_mismatch)

using namespace std;


LatticeMismatch::LatticeMismatch(const ModelOptions& options):BodyForceModel(options)
{
}




void
LatticeMismatch::do_init(void)
{


  RealGradient body_force(0);
  
 //Get reference lattice
  const std::string name = get_option("reference_material", "none"); 
 

  //  get_options().print_all();
  if (name == "none") throw InitFailedException("Lattice mismatch: reference material is not defined");

  //  cout<<get_material()->get_name()<<endl;

  Material* ref_mat = Material::create (name,get_options());
  ref_mat->init();

  const RotatedCrystal* ref_crystal =  &(ref_mat->get_rotated_crystal());
  
  double ref_lat_const[3];
   
  ref_crystal->get_lat_const(ref_lat_const);
 
  //Get local lattice

  const Material* mat = get_material();
  
  const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());
       
  Tensor2Sym _eps0 = crystal_el->get_eps0(ref_lat_const);


  cout<<_eps0(1,1)<<endl;
  cout<<_eps0(2,2)<<endl;
  cout<<_eps0(3,3)<<endl;
  cout<<endl;

  //double loc_lat_const[3];
  //crystal_el->get_lat_const(loc_lat_const);

  RealTensor eps0(0);
  for (ID i = 0; i<3; i ++)
    for (ID j = i; j<=i; j ++)
    {
      eps0(i,j) = _eps0(i+1,j+1);
      eps0(j,i) =  eps0(i,j);
    }


  double relax = get_options().get_option("relaxation_factor", 1.0);

  eps0 *=relax;
  //eps0 *= -1.0;

  cout<<eps0<<endl;
  set_strain_source(eps0);

  RealTensor dummy_tens(0);
  set_stress_source(dummy_tens);

  RealGradient dummy_grad(0);
  set_force_source(dummy_grad);

 
}




 
