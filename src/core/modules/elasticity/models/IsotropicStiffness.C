// $Id$

#include "IsotropicStiffness.h"
#include "Material.h"
#include "Database.h"

// The first string is the class name, the second one
// is the type of the model (here it is a bulk model),
// the third one is the specific model implementation.
// The library name will then be bulk_default.so

TIBER_MODULE(IsotropicStiffness, stiffness, isotropic)

using namespace std;


IsotropicStiffness::IsotropicStiffness(const ModelOptions& options):StiffnessModel(options)
{
}


// void IsotropicStiffness::read_database( )
// {

//   Database& db = get_database();
//   db.set_section("stiffness/isotropic");

//   _E = db.get("young", 0.0, false);
//   _P = db.get("poisson", 0.0, false);
 

// }


void
IsotropicStiffness::do_init(void)
{

  _E = 0.0;
  _P = 0.0;

  get_parameter("young",_E);
  get_parameter("poisson",_P);


  Tensor4DSym stiffness(0);

  double A = _E / (1 + _P) / (1 - 2 * _P); 


  stiffness(1,1,1,1) = 1 - _P;
  stiffness(2,2,2,2) = 1 - _P;
  stiffness(3,3,3,3) = 1 - _P;
  stiffness(2,2,1,1) = _P;
  stiffness(3,3,1,1) = _P;
  stiffness(3,3,2,2) = _P;

  stiffness(3,2,3,2) = (1 - 2.0 * _P)/2.0;
  stiffness(3,1,3,1) = (1 - 2.0 * _P)/2.0;
  stiffness(2,1,2,1) = (1 - 2.0 * _P)/2.0;

  stiffness *= A;

  set_stiffness_constant(stiffness);

  //rotate();
  
}




 
