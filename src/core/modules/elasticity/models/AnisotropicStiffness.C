// $Id$

#include "AnisotropicStiffness.h"
#include "Database.h"

// The first string is the class name, the second one
// is the type of the model (here it is a bulk model),
// the third one is the specific model implementation.
// The library name will then be bulk_default.so

TIBER_MODULE(AnisotropicStiffness, Stiffness, anisotropic)

using namespace std;


AnisotropicStiffness::AnisotropicStiffness(const ModelOptions& options):StiffnessModel(options)
{
}


void AnisotropicStiffness::read_database( )
{

  Database& db = get_database();
  db.set_section("stiffness");

  _c11 = db.get("C11", 0.0,false);
  _c12 = db.get("C12", 0.0, false);
  _c13 = db.get("C13", 0.0, false);
  _c33 = db.get("C33", 0.0, false);
  _c44 = db.get("C44", 0.0, false);

}


void
AnisotropicStiffness::do_init(void)
{



  _c11 = get_option("C11", _c11);
  _c12 = get_option("C12", _c12);
  _c13 = get_option("C13", _c13);
  _c33 = get_option("C33", _c33);
  _c44 = get_option("C44", _c44);


  Tensor4DSym stiffness(0);
  if (get_material()->get_structure() == "wz")
  {
    
    stiffness(1,1,1,1) = _c11;
    stiffness(2,2,2,2) = _c11;
    stiffness(3,3,3,3) = _c33;
    stiffness(2,2,1,1) = _c12;
    stiffness(3,3,1,1) = _c13;
    stiffness(3,3,2,2) = _c13;
    stiffness(3,2,3,2) = _c44;
    stiffness(3,1,3,1) = _c44;
    stiffness(2,1,2,1) = 0.5*(_c11-_c12);
   
  }
  else
  {

    stiffness(1,1,1,1) = _c11;
    stiffness(2,2,2,2) = _c11;
    stiffness(3,3,3,3) = _c11;
    stiffness(2,2,1,1) = _c12;
    stiffness(3,3,1,1) = _c12;
    stiffness(3,3,2,2) = _c12;
    stiffness(3,2,3,2) = _c44;
    stiffness(3,1,3,1) = _c44;
    stiffness(2,1,2,1) = _c44;
    
  }
   
  set_stiffness_constant(stiffness);

  rotate();

}




 
