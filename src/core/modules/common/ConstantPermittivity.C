// $Id$

#include "ConstantPermittivity.h"
#include "Database.h"
#include "Material.h"
#include "tensor.h"

// The first string is the class name, the second one
// is the type of the model (here it is a polarization model),
// the third one is the specific model implementation.
// The library name will then be permittivity_constant.so
TIBER_MODULE(ConstantPermittivity, permittivity, constant)

using namespace std;

void
ConstantPermittivity::do_init(void)
{
 
  get_parameter("permittivity", _permittivity_diag, true);
  set_permittivity(_permittivity_diag);
 
  rotate();
}



void
ConstantPermittivity::read_database(void)
{
   Database& db = get_database();
   db.set_section("permittivity");

   _permittivity_diag(0) = _permittivity_diag(1) = _permittivity_diag(2) = 1;
   db.get("permittivity", _permittivity_diag, true);
  
}
