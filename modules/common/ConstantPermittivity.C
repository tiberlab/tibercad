// $Id$

#include "ConstantPermittivity.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"


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
   const Database& db = get_database();
   db.set_section("permittivity");

   _permittivity_diag(0) = _permittivity_diag(1) = _permittivity_diag(2) = 1;
   db.get("permittivity", _permittivity_diag, true);
  
}
