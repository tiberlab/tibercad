// $Id$

#include "Pyropolarization.h"
#include "Database.h"
#include "Material.h"

#include "TiberModule.h"



using namespace std;

void
Pyropolarization::do_init(void)
{
 
  // read the Pz componenents in the crystal system
  get_parameter("Pz", _Pz, true, initializer(&Pyropolarization::_initP));
  _initP();
  
  if (has_parameter("P"))
  {
    // read the polarization vector in the calc system
    get_parameter("P", polarization());
  }

}


void
Pyropolarization::read_database(void)
{

  // Read the Pz componenents in the crystal system
  const Database& db = get_database();
  db.set_section("pyroelectricity");
  _Pz = db.get("Pz", _Pz);

}


void
Pyropolarization::_initP(void)
{
  polarization() = 0;
  polarization()(2) = _Pz;
  rotate();
}
