// $Id$

#include "Pyropolarization.h"
#include "Database.h"
#include "Material.h"

// The first string is the class name, the second one
// is the type of the model (here it is a polarization model),
// the third one is the specific model implementation.
// The library name will then be polarization_piezo.so
TIBER_MODULE(Pyropolarization, polarization, pyro)

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
  Database& db = get_database();
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
