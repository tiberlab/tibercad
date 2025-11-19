// $Id$

#include "Pyropolarization.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"



using namespace std;

void
Pyropolarization::do_init(void)
{
 
  // read the Pz componenents in the crystal system
  get_parameter("Pz", _Pz, true, initializer(&Pyropolarization::_initP));
  _initP();
  
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
  PolarizationModel::do_init();

  set_polarization(libMesh::RealVectorValue(0.0, 0.0, _Pz));
  rotate();
}
