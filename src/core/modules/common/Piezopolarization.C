// $Id$

#include "Piezopolarization.h"
#include "SimulationInterface.h"
#include "Database.h"
#include "Material.h"
#include "TensorOperators.h"


// The first string is the class name, the second one
// is the type of the model (here it is a polarization model),
// the third one is the specific model implementation.
// The library name will then be polarization_piezo.so
TIBER_MODULE(Piezopolarization, polarization, piezo)

using namespace std;


void
Piezopolarization::do_init(void)
{

  std::string sim_name = "";
  get_parameter("simulation_name", sim_name);
  _strain.set_simulation(sim_name);

}


void
Piezopolarization::read_database(void)
{
  // get piezoelectric_coefficients
  Database& db = get_database();
  db.set_section("piezoelectricity");
      
  if (get_material()->get_structure() == "wz")
  {
    _e33 = db.get("e33", 0.0, true);
    _e31 = db.get("e31", 0.0, true);
    _e15 = db.get("e15", 0.0, true);
  } 
  else if (get_material()->get_structure() == "zb")
    db.get("e14",0.0, true);

}


void
Piezopolarization::calculate(const Elem* elem, const Point& point)
{  

  RealVectorValue polarization(0);

  Tensor2Sym& strain = get_strain();
  _strain.get_crystal_strain(elem, point, strain);
  

  // compute polarization
  if (get_material()->get_structure() == "wz")
  {
    polarization(0) = 2.0 * _e15 * strain(3,1);
    polarization(1) = 2.0 * _e15 * strain(3,2);
    polarization(2) = _e31 * strain(1,1) + _e31 * strain(2,2) + _e33 * strain(3,3);
  }
  else
  {
    polarization(0) = 2.0 * _e14 * strain(3,2);
    polarization(1) = 2.0 * _e14 * strain(3,1);
    polarization(2) = 2.0 * _e14 * strain(2,1);
  }
  set_polarization(polarization);

  // rotate polarization to calculation system
  rotate();

}
