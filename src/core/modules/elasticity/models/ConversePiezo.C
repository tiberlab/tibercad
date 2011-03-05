// $Id$

#include "ConversePiezo.h"
#include "Material.h"
#include "RotatedCrystal.h"
#include "Database.h"
#include "TensorOperators.h"

// The first string is the class name, the second one
// is the type of the model (here it is a bulk model),
// the third one is the specific model implementation.
// The library name will then be bulk_default.so

TIBER_MODULE(ConversePiezo,body_force, converse_piezo)

using namespace std;


ConversePiezo::ConversePiezo(const ModelOptions& options):BodyForceModel(options)
{
}

void
ConversePiezo::calculate(const Elem* elem, const Point& point)
{
  //Get ElectricField
  RealVectorValue ElField = 0;
  std::vector<Point> p(1);
  p[0] = point;
  std::vector<double> values(3);
  _simul->get_solution(elem, ElFieldID, values, p, true);
  ElField(0) = values[0];
  ElField(1) = values[1];
  ElField(2) = values[2];

  //From V/cm to V/m, and from Pa to GPa
  ElField *= 100.0 * 1e-9;
  //Rotate to crystal system
  Material* mat = get_material();
  const RotatedCrystal&   cr = mat->get_rotated_crystal ();
  ElField = cr.RotMatrix.transpose() * ElField;


  //Compute the converse piezo stress
  RealTensor stress(0);
  stress(0,2) = ElField(0) * _e15;
  stress(2,0) =  stress(0,2);
  stress(2,1) = ElField(1) * _e15;
  stress(1,2) = stress(2,1);
  stress(0,0) = ElField(2) * _e31;
  stress(1,1) = ElField(2) * _e31;
  stress(2,2) = ElField(2) * _e33;

  //Rotate the converse piezo stress
  stress =  cr.RotMatrix * (stress * cr.RotMatrix.transpose());

  //From Pa tp GPa
  //stress *= 1e-9;
  set_stress_source(stress);

  //RealTensor dummy_tens(0);
  //set_strain_source(dummy_tens);

  //RealGradient dummy_grad(0);
  //set_force_source(dummy_grad);

}





void
ConversePiezo::do_init(void)
{

  _e33 = get_option("e33", _e33);
  _e31 = get_option("e31", _e31);
  _e15 = get_option("e15", _e15);

  std::string _sim_name = "null";
  get_parameter("poisson_simulation", _sim_name);
  _simul = SimulationInterface::find_simulation(_sim_name);


  if ( _simul == NULL)
    throw InitFailedException("Could not find " + _sim_name);

  ElFieldID = _simul->get_solution_id("ElField");


}


void
ConversePiezo::read_database(void)
{

  Database& db = get_database();
  db.set_section("piezoelectricity");

  _e33 = db.get("e33", 0.0, true);
  _e31 = db.get("e31", 0.0, true);
  _e15 = db.get("e15", 0.0, true);

//---------------------------------------------------------//
}



