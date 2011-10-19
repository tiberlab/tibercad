// $Id: ThermalStress.C 2451 2011-03-05 14:45:46Z maufder $

#include "ThermalStress.h"
#include "Material.h"
#include "RotatedCrystal.h"
#include "Database.h"
#include "SimulationOptions.h"
#include "TensorOperators.h"

#include "TiberModule.h"


using namespace std;


ThermalStress::ThermalStress(const ModelOptions& options) :
    BodyForceModel(options)
{
}

ThermalStress::~ThermalStress(void)
{

}


void
ThermalStress::do_init(void)
{

  RealGradient body_force(0);
  
 //Get reference lattice
  _ref_temp = get_option("reference_temperature", SimulationOptions::temperature);
 
  get_parameter("thermal_expansion_coefficient", _alpha);

 
  std::string temp_sim = get_option("thermal_simulation", "");

  _temp.set_simulation(temp_sim);
}



void
ThermalStress::read_database(void)
{

  const Database& db = get_database();
  db.set_section("lattice");

  vector<double> alpha;
  db.get("thermal_coefficient", alpha);

  switch (alpha.size())
  {
    case 1:
      _alpha(0) = _alpha(1) = _alpha(2) = alpha[0];
      break;
    case 2:
      _alpha(0) = _alpha(1) = alpha[0];
      _alpha(2) = alpha[1];
      break;
    case 3:
      _alpha(0) = alpha[0];
      _alpha(1) = alpha[1];
      _alpha(2) = alpha[2];
      break;
    default:
      break;
  }

}

 

void
ThermalStress::calculate(const Elem* elem, const Point& point)
{
  // get temperature
  double deltaT = _temp.get_temperature(elem, point) - _ref_temp;

  const RotatedCrystal& cr = get_material()->get_rotated_crystal();


  // compute thermally induced strain
  RealTensor strain(0);
  strain(0,0) = -_alpha(0) * deltaT;
  strain(1,1) = -_alpha(1) * deltaT;
  strain(2,2) = -_alpha(2) * deltaT;

  // rotate
  strain = cr.RotMatrix * (strain * cr.RotMatrix.transpose());

  set_strain_source(strain);
}

