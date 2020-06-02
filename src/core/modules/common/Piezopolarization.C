// $Id$

#include "Piezopolarization.h"
#include "SimulationInterface.h"
#include "Database.h"
#include "Material.h"
#include "TensorOperators.h"

#include "TiberModule.h"


using namespace std;


void
Piezopolarization::do_init(void)
{
  PolarizationModel::do_init();

  std::string sim_name = "";
  get_parameter("strain_simulation", sim_name);
  _strain = SimulationInterface::find_solution_provider(sim_name, "Strain");

}


void
Piezopolarization::read_database(void)
{
  // get piezoelectric_coefficients
  const Database& db = get_database();
  db.set_section("piezoelectricity");
      
  if (get_material()->get_structure() == "wz")
  {
    _e33 = db.get("e33", 0.0);
    _e31 = db.get("e31", 0.0);
    _e15 = db.get("e15", 0.0);
  } 
  else if (get_material()->get_structure() == "zb")
    _e14 = db.get("e14", 0.0);

}


void
Piezopolarization::do_calculate(const Elem* elem, const Point& point)
{  

  RealVectorValue polarization(0);

  Tensor2Sym& strain = get_strain();
  if (_strain.is_valid())
  {
    std::vector<Point> p(1);
    p[0] = point;
    std::vector<double> values(6);

    if (_strain.simulation()->get_solution(elem, _strain.id(), values, p))
    {
      strain(1,1) = values[0];
      strain(2,2) = values[1];
      strain(3,3) = values[2];
      strain(2,1) = values[3];
      strain(3,2) = values[4];
      strain(3,1) = values[5];
    }

    const Material* mat = get_material();
    const RotatedCrystal& cr = mat->get_rotated_crystal();
    const Tensor2Gen& rotate = cr.RotMatrix;
    strain = multAtSA(rotate, strain);
  }
  

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
