// $Id$


#include "Material.h"
#include "Database.h"
#include "ConstantThermalConductivity.h"
#include "RotatedCrystal.h"


TIBER_MODULE(ConstantThermalConductivity, thermal_conductivity, constant)

using namespace std;


ConstantThermalConductivity::ConstantThermalConductivity(const ModelOptions& options):ThermalConductivityModel(options)
{
 
}

void
ConstantThermalConductivity::read_database(void)
{

  Database& db = get_database();
  db.set_section("thermal_conductivity/constant");

  _kappa.resize(3,0.0);
  db.get("ThermCond", _kappa,true);

}

void
ConstantThermalConductivity::do_init(void)
{

  get_parameter("ThermCond",_kappa,true);


  RealTensor kappa_tens(0);
  kappa_tens(0,0) = _kappa[0];
  kappa_tens(1,1) = _kappa[1];
  kappa_tens(2,2) = _kappa[2];

  set_thermal_conductivity(kappa_tens);


  if (get_material()->get_structure() == "wz")
  {
    const RotatedCrystal&   cr = get_material()->get_rotated_crystal();
    rotate_to_calculation_system(cr.RotMatrix);
  }

}

