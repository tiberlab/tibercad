// $Id$


#include "Material.h"
#include "Database.h"
#include "ConstantThermalConductivity.h"
#include "RotatedCrystal.h"


TIBER_MODULE(ConstantThermalConductivity, ThermalConductivity, constant)

using namespace std;


ConstantThermalConductivity::ConstantThermalConductivity(const ModelOptions& options):ThermalConductivityModel(options)
{
 
}

void
ConstantThermalConductivity::read_database(void)
{

  Database& db = get_database();
  db.set_section("thermal_conductivity/constant");

  _kappa(0);
  db.get("ThermCond", _kappa,false);

}

void
ConstantThermalConductivity::do_init(void)
{

  get_parameter("ThermCond",_kappa,true);

  set_thermal_conductivity(_kappa);

  rotate();

 //  if (get_material()->get_structure() == "wz")
//   {
//     const RotatedCrystal&   cr = get_material()->get_rotated_crystal();
//     rotate_to_calculation_system(cr.RotMatrix);
//   }

}

