// $Id: ConstantThermalConductivity.C 2457 2011-03-06 23:52:12Z gromano $


#include "Material.h"
#include "Database.h"
#include "ConstantThermalConductivity.h"

#include "TiberModule.h"


using namespace std;


ConstantThermalConductivity::ConstantThermalConductivity(const ModelOptions& options) :
    ThermalConductivityModel(options),
    _kappa(0),
    _temp_coeff(0.0),
    _ref_temp(300)
{
 
}

void
ConstantThermalConductivity::read_database(void)
{

  const Database& db = get_database();
  db.set_section("thermal_conductivity/constant");

  db.get("ThermCond", _kappa, false);
  db.get("TempCoeff", _temp_coeff, false);
  db.get("RefTemp", _ref_temp, false);

}

void
ConstantThermalConductivity::do_init(void)
{

  get_parameter("ThermCond", _kappa, true);

  get_parameter("TempCoeff", _temp_coeff);

  set_thermal_conductivity(_kappa);

  rotate();



}


void
ConstantThermalConductivity::calculate(const Elem* elem, const Point& point, double temperature)
{
  libMesh::RealGradient k(_kappa);
  if (((_temp_coeff(0) != 0) || (_temp_coeff(1) != 0) ||
       (_temp_coeff(2) != 0)) && (temperature != _ref_temp))
  {
    k(0) *= pow(_ref_temp / temperature, _temp_coeff(0));
    k(1) *= pow(_ref_temp / temperature, _temp_coeff(1));
    k(2) *= pow(_ref_temp / temperature, _temp_coeff(2));

    set_thermal_conductivity(k);
    rotate();
  }

}
