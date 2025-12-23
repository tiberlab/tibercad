/*  
 * This file is part of the tiberCAD module common.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file ConstantThermalConductivity.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */



#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"
#include "ConstantThermalConductivity.h"

#include "tibercad/module/TiberModule.h"


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
