/*  
 * This file is part of the tiberCAD module masstransport.
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
 * \file BCPressure.C
 * \brief tiberCAD masstransport module implementation.
 *
 * \note This file is part of module masstransport.
 */


#include "BCPressure.h"
#include "WIUtils.h"
#include "tibercad/base/SimulationOptions.h"

#include "tibercad/module/TiberModule.h"

using namespace libMesh;


void
BCPressure::do_init(void)
{
  _relative_pressure = get_option("relative_pressure", _relative_pressure); 
  _mass_tran_spec = get_option("molecule", _mass_tran_spec);  
}


void
BCPressure::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{

 double pressure = 101324.6;   // Pa, standard atmospheric pressure; 

 if (_mass_tran_spec == "H2O")
 {
    double temp = SimulationOptions::temperature;
    pressure = WIUtils::goff_gratch(temp);
 }

  pressure *= _relative_pressure / 100; 

  set_coefficients(1, 0, pressure);
}


