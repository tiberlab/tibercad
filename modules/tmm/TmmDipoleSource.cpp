/*  
 * This file is part of the tiberCAD module tmm.
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
 * \file TmmDipoleSource.C
 * \brief tiberCAD tmm module implementation.
 *
 * \note This file is part of module tmm.
 */


#include "TmmDipoleSource.h"
#include "tibercad/io/Database.h"
#include "tibercad/io/Messages.h"

#include <boost/filesystem/operations.hpp>

using std::string;
using namespace libMesh;


TmmDipoleSource::TmmDipoleSource(const ModelOptions& options) :
  PhysicalModel(options),
  _emission_power(0.0)
{


}


const double&
TmmDipoleSource::get_emission_power(void) const
{

  return(_emission_power);
}


void
TmmDipoleSource::set_emission_power(const double& emission_power)
{
  _emission_power = emission_power;
}




