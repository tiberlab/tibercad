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
 * \file ConstantHeatSource.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */


#include "ConstantHeatSource.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"




using namespace std;


ConstantHeatSource::ConstantHeatSource(const ModelOptions& options):HeatSourceModel(options)
{
}

void
ConstantHeatSource::do_init(void)
{
  double heat_source = 0.0;
 
  get_parameter("H",heat_source);  

  set_heat_source(heat_source);

}




 
