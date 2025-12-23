/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file WzDDsemiconductor.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/semiconductormodels/WzDDsemiconductor.h"
#include "tibercad/physics/semiconductormodels/WzSemiconductor.h"
#include "tibercad/physics/Constants.h"

using namespace std;
using namespace Constants;







//-------------------------------------------------------------/
void WzDDsemiconductor::do_calculate_conduction_band_extremum(void)
{
  vector<DDsemiconductor::band_extremum> result;

  const WzSemiconductor::WzDDparameters& par = (dynamic_cast<WzSemiconductor*> (semiconductor))->get_parameters ();


  double energy = par.Ev + par.EgGamma;

  if (strained) energy += (strain(1,1) + strain(2,2))* par.a_x + par.a_z *  strain(3,3);

  DDsemiconductor::band_extremum band_ext;
  band_ext.energy = energy;
  band_ext.degeneracy = 2;
  band_ext.mass_DOS = pow(par.m_c_xx * par.m_c_xx * par.m_c_zz, 1.0/3.0) ;

  result.push_back(band_ext);

  conduction_band = result;
}




