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
 * \file FowlerNordheim.C
 * \brief tiberCAD API implementation.
 */



#include "tibercad/physics/misc/FowlerNordheim.h"

#include <cmath>
#include <iostream>

using namespace std;

double
FowlerNordheim::get_emission_current(double F)
{
  double J = 0.0;
  // the formula makes sense only for negative field (E.n)
  if (F < -1e-3)
  {
    double t_square = 1.1164;
    double nu_0 = 0.93685;
    double sqrt_wf = sqrt(_workfunction);
    double A = 1.5414e-6 / (t_square * _workfunction) * exp(9.83596 / sqrt_wf);
    double B = nu_0 * 6.8309e7 * sqrt_wf * sqrt_wf * sqrt_wf;

    // the field is opposite to that assumed in the original formula
    // due to the projection onto the normal
    J = A * F * F * exp(B / F);
  }

  return J;
}


