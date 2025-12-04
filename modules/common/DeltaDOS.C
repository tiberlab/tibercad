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
 * \file DeltaDOS.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */


#include "DeltaDOS.h"

#include "tibercad/module/TiberModule.h"

using namespace std;

DeltaDOS::DeltaDOS(const ModelOptions& options) :
  DensityOfStates(options),
  _N0(1.0)
{
}


void
DeltaDOS::do_init(void)
{
  get_parameter("level", reference_energy());
  _N0 = get_option("N0", _N0);
  effective_mass()[0] = 1.0;
  effective_dos() = _N0;
  total_state_density() = _N0;
}

void
DeltaDOS::calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
    double kT, double , const Elem* , const Point& ) const
{

  double dens, der, der2;

  double expf = exp(-(Ef - get_reference_energy()[0] - Epot) / kT);
  dens = _N0 / (1.0 + expf);

  if (result.size() > 1)
  {
    der = dens;
    der /= kT * (1.0 + expf);
    der *= expf;
  }
  if (result.size() > 1)
  {
    der2 = der + dens/kT * expf / (1.0 + expf);
    der2 *= expf;
    der2 /= kT * (1.0 + expf);
  }

  result[0] = dens;
  if (result.size() > 1)
    result[1] = der;
  if (result.size() > 2)
    result[2] = der2;
}
