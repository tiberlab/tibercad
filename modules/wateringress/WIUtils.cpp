/*  
 * This file is part of the tiberCAD module wateringress.
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
 * \file WIUtils.C
 * \brief tiberCAD wateringress module implementation.
 *
 * \note This file is part of module wateringress.
 */

#include "WIUtils.h"

#include <cmath>

double
WIUtils::goff_gratch(double T)
{
  const double Tb = 373.16;     // K, boiling point of water
  const double Pb = 1013.246;   // hPa, standard atmospheric pressure

  double term1 = -7.90298 * (Tb / T - 1.0);
  double term2 = 5.02808 * std::log10(Tb / T);
  double term3 = -1.3816e-7 * (std::pow(10.0, 11.344 * (1.0 - T / Tb)) - 1.0);
  double term4 = 8.1328e-3 * (std::pow(10.0, -3.49149 * (Tb / T - 1.0)) - 1.0);
  double log10_es = term1 + term2 + term3 + term4 + std::log10(Pb);

  return 100 * std::pow(10.0, log10_es);  // e_s in hPa
}
