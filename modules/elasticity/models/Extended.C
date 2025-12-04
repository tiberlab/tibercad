/*  
 * This file is part of the tiberCAD module elasticity.
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
 * \file Extended.C
 * \brief tiberCAD elasticity module implementation.
 *
 * \note This file is part of module elasticity.
 */


#include "Extended.h"

#include "tibercad/module/TiberModule.h"



void
Extended::do_init(void)
{
  
  libMesh::RealTensor H(0);
  libMesh::RealGradient  R(0);
  double A(0);

  set_is_extended(true);
  
  set_coefficients(H,A,R);

}



void
Extended::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  libMesh::RealTensor H(0);
  libMesh::RealGradient  R(0);
  double A(0);

  H(0,0) = _normal(0);
  H(1,1) = _normal(1);
  H(2,2) = _normal(2);
  set_coefficients(H,A,R);

}
