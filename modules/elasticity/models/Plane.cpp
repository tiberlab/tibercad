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
 * \file Plane.C
 * \brief tiberCAD elasticity module implementation.
 *
 * \note This file is part of module elasticity.
 */


#include "Plane.h"

#include "tibercad/module/TiberModule.h"



void Plane::calculate(const libMesh::Elem* elem, unsigned int side,
			   const libMesh::Point& point)
{

  libMesh::RealTensor H(0);
  libMesh::RealGradient  R(0);
  double A(0);
  set_is_extended(false);

  double x = _normal(0);
  double y = _normal(1);
  double z = _normal(2);

  H(0,0) = _normal(0);
  H(1,1) = _normal(1);
  H(2,2) = _normal(2);

  /*
  H(0,0) = 0.0;
  H(0,1) = -z;
  H(0,2) = y;

  H(1,0) = z;
  H(1,1) = 0.0;
  H(1,2) = -x;

  H(2,0) = -y;
  H(2,1) = x;
  H(2,2) = 0.0;
  */

  set_coefficients(H,A,R);

}
