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
 * \file Tensor1.C
 * \brief tiberCAD API implementation.
 */

#include "tibercad/math/Tensor1.h"

#include "libmesh/point.h"

#include <cmath>

Tensor1::Tensor1(const libMesh::Point& p)
{
  _comp[0] = p(0);
  _comp[1] = p(1);
  _comp[2] = p(2);
}

libMesh::Point
Tensor1::get_point(void) const
{
  return(libMesh::Point(_comp[0], _comp[1], _comp[2]));
}

    
double
Tensor1::norm(void) const
{
  double d = _comp[0]*_comp[0] + _comp[1]*_comp[1] + _comp[2]*_comp[2];

  return std::sqrt(d);
}

