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
 * \file Scaling.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/math/Scaling.h"

Scaling::Scaling(void)
  : _type(NONE),
    _potential(1.0),
    _length(1.0),
    _mesh_units(1.0),
    _mobility(1.0),
    _density(1.0)
{
}

Scaling::Scaling(const Scaling& scaling)
  : _type(scaling._type),
    _potential(scaling._potential),
    _length(scaling._length),
    _mesh_units(scaling._mesh_units),
    _mobility(scaling._mobility),
    _density(scaling._density)
{
}

