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
 * \file libMeshDefs.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _LIBMESHDEFS_H_
#define _LIBMESHDEFS_H_


//! Definitions for easier use of libMesh stuff

#define LIBMESHCLASS(name) namespace libMesh { class name; } using libMesh::name
#define USELIBMESHTYPE(name) using libMesh::name

// these are needed almost everywhere
LIBMESHCLASS(Elem);
LIBMESHCLASS(Point);
LIBMESHCLASS(Node);
LIBMESHCLASS(MeshBase);

// they have to correspond to the libMesh typedefs!
namespace libMesh
{
typedef double Real;
typedef Real Number;
}

USELIBMESHTYPE(Real);
USELIBMESHTYPE(Number);

#endif // _LIBMESHDEFS_H_
