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
 * \file TypeDefs.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include "libmesh/id_types.h"
#include <limits.h>


//! A typedef to be used for numerical identifiers
/*!
 * We use the same type as the libMesh subdomain id type
 */
typedef libMesh::subdomain_id_type ID;

#ifndef INVALID_ID
# define INVALID_ID std::numeric_limits<ID>::max()
#endif



//! To ignore unused variables
/*!
 * If you have a variable \a which is unused you can call
 * \c ignore_unused_variable(a) so that the compiler gives no warning
 */
template <typename T>
inline
void ignore_unused_variable(T)
{
}

#endif // _TYPEDEFS_H_
