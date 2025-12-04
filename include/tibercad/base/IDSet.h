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
 * \file IDSet.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _IDSET_H_
#define _IDSET_H_

#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/HashSet.h"
#include <set>

//! A convenient typedef for a set of IDs
typedef std::set<ID> IDSet;

//! A convenient typedef for a hash set of IDs
typedef HashSet<ID>::Type IDHashSet;


#endif // _IDSET_H_
