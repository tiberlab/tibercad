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
 * \file WIUtils.h
 * \brief tiberCAD wateringress module header.
 *
 * \note This file is part of module wateringress.
 */

#ifndef TC_WIUTILS_H
#define TC_WIUTILS_H

#include "tibercad/base/tiber_dll.h"

//! A namespace with some utility functions
namespace WIUtils
{

  //! Provide the saturation water pressure using Goff-Gratch equation
  TC_DLEXPORT double goff_gratch(double T);

};

#endif // TC_WIUTILS_H
