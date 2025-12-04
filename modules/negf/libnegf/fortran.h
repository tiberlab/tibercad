/*  
 * This file is part of the tiberCAD module negf.
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
 * \file fortran.h
 * \brief tiberCAD negf module header.
 *
 * \note This file is part of module negf.
 */

#ifndef _FORTRAN_H
#define _FORTRAN_H
// C++ types corresponding to Fortran 77 types

#include <complex>

typedef int     f77_int;          // integer
typedef float   f77_real;         // real(4)
typedef double  f77_double;       // real(8)
typedef char    f77_char;         // character
typedef std::complex<double> f77_complex;      // complex(8)
//typedef int     f77_logical;      // not used, use f77_int instead

#endif
