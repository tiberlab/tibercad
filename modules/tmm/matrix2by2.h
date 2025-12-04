/*  
 * This file is part of the tiberCAD module tmm.
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
 * \file matrix2by2.h
 * \brief tiberCAD tmm module header.
 *
 * \note This file is part of module tmm.
 */

#ifndef MATRIX2BY2_H
#define MATRIX2BY2_H
#include <complex>


class matrix2by2
{
public:
  matrix2by2();
  matrix2by2(double, double, double, double);
  ~matrix2by2();
  
  void inv();
  void unit_matrix();
  void print();
  matrix2by2 operator*(matrix2by2 const& old_matrix2by2);
  
  void set(int elm, std::complex<double> a);
  std::complex<double> get(int elm);


private:
  std::complex<double> _m00;
  std::complex<double> _m01;
  std::complex<double> _m10;
  std::complex<double> _m11;
};
#endif
