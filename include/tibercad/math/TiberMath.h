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
 * \file TiberMath.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _TIBERMATH_H_
#define _TIBERMATH_H_

#include "tibercad/base/tiber_dll.h"

#include <utility>
#include <vector>

namespace {

  //! The fermi integral of order -1/2
  /*!
   * 
   */
  extern "C" double fdm0p5_(double&);


  //! The fermi integral of order +1/2
  /*!
   * 
   */
  extern "C" double fdp0p5_(double&);

  //! The second derivative of fermi integral of order +1/2
  /*!
   *
   */
  extern "C" double d2_fd_(double&);
}

//! Mathematical functions
/*!
 * \note We do not export inlined functions.
 */
namespace TiberMath
{

  //! The power of 2
  double pow_2(double x);

  //! The fermi integral of order +1/2
  double fermidirac_half(double x);

  //! The fermi integral of order -1/2
  double fermidirac_mhalf(double x);

  //! The second derivative of fermi integral of order +1/2
  double d2_fermidirac(double x);


  //! Calculate SVD of a matrix
  //void svd(DenseMatrix<double>& matrix, DenseVector<double>& sigma)
}


//! Distribution functions
namespace Distributions
{

  //! The Fermi-Dirac distribution
  /*!
   * \param E the argument \f$E=E_F-E_0\f$
   * \param kT the thermal energy
   * \return the value and the derivative with respect to \f$E_0\f$
   */
  std::pair<double, double> fermi_dirac(double E, double kT);

  //! The Fermi-Dirac distribution
  /*!
   * \param E the argument \f$E=E_F-E_0\f$
   * \param kT the thermal energy
   * \return the value and the derivatives with respect to \f$E_0\f$
   *
   * \c result will contain the value, its first derivative and
   * its second derivative, depending on the size of the vector
   */
  void fermi_dirac(std::vector<double>& result, double E, double kT);

}


inline
double
TiberMath::pow_2(double x)
{
  return x * x;
}


inline
double
TiberMath::fermidirac_half(double x)
{
  return fdp0p5_(x);
}


inline
double
TiberMath::fermidirac_mhalf(double x)
{
  return fdm0p5_(x);
}

// TODO
inline
double
TiberMath::d2_fermidirac(double x)
{
  return d2_fd_(x);
}

#endif // _TIBERMATH_H_
