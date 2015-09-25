// $Id$

#ifndef _TIBERMATH_H_
#define _TIBERMATH_H_

#include "tiber_dll.h"

#include <utility>

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
}

//! Mathematical functions
/*!
 * \note We do not export inlined functions.
 */
namespace TiberMath
{

  //! The power of 2
  inline double pow_2(double x);

  //! The fermi integral of order +1/2
  inline double fermidirac_half(double x);

  //! The fermi integral of order -1/2
  inline double fermidirac_mhalf(double x);

  //! Calculate SVD of a matrix
  //void svd(DenseMatrix<double>& matrix, DenseVector<double>& sigma)
}


//! Distribution functions
namespace Distributions
{

  //! The Fermi-Dirac distribution
  std::pair<double, double> fermi_dirac(double E, double kT);

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



#endif // _TIBERMATH_H_
