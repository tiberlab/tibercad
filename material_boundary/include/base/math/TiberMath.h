// $Id$

#ifndef _TIBERMATH_H_
#define _TIBERMATH_H_

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


namespace TiberCad
{

  //! The power of 2
  extern inline double pow_2(double x) { return x * x; }

  //! The fermi integral of order +1/2
  extern inline double fermidirac_half(double x);

  //! The fermi integral of order -1/2
  extern inline double fermidirac_mhalf(double x);

}


extern inline
double
TiberCad::fermidirac_half(double x)
{
  return fdp0p5_(x);
}


extern inline
double
TiberCad::fermidirac_mhalf(double x)
{
  return fdm0p5_(x);
}



#endif // _TIBERMATH_H_
