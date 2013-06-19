// $Id: Clamp.C 2338 2011-02-15 21:53:25Z maufder $

#include "Extended.h"

#include "TiberModule.h"



void
Extended::do_init(void)
{
  
  RealTensor H(0);
  RealGradient  R(0);
  double A(0);

  set_is_extended(true);
  
  set_coefficients(H,A,R);

}



void
Extended::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  RealTensor H(0);
  RealGradient  R(0);
  double A(0);

  H(0,0) = _normal(0);
  H(1,1) = _normal(1);
  H(2,2) = _normal(2);
  set_coefficients(H,A,R);

}
