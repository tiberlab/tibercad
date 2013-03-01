// $Id$

#include "Clamp.h"

#include "TiberModule.h"





void
Clamp::do_init(void)
{
  
  RealTensor H(0);
  RealGradient  R(0);
  double A(0);
  set_is_extended(false);

  H(0,0) = 1;
  H(1,1) = 1;
  H(2,2) = 1;

  set_coefficients(H,A,R);
}
