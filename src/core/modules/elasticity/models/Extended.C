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
