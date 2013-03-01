// $Id: Clamp.C 2124 2010-10-22 14:00:17Z gromano $

#include "Custom.h"

#include "TiberModule.h"




void
Custom::do_init(void)
{
  
  set_is_extended(false);

  RealTensor H(0);
  RealGradient  R(0);
  double A(0);
 
  get_parameter("H",H,false);
  get_parameter("R",R,false);

  
  set_coefficients(H,A,R);


}


