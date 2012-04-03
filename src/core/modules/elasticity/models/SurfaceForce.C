// $Id$

#include "SurfaceForce.h"

#include "TiberModule.h"





void
SurfaceForce::do_init(void)
{
  
  RealGradient  R(0);
  RealTensor H(0);
  double A(1.0);
  set_is_extended(false);

  get_parameter("force",R,false);

  set_coefficients(H,A,R);

}
