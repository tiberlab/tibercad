// $Id$

#include "SurfaceForce.h"

#include "tibercad/module/TiberModule.h"





void
SurfaceForce::do_init(void)
{
  
  libMesh::RealGradient  R(0);
  libMesh::RealTensor H(0);
  double A(1.0);
  set_is_extended(false);

  get_parameter("force",R,false);

  set_coefficients(H,A,R);

}
