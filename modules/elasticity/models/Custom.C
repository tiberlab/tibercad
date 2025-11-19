// $Id: Clamp.C 2124 2010-10-22 14:00:17Z gromano $

#include "Custom.h"

#include "tibercad/module/TiberModule.h"




void
Custom::do_init(void)
{
  
  set_is_extended(false);

  libMesh::RealTensor H(0);
  libMesh::RealGradient  R(0);
  double A(0);
 
  get_parameter("H",H,false);
  get_parameter("R",R,false);

  
  set_coefficients(H,A,R);


}


