// $Id$

#include "SurfaceForce.h"



// The first string is the class name, the second one
// is the type of the model (here it is a contact model),
// the third one is the specific model implementation.
// The library name will then be contact_dirichlet.so
TIBER_MODULE(SurfaceForce, boundary, surface_force )



void
SurfaceForce::do_init(void)
{
  
  RealGradient  R(0);
  RealTensor H(0);
  
  get_parameter("force",R,false);

  set_coefficients(H,R);

}
