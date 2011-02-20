// $Id: Clamp.C 2338 2011-02-15 21:53:25Z maufder $

#include "Extended.h"



// The first string is the class name, the second one
// is the type of the model (here it is a contact model),
// the third one is the specific model implementation.
// The library name will then be contact_dirichlet.so
TIBER_MODULE(Extended, ebnd, extended)



void
Extended::do_init(void)
{
  
  RealTensor H(0);
  RealGradient  R(0);
  double A = 1;

  set_coefficients(H,R,A);
}
