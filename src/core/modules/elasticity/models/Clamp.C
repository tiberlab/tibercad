// $Id$

#include "Clamp.h"



// The first string is the class name, the second one
// is the type of the model (here it is a contact model),
// the third one is the specific model implementation.
// The library name will then be contact_dirichlet.so
TIBER_MODULE(Clamp, Boundary, clamp )



void
Clamp::do_init(void)
{
 
  RealTensor H(0);
  RealGradient  R(0);


  H(0,0) = 1;
  H(1,1) = 1;
  H(2,2) = 1;

  H *=1e25;
  set_coefficients(H,R);
}
