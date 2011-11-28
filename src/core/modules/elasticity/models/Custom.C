// $Id: Clamp.C 2124 2010-10-22 14:00:17Z gromano $

#include "Custom.h"



// The first string is the class name, the second one
// is the type of the model (here it is a contact model),
// the third one is the specific model implementation.
// The library name will then be contact_dirichlet.so
TIBER_MODULE(Custom, ebnd, custom)


void
Custom::do_init(void)
{
  
  set_is_extended(false);

  RealTensor H(0);
  RealGradient  R(0);
  double A(0);
 
  get_parameter("H",H,false);
  get_parameter("R",R,false);

  
  set_coefficients(H*1e20,A,R*1e20);


}


