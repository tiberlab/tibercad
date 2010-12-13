// $Id: Clamp.C 2124 2010-10-22 14:00:17Z gromano $

#include "Substrate.h"



// The first string is the class name, the second one
// is the type of the model (here it is a contact model),
// the third one is the specific model implementation.
// The library name will then be contact_dirichlet.so
TIBER_MODULE(Substrate, boundary, substrate)



void Substrate::calculate(const Elem* elem, unsigned int side,
			   const Point& point)
{

  RealTensor H(0);
  RealGradient  R(0);

  double x = _normal(0);
  double y = _normal(1);
  double z = _normal(2);


  H(0,0) = 0.0;
  H(0,1) = -z;
  H(0,2) = y;

  H(1,0) = z;
  H(1,1) = 0.0;
  H(1,2) = -x;

  H(2,0) = -y;
  H(2,1) = x;
  H(2,2) = 0.0;

  H *=1e15;
 
  set_coefficients(H,R);

}
