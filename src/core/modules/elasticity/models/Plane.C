// $Id: Clamp.C 2124 2010-10-22 14:00:17Z gromano $

#include "Plane.h"

#include "TiberModule.h"



void Plane::calculate(const Elem* elem, unsigned int side,
			   const Point& point)
{

  RealTensor H(0);
  RealGradient  R(0);
  double A(0);
  set_is_extended(false);

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

  set_coefficients(H,A,R);

}
