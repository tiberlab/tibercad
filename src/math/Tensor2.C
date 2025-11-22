#include "tibercad/math/Tensor2.h"

#include <algorithm>
#include <cassert>
#include <cmath>

Tensor2
Tensor2::inverse(void) const
{
  Tensor2 T(*this);

  T.invert();
  return T;
}

void
Tensor2::invert(void)
{
   double sub1 = _data[4]*_data[8] - _data[7]*_data[5];
   double sub2 = _data[3]*_data[8] - _data[6]*_data[5];
   double sub3 = _data[3]*_data[7] - _data[6]*_data[4];
   double detA = _data[0]*sub1 - _data[1]*sub2 + _data[2]*sub3;
   
   assert(std::fabs(detA) > 1.e-15);

   Tensor2 r;

   r._data[0] =  sub1/detA;
   r._data[1] = (-_data[1]*_data[8]+_data[2]*_data[7])/detA;
   r._data[2] = ( _data[1]*_data[5]-_data[2]*_data[4])/detA;
   r._data[3] = -sub2/detA;
   r._data[4] = ( _data[0]*_data[8]-_data[2]*_data[6])/detA;
   r._data[5] = (-_data[0]*_data[5]+_data[2]*_data[3])/detA;
   r._data[6] =  sub3/detA;
   r._data[7] = (-_data[0]*_data[7]+_data[1]*_data[6])/detA;
   r._data[8] = ( _data[0]*_data[4]-_data[1]*_data[3])/detA;

   _data.swap(r._data);
}

void
Tensor2::inplace_transpose(void)
{
  std::swap(_data[1], _data[3]);
  std::swap(_data[2], _data[6]);
  std::swap(_data[5], _data[7]);
}


double
Tensor2::det(void) const
{
  auto& T = _data;
  double d = T[0] * (T[4]*T[8] - T[7]*T[5])
           - T[1] * (T[3]*T[8] - T[6]*T[5])
           + T[2] * (T[3]*T[7] - T[6]*T[4]);

  return d;
}


double
Tensor2::norm(void) const
{
  double d = 0.0;
  for (unsigned int i = 0; i < 9; ++i)
    d += _data[i]*_data[i];

  return std::sqrt(d);
}
