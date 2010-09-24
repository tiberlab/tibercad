// $Id$

#ifndef _TENSOROPERATORS_H_
#define _TENSOROPERATORS_H_


#include <tensor.h>
#include <tensor_value.h>


//! Multiplication of symmetric 2-rank tensor with vector
inline
RealVectorValue operator*(const Tensor2Sym& A, const RealVectorValue& x)
{
  RealVectorValue v(0);
  for (unsigned int i = 0; i < 3; ++i)
  {
    v(0) += A(i+1, 1) * x(i);
    v(1) += A(i+1, 2) * x(i);
    v(2) += A(i+1, 3) * x(i);
  }

  return v;
}


#endif /* _TENSOROPERATORS_H_ */
