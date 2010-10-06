// $Id$

#ifndef _TENSOROPERATORS_H_
#define _TENSOROPERATORS_H_


#include <tensor.h>
#include <vector_value.h>
#include <tensor_value.h>


//! Multiplication of symmetric 2-rank tensor with vector
inline
RealVectorValue operator*(const Tensor2Sym& A, const RealVectorValue& x)
{
  RealVectorValue v(0);
  v(0) += A(1, 1) * x(0) + A(2, 1) * x(1) + A(3, 1) * x(2);
  v(1) += A(1, 2) * x(0) + A(2, 2) * x(1) + A(3, 2) * x(2);
  v(2) += A(1, 3) * x(0) + A(2, 3) * x(1) + A(3, 3) * x(2);

  return v;
}


//! Multiplication of generic 2-rank tensor with vector
inline
RealVectorValue operator*(const Tensor2Gen& A, const RealVectorValue& x)
{
  RealVectorValue v(0);
  for (unsigned int i = 0; i < 3; ++i)
  {
    v(0) += A(1, i+1) * x(i);
    v(1) += A(2, i+1) * x(i);
    v(2) += A(3, i+1) * x(i);
  }

  return v;
}


//! Multiplication of Tenor types
inline
RealTensor operator*(const RealTensor& A, const Tensor2Gen& B)
{
  RealTensor R(0);
  for (unsigned int j = 0; j < 3; j++)
  {
    unsigned int J = j + 1;
    for (unsigned int i = 0; i < 3; i++)
    {
      R(i,0) += A(i,j) * B(J,1);
      R(i,1) += A(i,j) * B(J,2);
      R(i,2) += A(i,j) * B(J,3);
    }
  }

  return R;
}


//! Multiplication of Tenor types
inline
RealTensor operator*(const Tensor2Gen& A, const RealTensor& B)
{
  RealTensor R(0);
  for (unsigned int i = 0; i < 3; i++)
  {
    unsigned int I = i + 1;
    for (unsigned int j = 0; j < 3; j++)
    {
      R(i,j) += A(I,1) * B(0,j);
      R(i,j) += A(I,2) * B(1,j);
      R(i,j) += A(I,3) * B(2,j);
    }
  }

  return R;
}


//! Double contraction \f$A_{ijkl}B_{kl}\f$
RealTensorValue doubleContraction(const Tensor4DSym& A, const RealTensorValue& B);

//! Double contraction \f$A_{ijkl}B_{kl}\f$
/*!
 * The result is written to B
 */
RealTensorValue doubleContraction(const Tensor4DSym& A, const RealTensorValue& B);


inline
RealTensorValue operator*(const Tensor4DSym& A, const RealTensorValue& B)
{
  return doubleContraction(A, B);
}


#endif /* _TENSOROPERATORS_H_ */
