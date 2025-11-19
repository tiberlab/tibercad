// $Id$

#ifndef _TENSOROPERATORS_H_
#define _TENSOROPERATORS_H_


#include "tibercad/math/Tensor1.h"
#include "tibercad/math/Tensor2.h"
#include "tibercad/math/Tensor4DSym.h"

#include <libmesh/vector_value.h>
#include <libmesh/tensor_value.h>

#include <iostream>

//! Multiplication of generic 2-rank tensor with vector
inline
libMesh::RealVectorValue operator*(const Tensor2& A, const libMesh::RealVectorValue& x)
{
  libMesh::RealVectorValue v(0);
  for (unsigned int i = 0; i < 3; ++i)
  {
    v(0) += A(1, i+1) * x(i);
    v(1) += A(2, i+1) * x(i);
    v(2) += A(3, i+1) * x(i);
  }

  return v;
}


//! Multiplication of Tensor types
inline
libMesh::RealTensor operator*(const libMesh::RealTensor& A, const Tensor2& B)
{
  libMesh::RealTensor R(0);
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


//! Multiplication of Tensor types
inline
libMesh::RealTensor operator*(const Tensor2& A, const libMesh::RealTensor& B)
{
  libMesh::RealTensor R(0);
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


//! Transform RealTensor into Tensor2
inline
void transform_tensor_format(const libMesh::RealTensor& in, Tensor2& out)
{
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      out(i+1, j+1) = in(i, j);
}


//! Double contraction \f$A_{ijkl}B_{kl}\f$
libMesh::RealTensorValue doubleContraction(const Tensor4DSym& A, const libMesh::RealTensorValue& B);

//! Double contraction \f$A_{ijkl}B_{kl}\f$
/*!
 * The result is written to B
 */
void doubleContraction(const Tensor4DSym& A, libMesh::RealTensorValue& B);

//! Double contraction \f$A_{ijkl}B_{kl}\f$
Tensor2 doubleContraction(const Tensor4DSym& A, const Tensor2& B);

inline
libMesh::RealTensorValue operator*(const Tensor4DSym& A, const libMesh::RealTensorValue& B)
{
  return doubleContraction(A, B);
}

//! norm of Tensor1
inline
double norm(const Tensor1& t)
{
  return t.norm();
}

//! norm of Tensor2
inline
double norm(const Tensor2& t)
{
  return t.norm();
}

//! det of Tensor2
inline
double det(const Tensor2& t)
{
  return t.det();
}

//! trace of Tensor2
inline
double trace(const Tensor2& t)
{
  return t.trace();
}

//! Inverse of Tensor2
inline
Tensor2 inv(const Tensor2& t)
{
  Tensor2 r(t);
  r.invert();

  return r;
}


// arithmetic operators

inline
Tensor1 operator+(const Tensor1& a, const Tensor1& b)
{
  Tensor1 r(a);
  return r += b;
}

inline
Tensor1 operator-(const Tensor1& a, const Tensor1& b)
{
  Tensor1 r(a);
  return r -= b;
}

inline
Tensor1 operator*(const Tensor1& a, double b)
{
  Tensor1 r(a);
  return r *= b;
}

inline
Tensor1 operator*(double b, const Tensor1& a)
{
  Tensor1 r(a);
  return r *= b;
}

inline
Tensor1 operator/(const Tensor1& a, double b)
{
  Tensor1 r(a);
  return r /= b;
}

inline
Tensor1 operator-(const Tensor1& a)
{
  return -1.0*a;
}

inline
Tensor1 operator+(const Tensor1& a)
{
  return a;
}


inline
double operator*(const Tensor1& a, const Tensor1& b)
{
   return (a(1)*b(1) + a(2)*b(2) + a(3)*b(3));
}

inline
Tensor1 operator^(const Tensor1& a, const Tensor1& b)
{
   Tensor1 r;
   r(1) = a(2)*b(3) - a(3)*b(2);
   r(2) = a(3)*b(1) - a(1)*b(3);
   r(3) = a(1)*b(2) - a(2)*b(1);

   return r;
}

inline
Tensor1 vectorProduct(const Tensor1& a, const Tensor1& b)
{
   return a ^ b;
}


inline
Tensor2 operator+(const Tensor2& A, const Tensor2& B)
{
  Tensor2 R(A);
  return R += B;
}

inline
Tensor2 operator-(const Tensor2& A, const Tensor2& B)
{
  Tensor2 R(A);
  return R -= B;
}

inline
Tensor2 operator*(const Tensor2& A, double s)
{
  Tensor2 R(A);
  return R *= s;
}

inline
Tensor2 operator*(double s, const Tensor2& A)
{
  Tensor2 R(A);
  return R *= s;
}

inline
Tensor2 operator/(const Tensor2& A, double s)
{
  Tensor2 R(A);
  return R /= s;
}

inline
Tensor2 operator-(const Tensor2& A)
{
  return -1.0*A;
}

inline
Tensor2 operator+(const Tensor2& A)
{
  return A;
}

inline
Tensor4DSym operator+(const Tensor4DSym& A, const Tensor4DSym& B)
{
  Tensor4DSym R(A);
  return R += B;
}

inline
Tensor4DSym operator-(const Tensor4DSym& A, const Tensor4DSym& B)
{
  Tensor4DSym R(A);
  return R -= B;
}

inline
Tensor4DSym operator*(const Tensor4DSym& A, double s)
{
  Tensor4DSym R(A);
  return R *= s;
}

inline
Tensor4DSym operator*(double s, const Tensor4DSym& A)
{
  Tensor4DSym R(A);
  return R *= s;
}

inline
Tensor4DSym operator/(const Tensor4DSym& A, double s)
{
  Tensor4DSym R(A);
  return R /= s;
}

inline
Tensor1
operator*(const Tensor2 &A, const Tensor1 &b)
{
  Tensor1 r;
  double b1 = b(1);
  double b2 = b(2);
  double b3 = b(3);

  r(1) = A(1, 1) * b1 + A(1, 2) * b2 + A(1, 3) * b3;
  r(2) = A(2, 1) * b1 + A(2, 2) * b2 + A(2, 3) * b3;
  r(3) = A(3, 1) * b1 + A(3, 2) * b2 + A(3, 3) * b3;

  return r;
}


inline
Tensor2
operator*(const Tensor2& A, const Tensor2& B)
{
  Tensor2 L;
  double b;

  for (unsigned int i = 1; i <= 3; i++)
    for (unsigned int j = 1; j <= 3; j++)
    {
      b = 0;
      for (unsigned int k = 1; k <= 3; k++)
        b += A(i, k) * B(k, j);

      L(i, j) = b;
    }

  return L;
}

//! Push forward with linear map F
Tensor4DSym push_forward(const Tensor4DSym& A, const Tensor2& F);


// Writing to output streams

std::ostream& operator<<(std::ostream &s, const Tensor1 &T);

std::ostream& operator<<(std::ostream &s, const Tensor2 &T);

std::ostream& operator<<(std::ostream &s, const Tensor4DSym &T);


#endif // _TENSOROPERATORS_H_
