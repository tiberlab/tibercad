// $Id$

#include "TensorOperators.h"



RealTensorValue doubleContraction(const Tensor4DSym& A, const RealTensorValue& B)
{
  // we just call the Tensor implementation

  Tensor2Sym T;
  T(1, 1) = B(0, 0);
  T(2, 1) = B(1, 0); T(2, 2) = B(1, 1);
  T(3, 1) = B(2, 0); T(3, 2) = B(2, 1); T(3, 3) = B(2, 2);

  T = doubleContraction(A, T);

  RealTensorValue R;
  R(0, 0) = T(1, 1);
  R(1, 0) = R(0, 1) = T(2, 1);
  R(1, 1) = T(2, 2);
  R(2, 0) = R(0, 2) = T(3, 1);
  R(2, 1) = R(1, 2) = T(3, 2);
  R(2, 2) = T(3, 3);

  return R;
}



void doubleContraction(const Tensor4DSym& A, RealTensorValue& B)
{
  // we just call the Tensor implementation

  Tensor2Sym T;
  T(1, 1) = B(0, 0);
  T(2, 1) = B(1, 0); T(2, 2) = B(1, 1);
  T(3, 1) = B(2, 0); T(3, 2) = B(2, 1); T(3, 3) = B(2, 2);

  T = doubleContraction(A, T);

  B(0, 0) = T(1, 1);
  B(1, 0) = B(0, 1) = T(2, 1);
  B(1, 1) = T(2, 2);
  B(2, 0) = B(0, 2) = T(3, 1);
  B(2, 1) = B(1, 2) = T(3, 2);
  B(2, 2) = T(3, 3);
}
