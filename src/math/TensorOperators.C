// $Id$

#include "tibercad/math/TensorOperators.h"

#include <iomanip>



libMesh::RealTensorValue doubleContraction(const Tensor4DSym& A, const libMesh::RealTensorValue& B)
{
  // we just call the Tensor implementation

  Tensor2 T;
  T(1, 1) = B(0, 0); T(1, 2) = B(0, 1); T(1, 3) = B(0, 2);
  T(2, 1) = B(1, 0); T(2, 2) = B(1, 1); T(2, 3) = B(1, 2);
  T(3, 1) = B(2, 0); T(3, 2) = B(2, 1); T(3, 3) = B(2, 2);

  T = doubleContraction(A, T);

  libMesh::RealTensorValue R;
  R(0, 0) = T(1, 1);
  R(1, 0) = R(0, 1) = T(2, 1);
  R(1, 1) = T(2, 2);
  R(2, 0) = R(0, 2) = T(3, 1);
  R(2, 1) = R(1, 2) = T(3, 2);
  R(2, 2) = T(3, 3);

  return R;
}



void doubleContraction(const Tensor4DSym& A, libMesh::RealTensorValue& B)
{
  // we just call the Tensor implementation

  Tensor2 T;
  T(1, 1) = B(0, 0); T(1, 2) = B(0, 1); T(1, 3) = B(0, 2);
  T(2, 1) = B(1, 0); T(2, 2) = B(1, 1); T(2, 3) = B(1, 2);
  T(3, 1) = B(2, 0); T(3, 2) = B(2, 1); T(3, 3) = B(2, 2);

  T = doubleContraction(A, T);

  B(0, 0) = T(1, 1);
  B(1, 0) = B(0, 1) = T(2, 1);
  B(1, 1) = T(2, 2);
  B(2, 0) = B(0, 2) = T(3, 1);
  B(2, 1) = B(1, 2) = T(3, 2);
  B(2, 2) = T(3, 3);
}


Tensor2 doubleContraction(const Tensor4DSym& A, const Tensor2& B)
{
  Tensor2 R;

  for (unsigned int i = 1; i <= 3; i++)
  {
    for (unsigned int j = 1; j <= 3; j++)
    {
      double r = 0;
      for (unsigned int k = 1; k <= 3; k++)
        for (unsigned int l = 1; l <= 3; l++)
          r += A(i, j, k, l) * B(k, l);

      R(i, j) = r;
    }
  }

  return R;
}


Tensor4DSym push_forward(const Tensor4DSym& A, const Tensor2& F)
{
  Tensor4DSym T(A);
  T.push_forward(F);

  return T;
}


std::ostream& operator<<(std::ostream &s, const Tensor1 &T)
{
  int w = s.width();
  s << std::setw(w) << T(1) << " " << std::setw(w) << T(2) << " " << std::setw(w) << T(3) << std::endl;

  return s;
}

std::ostream& operator<<(std::ostream &s, const Tensor2 &T)
{
  int w = s.width();
  s << std::setw(w) << T(1,1) << " " << std::setw(w) << T(1,2) << " " << std::setw(w) << T(1,3) << std::endl;
  s << std::setw(w) << T(2,1) << " " << std::setw(w) << T(2,2) << " " << std::setw(w) << T(2,3) << std::endl;
  s << std::setw(w) << T(3,1) << " " << std::setw(w) << T(3,2) << " " << std::setw(w) << T(3,3) << std::endl;

  return s;
}

std::ostream& operator<<(std::ostream &s, const Tensor4DSym &T)
{
  // use Voigt notation:
  // 11 -> 1, 22 -> 2, 33 ->3, 23 -> 4, 13 -> 5, 12 -> 6
  int w = s.width();
  s << std::setw(w) << T(1,1,1,1) << " " << std::setw(w) << T(1,1,2,2) << " " << std::setw(w) << T(1,1,3,3) << " ";
  s << std::setw(w) << T(1,1,2,3) << " " << std::setw(w) << T(1,1,1,3) << " " << std::setw(w) << T(1,1,1,2) << std::endl;
  s << std::setw(w) << T(2,2,1,1) << " " << std::setw(w) << T(2,2,2,2) << " " << std::setw(w) << T(2,2,3,3) << " ";
  s << std::setw(w) << T(2,2,2,3) << " " << std::setw(w) << T(2,2,1,3) << " " << std::setw(w) << T(2,2,1,2) << std::endl;
  s << std::setw(w) << T(3,3,1,1) << " " << std::setw(w) << T(3,3,2,2) << " " << std::setw(w) << T(3,3,3,3) << " ";
  s << std::setw(w) << T(3,3,2,3) << " " << std::setw(w) << T(3,3,1,3) << " " << std::setw(w) << T(3,3,1,2) << std::endl;
  s << std::setw(w) << T(2,3,1,1) << " " << std::setw(w) << T(2,3,2,2) << " " << std::setw(w) << T(2,3,3,3) << " ";
  s << std::setw(w) << T(2,3,2,3) << " " << std::setw(w) << T(2,3,1,3) << " " << std::setw(w) << T(2,3,1,2) << std::endl;
  s << std::setw(w) << T(1,3,1,1) << " " << std::setw(w) << T(1,3,2,2) << " " << std::setw(w) << T(1,3,3,3) << " ";
  s << std::setw(w) << T(1,3,2,3) << " " << std::setw(w) << T(1,3,1,3) << " " << std::setw(w) << T(1,3,1,2) << std::endl;
  s << std::setw(w) << T(1,2,1,1) << " " << std::setw(w) << T(1,2,2,2) << " " << std::setw(w) << T(1,2,3,3) << " ";
  s << std::setw(w) << T(1,2,2,3) << " " << std::setw(w) << T(1,2,1,3) << " " << std::setw(w) << T(1,2,1,2) << std::endl;

  return s;
}
