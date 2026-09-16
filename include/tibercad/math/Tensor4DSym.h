/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file Tensor4DSym.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef TC_TENSOR4DSYM_H
#define TC_TENSOR4DSYM_H

#include <array>
#include <cassert>

class Tensor2;

/*!
 * \brief A simple symmetric rank-4 tensor class
 *
 * The class represents a rank-4 tensor which is symmetric
 * in both index pairs: A_ijkl = A_klij = A_jikl
 */
class Tensor4DSym 
{
  public:

    //! Default constructor
    Tensor4DSym(void) = default;

    //! Constructor to set a specific form
      explicit Tensor4DSym(double x);

      //! Read-write access elements
      double& operator()(unsigned int i, unsigned int j,
                         unsigned int k, unsigned int l);
      

      //! Read-only access
      double operator()(unsigned int i, unsigned int j,
                        unsigned int k, unsigned int l) const;
      

      Tensor4DSym& operator=(double x);

      Tensor4DSym& operator+=(const Tensor4DSym &A);

      Tensor4DSym& operator-=(const Tensor4DSym &A);

      Tensor4DSym& operator*=(double s);

      Tensor4DSym& operator/=(double s);

      //! Apply coordinate transformation
      void push_forward(const Tensor2& F);


  private:

    //! The internal data structure
    std::array<double, 21> _data = {};

};


inline
double&
Tensor4DSym::operator()(unsigned int i, unsigned int j,
                        unsigned int k, unsigned int l)
{
  assert( i > 0 && i < 4);
  assert( j > 0 && j < 4);
  assert( k > 0 && k < 4);
  assert( l > 0 && l < 4);

  if (j > i)
  {
    unsigned int ii = i;
    i = j;
    j = ii;
  }
  
  if (l > k)
  {
    unsigned int kk = k;
    k = l;
    l = kk;
  }

  unsigned int I;
  switch (i - j)
  {
    case 0:
      I = i-1;
      break;
    case 1:
      I = i+1;
      break;
    case 2:
      I = 5;
      break;

    default:
      break;
  }

  unsigned int J;
  switch (k - l){
    case 0:
      J = k-1;
      break;
    case 1:
      J = k+1;
      break;
    case 2: 
      J = 5;
      break;

    default:
      break;
  }

  if (J > I)
  {
    int II = J;
    J=I;
    I=II;
  }

  switch (I - J){
    case 0:
      return _data[I];
    case 1:
      return _data[I+5];
    case 2:
      return _data[I+9];
    case 3:
      return _data[I+12];
    case 4:
      return _data[I+14];
    case 5:
      return _data[20];
    default:
      return _data[0];
  }
}

inline
double
Tensor4DSym::operator()(unsigned int i, unsigned int j,
                        unsigned int k, unsigned int l) const
{
  assert( i > 0 && i < 4);
  assert( j > 0 && j < 4);
  assert( k > 0 && k < 4);
  assert( l > 0 && l < 4);

  if (j > i)
  {
    unsigned int ii = i;
    i = j;
    j = ii;
  }
  
  if (l > k)
  {
    unsigned int kk = k;
    k = l;
    l = kk;
  }

  unsigned int I;
  switch (i - j)
  {
    case 0:
      I = i-1;
      break;
    case 1:
      I = i+1;
      break;
    case 2:
      I = 5;
      break;

    default:
      break;
  }

  unsigned int J;
  switch (k - l){
    case 0:
      J = k-1;
      break;
    case 1:
      J = k+1;
      break;
    case 2: 
      J = 5;
      break;

    default:
      break;
  }

  if (J > I)
  {
    int II = J;
    J=I;
    I=II;
  }

  switch (I - J){
    case 0:
      return _data[I];
    case 1:
      return _data[I+5];
    case 2:
      return _data[I+9];
    case 3:
      return _data[I+12];
    case 4:
      return _data[I+14];
    case 5:
      return _data[20];
    default:
      return _data[0];
  }
}


inline
Tensor4DSym&
Tensor4DSym::operator+=(const Tensor4DSym &A)
{
  for (unsigned int i=0; i<21; i++ )
    _data[i] += A._data[i];

  return *this;
}

inline
Tensor4DSym&
Tensor4DSym::operator-=(const Tensor4DSym &A)
{
  for (unsigned int i=0; i<21; i++ )
    _data[i] -= A._data[i];

  return *this;
}

inline
Tensor4DSym&
Tensor4DSym::operator*=(double s)
{
  for (unsigned int i=0; i<21; i++ )
    _data[i] *= s;

  return *this;
}

inline
Tensor4DSym&
Tensor4DSym::operator/=(double s)
{
  double m = 1.0/s;
  for (unsigned int i=0; i<21; i++ )
    _data[i] *= m;

  return *this;
}


#endif // TC_TENSOR4DSYM_H
