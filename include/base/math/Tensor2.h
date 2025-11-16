#ifndef _TENSOR2_H_
#define _TENSOR2_H_

#include "Tensor1.h"

#include <array>
#include <cassert>

/*!
 * \brief Small class for rank-2 tensors
 *
 * Internal data structure is column-major, access
 * is using Fortran-style 1-based indexing.
 */
class Tensor2
{

  public:

    //! Default constructor
    Tensor2(void) = default;

    //! Create diagonal tensor
    explicit Tensor2(double x);

    //! Copy constructor
    Tensor2(const Tensor2& t) = default;

    //! Destructor
    ~Tensor2(void) = default;

    //! Get values
    double operator()(unsigned int i, unsigned int j) const;
    
    //! Get or set values
    double& operator()(unsigned int i, unsigned int j);
 
    //! Set to diagonal tensor with value x
    Tensor2& operator=(double x);

    //! Add another Tensor
    Tensor2& operator+=(const Tensor2& b);

    //! Subtract another Tensor
    Tensor2& operator-=(const Tensor2& B);

    //! Multiply with a scalar
    Tensor2& operator*=(double s);

    //! Divie by a scalar
    Tensor2& operator/=(double s);

    //! Calculate an return the inverse
    Tensor2 inverse(void) const;

    //! Get the transpose
    Tensor2 transpose(void) const;

    //! Get the determinant
    double det(void) const;

    //! Get the trace
    double trace(void) const;

    //! Get the norm
    double norm(void) const;

    //! Extract a column as Tensor1
    Tensor1 operator()(unsigned int i) const;

    //! Inplace inversion
    void invert(void);

    //! Inplace transpose
    void inplace_transpose(void);


  protected:

    //! The internal data structure
    /*!
     * Order is: 11, 12, 13, 21, 22, 23, 31, 32, 33
     */
    std::array<double, 9> _data = {};

};

typedef Tensor2 Tensor2Gen;


inline
Tensor2::Tensor2(double x)
{
  _data[0] = _data[4] = _data[8] = x;
}


inline
Tensor2&
Tensor2::operator=(double x)
{
  _data[0] = _data[4] = _data[8] = x;
  _data[1] = _data[2] = 0;
  _data[3] = _data[5] = 0;
  _data[6] = _data[7] = 0;

  return *this;
}

inline
Tensor2&
Tensor2::operator+=(const Tensor2 &b)
{
  for (int i = 0; i < 9; i++)
    _data[i] += b._data[i];

  return *this;
}

inline
Tensor2&
Tensor2::operator-=(const Tensor2 &b)
{
  for (int i = 0; i < 9; i++)
    _data[i] -= b._data[i];

  return *this;
}

inline
Tensor2&
Tensor2::operator*=(double s)
{
  for (int i = 0; i < 9; i++)
    _data[i] *= s;

  return *this;
}

inline
Tensor2&
Tensor2::operator/=(double s)
{
  return (this->operator*=(1.0 / s));
}

inline
double&
Tensor2::operator()(unsigned int i, unsigned int j)
{
  assert((i > 0) && (i < 4));
  assert((j > 0) && (j < 4));
  return _data[(j-1)*3 + i-1];
}

inline
double
Tensor2::operator()(unsigned int i, unsigned int j) const
{
  assert((i > 0) && (i < 4));
  assert((j > 0) && (j < 4));
  return _data[(j-1)*3 + i-1];
}


inline
Tensor1
Tensor2::operator()(unsigned int i) const
{
  Tensor1 t;
  t(1) = (*this)(1, i);
  t(2) = (*this)(2, i);
  t(3) = (*this)(3, i);

  return t;
}

inline
Tensor2
Tensor2::transpose() const
{
  Tensor2 t(*this);
  t.inplace_transpose();

  return t;
}

inline
double
Tensor2::trace(void) const
{
  double tr = _data[0] + _data[4] + _data[8];

  return tr;
}



#endif // _TENSOR2_H_

