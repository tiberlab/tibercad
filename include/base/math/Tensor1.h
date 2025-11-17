#ifndef _TENSOR1_H_
#define _TENSOR1_H_

#include <cassert>
#include <array>

namespace libMesh
{
  class Point;
}

/*!
 * \brief A small class to represent rank-1 tensors
 *
 * Access is using Fortran-style 1-based indexing.
 */
class Tensor1
{

  public:

    //! Default constructor
    Tensor1(void) = default;

    //! Copy constructur
    Tensor1(const Tensor1& t) = default;

    //! Construct with a given value
    explicit Tensor1(double x)
    {
      _comp = {x, x, x};
    }

    //! Construct Tensor1 from a Point
    explicit Tensor1(const libMesh::Point& p);

    //! Access elements
    double& operator()(unsigned int i)
    {
      assert((i > 0) && (i < 4));
      return _comp[i-1];
    }

    //! Const access elements
    double operator()(unsigned int i) const
    {
      assert((i > 0) && (i < 4));
      return _comp[i-1];
    }

    //! Set to a fixed value
    Tensor1& operator=(double v)
    {
      _comp = {v, v, v};
      return *this;
    }

    //! Add another Tensor1
    Tensor1& operator+=(const Tensor1& b)
    {
      _comp[0] += b._comp[0];
      _comp[1] += b._comp[1];
      _comp[2] += b._comp[2];
      return *this;
    }

    //! Subtract another Tensor1
    Tensor1& operator-=(const Tensor1& b)
    {
      _comp[0] -= b._comp[0];
      _comp[1] -= b._comp[1];
      _comp[2] -= b._comp[2];
      return *this;
    }

    //! Scale with scalar
    Tensor1& operator*=(double s)
    {
      _comp[0] *= s;
      _comp[1] *= s;
      _comp[2] *= s;
      return *this;
    }

    //! Divide by scalar
    Tensor1& operator/=(double s)
    {
      return(this->operator*=(1.0/s));
    }

    //! Get the norm of the tensor
    double norm(void) const;

    //! Transfer to a libMesh::Point
    libMesh::Point get_point(void) const;


  private:

    //! The internal representation
    std::array<double, 3> _comp = {0.0, 0.0, 0.0};

};


#endif // _TENSOR1_H_
