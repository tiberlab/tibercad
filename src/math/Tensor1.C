#include "tibercad/math/Tensor1.h"

#include "libmesh/point.h"

#include <cmath>

Tensor1::Tensor1(const libMesh::Point& p)
{
  _comp[0] = p(0);
  _comp[1] = p(1);
  _comp[2] = p(2);
}

libMesh::Point
Tensor1::get_point(void) const
{
  return(libMesh::Point(_comp[0], _comp[1], _comp[2]));
}

    
double
Tensor1::norm(void) const
{
  double d = _comp[0]*_comp[0] + _comp[1]*_comp[1] + _comp[2]*_comp[2];

  return std::sqrt(d);
}

