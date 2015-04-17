// $Id$

#include "TensorGrid.h"



void
TensorGrid::setup(const Point& p0, const Point& p1, int nk, int nl, int nm)
{
  _dimension = 3;
  _p0 = p0;
  _p1 = p1;
  _nk = nk;
  _nl = nl;
  _nm = nm;

  double lx = _p1(0) - _p0(0);
  double ly = _p1(1) - _p0(1);
  double lz = _p1(2) - _p0(2);
  if (lz == 0)
    _dimension = 2;
  if (ly == 0)
    _dimension = 1;

  double eps = 1e-6;
  _p0(0) -= eps * (lx + 1);
  _p0(1) -= eps * (ly + 1);
  _p0(2) -= eps * (lz + 1);
  _p1(0) += eps * (lx + 1);
  _p1(1) += eps * (ly + 1);
  _p1(2) += eps * (lz + 1);

  _dx = (_p1(0) - _p0(0)) / _nk;
  _dy = (_p1(1) - _p0(1)) / _nl;
  _dz = (_p1(2) - _p0(2)) / _nm;
}


Point
TensorGrid::get_centroid(unsigned int i) const
{
  unsigned int k, l, m;
  element_to_index(i, k, l, m);

  Point centroid(_p0);
  centroid += Point((k+0.5)*_dx, (l+0.5)*_dy, (m+0.5)*_dz);

  return(centroid);
}
