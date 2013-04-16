// $Id$

#ifndef _TENSORGRID_H_
#define _TENSORGRID_H_


#include "point.h"

//! A class describing a homogeneous tensor grid
/*!
 * We use this class mainly as auxiliary object for fast
 * retrieval of elements in unstructured meshes.
 */
class TensorGrid
{

  public:

    //! Default constructor
    TensorGrid(void) {};

    //! Constructor
    TensorGrid(const Point& p0, const Point& p1, int nk = 50, int nl = 50, int nm = 50);


    //! Setup a tensor grid
    void setup(const Point& p0, const Point& p1, int nk = 50, int nl = 50, int nm = 50);


    //! Find the element for a given point
    /*!
     * The elements are numbered as m*nl*nk + l*nk + k.
     * If the given point is outside the bounding box, -1 is returned.
     */
    int find_element(const Point& p) const;
    

    //! Get the number of elements
    int num_elements(void) const;


  private:

    //! The origin
    Point _p0;

    //! The second point of the bounding box
    Point _p1;

    //! The number of elements in x
    int _nk;

    //! The number of elements in y
    int _nl;

    //! The number of elements in z
    int _nm;

    //! The interval in x
    double _dx;

    //! The interval in y
    double _dy;

    //! The interval in z
    double _dz;

};




inline
TensorGrid::TensorGrid(const Point& p0, const Point& p1, int nk, int nl, int nm)
{
  setup(p0, p1, nk, nl, nm);
}


inline
void
TensorGrid::setup(const Point& p0, const Point& p1, int nk, int nl, int nm)
{
  _p0 = p0;
  _p1 = p1;
  _nk = nk;
  _nl = nl;
  _nm = nm;
  _dx = (_p1(0) - _p0(0)) / _nk;
  _dy = (_p1(1) - _p0(1)) / _nl;
  _dz = (_p1(2) - _p0(2)) / _nm;
}



inline
int
TensorGrid::find_element(const Point& p) const
{
  int k = floor((p(0) - _p0(0)) / _dx);
  int l = floor((p(1) - _p0(1)) / _dy);
  int m = floor((p(2) - _p0(2)) / _dz);

  if ((k < 0) || (l < 0) || (m < 0))
    return -1;

  if ((k > _nk) || (l > _nl) || (m > _nm))
    return -1;

  return m*_nl*_nk + l*_nk + k;
}


inline
int
TensorGrid::num_elements(void) const
{
  return _nk * _nl * _nm;
}


#endif // _TENSORGRID_H_
