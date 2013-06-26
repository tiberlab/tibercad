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
    /*!
     * \param p0 lower left corner
     * \param p1 upper right corner
     * \param nk number of elements along x
     * \param nl number of elements along y
     * \param nm number of elements along z
     */
    TensorGrid(const Point& p0, const Point& p1, int nk = 50, int nl = 50, int nm = 50);


    //! Setup a tensor grid
    /*!
     * \param p0 lower left corner
     * \param p1 upper right corner
     * \param nk number of elements along x
     * \param nl number of elements along y
     * \param nm number of elements along z
     */
    void setup(const Point& p0, const Point& p1, int nk = 50, int nl = 50, int nm = 50);


    //! Find the element for a given point
    /*!
     * The elements are numbered as m*nl*nk + l*nk + k.
     * If the given point is outside the bounding box, -1 is returned.
     */
    int find_element(const Point& p) const;
    

    //! Find the element for a given point
    /*!
     * The element is identified by the triple (k,l,m).
     */
    void find_element(const Point& p, int indices[3]) const;


    //! Get the number of elements
    int num_elements(void) const;


    //! Get the element for given coordinate index triple (k,l,m)
    int index_to_element(int indices[3]) const;


    //! Get the element for given coordinate index triple (k,l,m)
    int index_to_element(unsigned int k, unsigned int l, unsigned int m) const;


  private:

    //! The mesh dimension
    unsigned int _dimension;

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
TensorGrid::find_element(const Point& p, int indices[3]) const
{
  indices[0] = floor((p(0) - _p0(0)) / _dx);
  indices[1] = (_dimension > 1) ? floor((p(1) - _p0(1)) / _dy) : 0;
  indices[2] = (_dimension > 2) ? floor((p(2) - _p0(2)) / _dz) : 0;
}



inline
int
TensorGrid::find_element(const Point& p) const
{
  int k = floor((p(0) - _p0(0)) / _dx);
  int l = (_dimension > 1) ? floor((p(1) - _p0(1)) / _dy) : 0;
  int m = (_dimension > 2) ? floor((p(2) - _p0(2)) / _dz) : 0;

  if ((k < 0) || (l < 0) || (m < 0))
    return -1;

  if ((k >= _nk) || (l >= _nl) || (m >= _nm))
    return -1;

  return(m*_nl*_nk + l*_nk + k);
}


inline
int
TensorGrid::num_elements(void) const
{
  return _nk * _nl * _nm ;
}


inline
int
TensorGrid::index_to_element(int indices[3]) const
{
  return(indices[2]*_nl*_nk + indices[1]*_nk + indices[0]);
}

inline
int
TensorGrid::index_to_element(unsigned int k, unsigned int l, unsigned int m) const
{
  return(m*_nl*_nk + l*_nk + k);
}

#endif // _TENSORGRID_H_
