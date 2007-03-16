// $Id$


#ifndef _FINITEELEMENT_H_
#define _FINITEELEMENT_H_


#include "TiberCad.h"


#include "fe.h"


//! A finite element extension which adds scaling and symmetry
/*!
 * This extension to the libMesh finite element families includes
 * length scaling and symmetries in the gradients of the shape functions
 * and the Jacobi determinants.
 *
 * Explain how ...
 */
template <unsigned int Dim, FEFamily T>
class FiniteElement : public FE<Dim, T>
{

  public:

    FiniteElement(const FEType& type);
    
    virtual void reinit(const Elem* elem,
        const std::vector<Point>* points = NULL);

    virtual void reinit(const Elem* elem, const unsigned int side);
    
    //! Set the symmetry
    void set_symmetry(TiberCad::Symmetry symmetry);

    //! Set the length scaling and the mesh units
    void set_scaling(double length_scaling, double mesh_units);


  private:

    //! The symmetry
    TiberCad::Symmetry _symmetry;

    //! The length scaling
    double _length_scaling;

    //! The mesh units
    double _mesh_units;
    
};



//
// inline methods
// 

template <unsigned int Dim, FEFamily T>
inline
FiniteElement<Dim, T>::FiniteElement(const FEType& type)
  : FE<Dim, T>(type),
    _symmetry(TiberCad::NONE),
    _length_scaling(1.0),
    _mesh_units(1.0)
{
}


template <unsigned int Dim, FEFamily T>
inline
void
FiniteElement<Dim, T>::set_symmetry(TiberCad::Symmetry symmetry)
{
  _symmetry = symmetry;
}


template <unsigned int Dim, FEFamily T>
inline
void
FiniteElement<Dim, T>::set_scaling(double length_scaling, double mesh_units)
{
  _length_scaling = length_scaling;
  _mesh_units = mesh_units;
}


#endif // _FINITEELEMENT_H_
