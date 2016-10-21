// $Id$


#ifndef _FINITEELEMENT_H_
#define _FINITEELEMENT_H_


#include "TiberCad.h"
#include "tiber_dll.h"
//#include "libMeshDefs.h"


#include "fe.h"
#include "elem.h"
#include "point.h"



//USELIBMESHTYPE(FEType);
//USELIBMESHTYPE(Point);
//USELIBMESHTYPE(Elem);


//! A finite element extension which adds scaling and symmetry
/*!
 * This extension to the libMesh finite element families includes
 * length scaling and symmetries in the gradients of the shape functions
 * and the Jacobi determinants.
 *
 * Explain how ...
 */
template <unsigned int Dim, libMesh::FEFamily T>
class FiniteElement : public libMesh::FE<Dim, T>
{

  public:

    FiniteElement(const libMesh::FEType& type);
    
    virtual void reinit(const libMesh::Elem* elem,
        const std::vector<libMesh::Point>* const points = NULL,
        const std::vector<double>* const weights = NULL);



    virtual void reinit(const libMesh::Elem* elem, const unsigned int side,
        const double tolerance = libMesh::TOLERANCE,
        const std::vector<libMesh::Point>* const points = NULL,
        const std::vector<double>* const weights = NULL);
    
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

template <unsigned int Dim, libMesh::FEFamily T>
inline
FiniteElement<Dim, T>::FiniteElement(const libMesh::FEType& type)
  : libMesh::FE<Dim, T>(type),
    _symmetry(TiberCad::NONE),
    _length_scaling(1.0),
    _mesh_units(1.0)
{
}


template <unsigned int Dim, libMesh::FEFamily T>
inline
void
FiniteElement<Dim, T>::set_symmetry(TiberCad::Symmetry symmetry)
{
  _symmetry = symmetry;
}


template <unsigned int Dim, libMesh::FEFamily T>
inline
void
FiniteElement<Dim, T>::set_scaling(double length_scaling, double mesh_units)
{
  _length_scaling = length_scaling;
  _mesh_units = mesh_units;
}


#endif // _FINITEELEMENT_H_
