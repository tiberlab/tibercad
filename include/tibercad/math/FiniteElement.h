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
 * \file FiniteElement.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_FINITEELEMENT_H
#define TC_FINITEELEMENT_H


#include "tibercad/base/TiberCad.h"
#include "tibercad/base/tiber_dll.h"


#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"




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
        const std::vector<libMesh::Point>* const points = nullptr,
        const std::vector<double>* const weights = nullptr);



    virtual void reinit(const libMesh::Elem* elem, const unsigned int side,
        const double tolerance = libMesh::TOLERANCE,
        const std::vector<libMesh::Point>* const points = nullptr,
        const std::vector<double>* const weights = nullptr);
    
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


#endif // TC_FINITEELEMENT_H
