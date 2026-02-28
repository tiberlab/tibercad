/*  
 * This file is part of the tiberCAD module elasticity.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file ElasticityBoundaryModel.h
 * \brief tiberCAD elasticity module header.
 *
 * \note This file is part of module elasticity.
 */


#ifndef _ELASTICITYBOUNDARYMODEL_H_tens_
#define _ELASTICITYBOUNDARYMODEL_H_tens_

#include "tibercad/physics/PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/base/libMeshDefs.h"

//class Elem;
//class Point;

//! This is the base class for the Poisson physical model
class ElasticityBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~ElasticityBoundaryModel(void) {};
  
    //! Calculate properties 
    virtual void calculate(const libMesh::Elem* elem, unsigned int side, const libMesh::Point& point) = 0;

    //! Creator function
    static ElasticityBoundaryModel* create(const MaterialBoundary* boundary,
        const ModelOptions& options);

  const bool is_extended(void);

  const void get_coefficients(libMesh::RealTensor& H, double& A,libMesh::RealGradient& R);

    void set_normal(const Point p);


  protected:

    //! Constructor
    ElasticityBoundaryModel(const ModelOptions& options);

  void set_coefficients(libMesh::RealTensor H,double A, libMesh::RealGradient R);
  
  void set_is_extended(bool is_extended);

  Point _normal; 

  private:


 

  //! constrain matrix
  libMesh::RealTensor _H_tens;
 
  //! is extended
  double _coeff;

  //! contrain vector
  libMesh::RealGradient _R_vec;

  //! is extended
  bool _is_extended;

  static TiberModelObject*  _create(const ModelOptions& options);
  
  static void  _destroy( TiberModelObject* p);
 
};

inline
const
void
ElasticityBoundaryModel::get_coefficients(libMesh::RealTensor& H, double& A, libMesh::RealGradient& R)
{
  H = _H_tens;
  R = _R_vec;
  A = _coeff;
}

inline
void
ElasticityBoundaryModel::set_normal(const Point normal)
{
  _normal = normal;
}


inline
void
ElasticityBoundaryModel::set_is_extended(bool is_extended)
{
  _is_extended = is_extended;
}

inline
const
bool 
ElasticityBoundaryModel::is_extended(void)
{
  return _is_extended;
}

inline
void
ElasticityBoundaryModel::set_coefficients(libMesh::RealTensor H, double A, libMesh::RealGradient R)
{
  _H_tens = H;
  _R_vec = R;
  _coeff = A;
}

inline
ElasticityBoundaryModel::ElasticityBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options),
  _H_tens(0),
  _R_vec(0),
  _is_extended(0)
{
}




#endif // _MYPOISSONMODEL_H_tens_
