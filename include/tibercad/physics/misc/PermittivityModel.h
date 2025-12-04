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
 * \file PermittivityModel.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _PERMITTIVITYMODEL_H_
#define _PERMITTIVITYMODEL_H_

#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/physics/Material.h"
#include "tibercad/base/libMeshDefs.h"

#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

USELIBMESHTYPE(RealTensor);


// Base class for charge density models
class TBDLEXPORT PermittivityModel : public PhysicalModel
{

  public:

  virtual ~PermittivityModel(void) {};
  
  const RealTensor& get_permittivity(void) const;
  
  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point) = 0;
  
protected:
  
    PermittivityModel(const ModelOptions& options);

    void set_permittivity(const libMesh::RealVectorValue& permittivity);

  // void set_permittivity(const std::vector<double>& permittivity);
   
    void rotate(void);


  private:

   RealTensor _permittivity;

};


inline
PermittivityModel::PermittivityModel(const ModelOptions& options) :
  PhysicalModel(options)
{
  _permittivity = 0;
  _permittivity(0,0) = 1.0;
  _permittivity(1,1) = 1.0;
  _permittivity(2,2) = 1.0;
}

inline
const RealTensor&
PermittivityModel::get_permittivity() const
{
 
  return _permittivity;
}

inline
void
PermittivityModel::set_permittivity(const libMesh::RealVectorValue& permittivity_diag)
{
  
  _permittivity(0,0) =  permittivity_diag(0);
  _permittivity(1,1) =  permittivity_diag(1);
  _permittivity(2,2) =  permittivity_diag(2);
  
}

// inline
// void
// PermittivityModel::set_permittivity(const std::vector<double>& permittivity)
// {
//   _permittivity = 0;
//   _permittivity(0,0) =  permittivity[0];
//   _permittivity(1,1) =  permittivity[1];
//   _permittivity(2,2) =  permittivity[2];
  
// }

inline
void 
PermittivityModel::rotate(void)
{

  const Material* mat = get_material();
  const libMesh::RealTensor& rotm = mat->get_rotation_matrix();

  _permittivity = rotm * (_permittivity * rotm.transpose());

}

#endif // _POLARIZATIONMODEL_H_
