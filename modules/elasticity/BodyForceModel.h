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
 * \file BodyForceModel.h
 * \brief tiberCAD elasticity module header.
 *
 * \note This file is part of module elasticity.
 */


#ifndef _BODYFORCEMODEL_H_
#define _BODYFORCEMODEL_H_

#include "tibercad/physics/PhysicalModel.h"
//#include "tibercad/base/libMeshDefs.h"

#include "tensor_value.h"
#include "vector_value.h"
#include "elem.h"

//class Elem;
//class Point;

using namespace std;

//! The base class for body force models
class BodyForceModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~BodyForceModel(void) {};

    //! Creator function
    static BodyForceModel* create(const ModelOptions& options);

    const libMesh::RealGradient& get_force_source(void) const;

    const libMesh::RealTensor& get_stress_source(void) const;

    const libMesh::RealTensor& get_strain_source(void) const;


    //! Calculate local body force
    /*!
     * \param elem pointer to the current element
     * \param point the coordinates in the reference element
     */
    virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point){};

  protected:

    //! Constructor
    BodyForceModel(const ModelOptions& options);

    void set_force_source(const libMesh::RealGradient& force_source);

    void set_strain_source(const libMesh::RealTensor& strain_source);

    void set_stress_source(const libMesh::RealTensor& stress_source);



  private:

    libMesh::RealGradient _force_source;

    libMesh::RealTensor _strain_source;

    libMesh::RealTensor _stress_source;

};

inline
const libMesh::RealGradient&
BodyForceModel::get_force_source(void) const
{
  return _force_source;
}

inline
const libMesh::RealTensor&
BodyForceModel::get_strain_source(void) const
{
  return _strain_source;
}

inline
const libMesh::RealTensor&
BodyForceModel::get_stress_source(void) const
{
  return _stress_source;
}

inline 
void 
BodyForceModel::set_force_source(const libMesh::RealGradient& force_source)
{
  _force_source = force_source;
}

inline 
void 
BodyForceModel::set_strain_source(const libMesh::RealTensor& strain_source)
{
  _strain_source = strain_source;
}


inline 
void 
BodyForceModel::set_stress_source(const libMesh::RealTensor& stress_source)
{
  _stress_source = stress_source;
}


inline
BodyForceModel::BodyForceModel(const ModelOptions& options) :
  PhysicalModel(options),
  _force_source(0),
  _strain_source(0),
  _stress_source(0)
{
}

#endif // _THERMALCONDUCTIVITYMODEL_H_
