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
 * \file StiffnessModel.h
 * \brief tiberCAD elasticity module header.
 *
 * \note This file is part of module elasticity.
 */


#ifndef TC_STIFFNESSMODEL_H
#define TC_STIFFNESSMODEL_H

#include "tibercad/physics/PhysicalModel.h"

#include "tibercad/physics/Material.h"
#include "tibercad/math/TensorOperators.h"

#include "libmesh/elem.h"
#include "libmesh/tensor_value.h"


using namespace std;

class StiffnessModel : public PhysicalModel
{

  public:

  //! Destructor
  ~StiffnessModel(void) {};
  
  //! Creator function
  static StiffnessModel* create(const ModelOptions& options);

  Tensor4DSym get_stiffness(void) const;

  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point){};

  protected:

    //! Constructor
  StiffnessModel(const ModelOptions& options);


  void set_stiffness_constant(const Tensor4DSym& C);

  void rotate(void);

  private:

  Tensor4DSym _stiffness;

};

inline
Tensor4DSym
StiffnessModel::get_stiffness(void) const
{
  return _stiffness;
}


inline 
void 
StiffnessModel::set_stiffness_constant(const Tensor4DSym& C)
{
 
  _stiffness = C;

}



inline
StiffnessModel::StiffnessModel(const ModelOptions& options) :
PhysicalModel(options)
{
}


inline 
void 
StiffnessModel::rotate()
{

  const Material* mat = get_material();
  const libMesh::RealTensor& rot = mat->get_rotation_matrix();
  Tensor2 rotm;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      rotm(i+1, j+1) = rot(i, j);

  _stiffness = push_forward(_stiffness, rotm);
}

#endif // TC_THERMALCONDUCTIVITYMODEL_H
