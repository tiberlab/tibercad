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
 * \file ElasticityModel.h
 * \brief tiberCAD elasticity module header.
 *
 * \note This file is part of module elasticity.
 */


#ifndef TC_ELASTICITYMODEL_H
#define TC_ELASTICITYMODEL_H

#include "tibercad/physics/PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "StiffnessModel.h"
#include "BodyForceModel.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/io/Messages.h"
class StiffnessModel;
//class Elem;


//! This is the base class for the Poisson physical model
class TC_DLLOCAL ElasticityModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~ElasticityModel(void) {};
  
    //! Calculate properties 
    void  calculate(const libMesh::Elem* elem, const libMesh::Point& point);


    //! Creator function
    static ElasticityModel* create(const Material* mat, const ModelOptions& options);
  
    const Tensor4DSym& get_stiffness(void);
  
    const libMesh::RealGradient& get_force_source(void);

    const libMesh::RealTensor& get_stress_source(void);

    const libMesh::RealTensor& get_strain_source(void);
 
  protected:

    //! Constructor
    ElasticityModel(const ModelOptions& options);

    virtual void do_init(void);

    //! Print some useful information
    virtual void do_print_info(void);

    virtual void prepare_submodels(void);

  private:


   static TiberModelObject* _create(const ModelOptions& options);

  static void  _destroy( TiberModelObject* p);

 //!Body force model map
  std::vector<BodyForceModel* > _bfm;


  Tensor4DSym _stiffness;

  libMesh::RealGradient _force;
  
  libMesh::RealTensor _strain;
  
  libMesh::RealTensor _stress;

 
   
};


inline 
const Tensor4DSym&
ElasticityModel::get_stiffness()
{

  return _stiffness;
}

inline 
const libMesh::RealGradient&
ElasticityModel::get_force_source()
{
  return _force;
}

inline 
const libMesh::RealTensor&
ElasticityModel::get_strain_source()
{
  return _strain;
}

inline 
const libMesh::RealTensor&
ElasticityModel::get_stress_source()
{
  return _stress;
}

inline
ElasticityModel::ElasticityModel(const ModelOptions& options) :
  PhysicalModel(options),
  _force(0),
  _stress(0),
  _strain(0),
  _stiffness(0)
{
}

inline
TiberModelObject*  ElasticityModel::_create(const ModelOptions& options)
{

  return new ElasticityModel(options);

}

inline
void  ElasticityModel::_destroy( TiberModelObject* p)
{

  delete p;

}




#endif // TC_MYPOISSONMODEL_H
