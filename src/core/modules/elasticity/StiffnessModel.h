// $Id$

#ifndef _STIFFNESSMODEL_H_
#define _STIFFNESSMODEL_H_

#include "PhysicalModel.h"

#include "libmesh/tensor_value.h"
#include "tensor.h"
#include "libmesh/elem.h"
#include "Material.h"
#include "TensorOperators.h"



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
  Tensor2Gen rotm;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      rotm(i+1, j+1) = rot(i, j);

  _stiffness = push_forward(_stiffness, rotm);
}

#endif // _THERMALCONDUCTIVITYMODEL_H_
