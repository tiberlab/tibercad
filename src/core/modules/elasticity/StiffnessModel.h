// $Id$

#ifndef _STIFFNESSMODEL_H_
#define _STIFFNESSMODEL_H_

#include "PhysicalModel.h"

#include "tensor_value.h"
#include "tensor.h"
#include "elem.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "TensorOperators.h"


//class Elem;
//class Point;

using namespace std;

//! The base class for Poisson boundary conditions
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

//! Calculate for a point on the given side
  //void calculate(const Elem* elem, const Point& point);

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
  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

  _stiffness = push_forward(_stiffness, cr.RotMatrix);
}

#endif // _THERMALCONDUCTIVITYMODEL_H_
