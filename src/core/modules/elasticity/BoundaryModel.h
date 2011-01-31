// $Id$

#ifndef _BOUNDARYMODEL_H_tensor_
#define _BOUNDARYMODEL_H_tensor_

#include "PhysicalModelInterface.h"

#include "tensor_value.h"
#include "tensor.h"
#include "point.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "TensorOperators.h"


class Elem;
class Point;

using namespace std;

//! The base class for Poisson boundary conditions
class BoundaryModel : public PhysicalModelInterface
{

  public:

  //! Destructor
  ~BoundaryModel(void) {};
  
  //! Creator function
  static BoundaryModel* create(const ModelOptions& options);


   virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point) = 0;

  void get_coefficients(RealTensor& H, RealGradient& R);


 void set_normal(const Point p);

  protected:

    //! Constructor
  BoundaryModel(const ModelOptions& options);

  Point _normal;

  void set_coefficients(RealTensor H, RealGradient R);

  private:

  RealTensor _H_tensor;

  RealGradient _R_vec;

};


inline
void
BoundaryModel::set_coefficients(RealTensor H, RealGradient R)
{
  _H_tensor = H;
  _R_vec = R;
}


inline
void
BoundaryModel::get_coefficients(RealTensor& H, RealGradient& R)
{
  H = _H_tensor;
  R = _R_vec;
}

inline
BoundaryModel::BoundaryModel(const ModelOptions& options) :
  PhysicalModelInterface(options),
  _H_tensor(0),
  _R_vec(0)
{
}

inline
void
BoundaryModel::set_normal(const Point normal)
{
  _normal = normal;
}



#endif // _THERMALCONDUCTIVITYMODEL_H_tensor_
