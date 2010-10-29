// $Id$

#ifndef _BOUNDARYMODEL_H_
#define _BOUNDARYMODEL_H_

#include "PhysicalModelInterface.h"

#include "tensor_value.h"
#include "tensor.h"
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

  const void get_coefficients(RealTensor& H, double& A, RealGradient& R);


  protected:

    //! Constructor
  BoundaryModel(const ModelOptions& options);

  void set_coefficients(RealTensor H, double A, RealGradient R);

  private:

  RealTensor _H;

  RealGradient _R;

  double _A;


};


inline
void
BoundaryModel::set_coefficients(RealTensor H, double A, RealGradient R)
{
  _H = H;
  _R = R;
  _A = A;
}


inline
const
void
BoundaryModel::get_coefficients(RealTensor& H, double& A, RealGradient& R)
{
  H = _H;
  R = _R;
  A = _A;
}

inline
BoundaryModel::BoundaryModel(const ModelOptions& options) :
  PhysicalModelInterface(options),
  _H(0),
  _R(0),
  _A(0)
{
}


#endif // _THERMALCONDUCTIVITYMODEL_H_
