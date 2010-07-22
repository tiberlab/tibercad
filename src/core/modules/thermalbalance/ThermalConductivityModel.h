// $Id$

#ifndef _THERMALCONDUCTIVITYMODEL_H_
#define _THERMALCONDUCTIVITYMODEL_H_

#include "PhysicalModel.h"

#include "tensor_value.h"

class Elem;
class Point;

using namespace std;

//! The base class for Poisson boundary conditions
class ThermalConductivityModel : public PhysicalModelInterface
{

  public:

    //! Destructor
  virtual ~ThermalConductivityModel(void) {};
  
  //! Creator function
  static ThermalConductivityModel* create(const ModelOptions& options);

  virtual void calculate(const Elem* elem, const Point& point) = 0;

  RealTensor& get_thermal_conductivity(void);

  protected:

    //! Constructor
  ThermalConductivityModel(const ModelOptions& options);

  public:
  void set_thermal_conductivity(RealTensor kappa);
  protected:

  void rotate_to_calculation_system(const Tensor2Gen& RotMatrix);

  RealTensor _kappa;

  private:

};

inline
RealTensor&
ThermalConductivityModel::get_thermal_conductivity(void) 
{
  return _kappa;
}

ThermalConductivityModel::ThermalConductivityModel(const ModelOptions& options) :
PhysicalModelInterface(options)
{
}


inline 
void 
ThermalConductivityModel::set_thermal_conductivity(RealTensor kappa)
{
  _kappa = kappa;
}



inline
void
ThermalConductivityModel::rotate_to_calculation_system(const Tensor2Gen& RotMatrix)
{

  Tensor2Gen kappa(0);
  for (ID i = 0; i<3; i ++)
    for (ID j = i; j<3; j ++)
      kappa(i+1,j+1) = _kappa(i,j);


  // generates conductivity matrix in calculation system
  kappa = sym(RotMatrix * ( kappa * (RotMatrix.transpose())));
  
  _kappa = 0.0;
  for (ID i = 0; i<3; i ++)
    for (ID j = i; j<3; j ++)
      _kappa(i,j) = kappa(i+1,j+1);

}

#endif // _THERMALCONDUCTIVITYMODEL_H_
