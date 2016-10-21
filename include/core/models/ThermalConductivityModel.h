// $Id: ThermalConductivityModel.h 2132 2010-10-29 08:19:25Z maufder $

#ifndef _THERMALCONDUCTIVITYMODEL_H_
#define _THERMALCONDUCTIVITYMODEL_H_

#include "PhysicalModel.h"
#include "TensorOperators.h"

#include "tensor_value.h"

#include "vector_value.h"
#include "RotatedCrystal.h"

#include "libMeshDefs.h"


//! The base class for Poisson boundary conditions
class ThermalConductivityModel : public PhysicalModelInterface
{

  public:

    //! Destructor
  virtual ~ThermalConductivityModel(void) {};
  
  //! Creator function
  static ThermalConductivityModel* create(const ModelOptions& options);

  virtual void calculate(const Elem* elem, const Point& point) = 0;

  const libMesh::RealTensor& get_thermal_conductivity(void);

  protected:

    //! Constructor
  ThermalConductivityModel(const ModelOptions& options);

  public:

  // void set_thermal_conductivity(RealTensor kappa);

  void set_thermal_conductivity(libMesh::RealGradient kappa);

  protected:

  void rotate(void);

  libMesh::RealTensor _kappa;

  private:

};

inline
const libMesh::RealTensor&
ThermalConductivityModel::get_thermal_conductivity(void) 
{
  return _kappa;
}

inline
ThermalConductivityModel::ThermalConductivityModel(const ModelOptions& options)
 : PhysicalModelInterface(options)
{
}



inline 
void 
ThermalConductivityModel::set_thermal_conductivity(libMesh::RealGradient kappa)
{

  _kappa = 0;
  _kappa(0,0) = kappa(0);
  _kappa(1,1) = kappa(1);
  _kappa(2,2) = kappa(2);

}


inline
void
ThermalConductivityModel::rotate(void)
{

  if (get_material()->get_structure() == "wz")
  {
    const RotatedCrystal&   cr = get_material()->get_rotated_crystal();
    Tensor2Gen RotMatrix = cr.RotMatrix; 
    _kappa = RotMatrix * ( _kappa * (RotMatrix.transpose()));
  }

}





// inline
// void
// ThermalConductivityModel::rotate_to_calculation_system(const Tensor2Gen& RotMatrix)
// {

//   Tensor2Gen kappa(0);
//   for (ID i = 0; i<3; i ++)
//     for (ID j = i; j<3; j ++)
//       kappa(i+1,j+1) = _kappa(i,j);


//   // generates conductivity matrix in calculation system
//   kappa = sym(RotMatrix * ( kappa * (RotMatrix.transpose())));
  
//   _kappa = 0.0;
//   for (ID i = 0; i<3; i ++)
//     for (ID j = i; j<3; j ++)
//       _kappa(i,j) = kappa(i+1,j+1);

// }

#endif // _THERMALCONDUCTIVITYMODEL_H_
