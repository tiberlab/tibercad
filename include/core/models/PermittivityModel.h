// $Id$

#ifndef _PERMITTIVITYMODEL_H_
#define _PERMITTIVITYMODEL_H_

#include "PhysicalModelInterface.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "tensor.h"
#include "TensorOperators.h"

class Elem;
class Point;

// Base class for charge density models
class TBDLEXPORT PermittivityModel : public PhysicalModelInterface
{

  public:

  virtual ~PermittivityModel(void) {};
  
  const RealTensor& get_permittivity(void) const;
  
  virtual void calculate(const Elem* elem, const Point& point) = 0;
  
protected:
  
    PermittivityModel(const ModelOptions& options);

    void set_permittivity(const RealVectorValue& permittivity);

  // void set_permittivity(const std::vector<double>& permittivity);
   
    void rotate(void);


  private:

   RealTensor _permittivity;

};


inline
PermittivityModel::PermittivityModel(const ModelOptions& options) :
  PhysicalModelInterface(options)
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
PermittivityModel::set_permittivity(const RealVectorValue& permittivity_diag)
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

  Material* mat = get_material();
  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

 Tensor2Gen rotate = cr.RotMatrix;

 _permittivity = rotate * (_permittivity * rotate.transpose());


}

#endif // _POLARIZATIONMODEL_H_
