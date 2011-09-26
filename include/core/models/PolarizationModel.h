// $Id$

#ifndef _POLARIZATIONMODEL_H_
#define _POLARIZATIONMODEL_H_

#include "PhysicalModelInterface.h"
#include "vector_value.h"
#include "tensor_value.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "tensor.h"
#include "TensorOperators.h"

class Elem;
class Point;

// Base class for charge density models
class PolarizationModel : public PhysicalModelInterface
{

  public:

   virtual ~PolarizationModel(void) {};
  
   //! Get the polarization vector
   /*!
    * \return Polarization in calculation coordinate system
    */
   const RealVectorValue& get_polarization(void) const;

   //! Set the strain
   /*!
    * Strain tensor has to be provided in crystal coordinate system
    */
   void set_strain(const Tensor2Sym& strain);
  
   virtual void calculate(const Elem* elem, const Point& point) = 0;
  

  protected:
  
   PolarizationModel(const ModelOptions& options);

   void set_polarization(const RealVectorValue& polarization);
   
   RealVectorValue& polarization(void);

   void rotate(void);

   Tensor2Sym& get_strain(void);


  private:

   RealVectorValue _polarization;

   Tensor2Sym _strain;

};


inline
void 
PolarizationModel::set_strain(const Tensor2Sym& strain)
{
  _strain = strain;

}


inline
Tensor2Sym&
PolarizationModel::get_strain(void)
{
  return _strain;
}


inline
PolarizationModel::PolarizationModel(const ModelOptions& options) :
  PhysicalModelInterface(options)
{
}


inline
const RealVectorValue&
PolarizationModel::get_polarization() const
{
  return _polarization;
}


inline
void
PolarizationModel::set_polarization(const RealVectorValue& polarization)
{
   _polarization = polarization;
}


inline
RealVectorValue&
PolarizationModel::polarization(void)
{
  return _polarization;
}


inline
void 
PolarizationModel::rotate(void)
{

  const Material* mat = get_material();
  const RotatedCrystal&   cr = mat->get_rotated_crystal();

  _polarization = cr.RotMatrix * _polarization;

}


#endif // _POLARIZATIONMODEL_H_
