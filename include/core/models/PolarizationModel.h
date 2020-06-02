// $Id$

#ifndef _POLARIZATIONMODEL_H_
#define _POLARIZATIONMODEL_H_

#include "PhysicalModel.h"
#include "SolutionProvider.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "tensor.h"
#include "TensorOperators.h"
#include "libMeshDefs.h"

#include "vector_value.h"
#include "tensor_value.h"

USELIBMESHTYPE(RealVectorValue);


// Base class for charge density models
class PolarizationModel : public PhysicalModel
{

  public:

   virtual ~PolarizationModel(void) {};

   //! Creator for the most basic model
   static PolarizationModel* create(const ModelOptions& options);
  
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

   //! Calculate the polarization
   /*!
    * If a constant value or module name is given in the input file,
    * this will override any internally calculated polarization.
    */
   void calculate(const libMesh::Elem* elem, const libMesh::Point& point);
  

  protected:
  
   PolarizationModel(const ModelOptions& options);
  
   //! Calculate the polarization in model specific way
   virtual void do_calculate(const libMesh::Elem* elem, const libMesh::Point& point);

   //! Initialize the bas class model
   /*!
    * You \c must call this method from derived classes if you want to use
    * default behaviour.
    */
   virtual void do_init(void);

   //! Print some info
   virtual void do_print_info(void);

   //! Set the polarization vector
   void set_polarization(const libMesh::RealVectorValue& polarization);
   
   //! Rotate the polarization vector to calculation system
   void rotate(void);

   Tensor2Sym& get_strain(void);


  private:

   libMesh::RealVectorValue _polarization;

   Tensor2Sym _strain;

   //! We may take it from some other module
   SolutionProvider _polarization_source;

   //! \c true if polarization is given from input or from other module
   bool _fixed_or_external;

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
const RealVectorValue&
PolarizationModel::get_polarization() const
{
  return _polarization;
}


inline
void
PolarizationModel::set_polarization(const libMesh::RealVectorValue& polarization)
{
  if (!_fixed_or_external)
   _polarization = polarization;
}

/*
inline
RealVectorValue&
PolarizationModel::polarization(void)
{
  return _polarization;
}
*/

inline
void 
PolarizationModel::rotate(void)
{

  if (!_fixed_or_external)
  {
    const Material* mat = get_material();
    const RotatedCrystal& cr = mat->get_rotated_crystal();

    _polarization = cr.RotMatrix * _polarization;
  }

}


#endif // _POLARIZATIONMODEL_H_
