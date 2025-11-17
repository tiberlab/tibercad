// $Id$

#ifndef _POLARIZATIONMODEL_H_
#define _POLARIZATIONMODEL_H_

#include "PhysicalModel.h"
#include "SolutionProvider.h"
#include "Material.h"
#include "Tensor2.h"
#include "libMeshDefs.h"

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

USELIBMESHTYPE(RealVectorValue);


//! Base class for electrical polarization models
/*!
 * Polarization is assumed in units of C/m^2
 *
 * A fixed polarisation vector can be specified as
 * \c polarization = (Px, Py, Pz)
 * In this case, it assumed in the simulation coordinate system.
 */
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
   void set_strain(const Tensor2& strain);

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

   Tensor2& get_strain(void);

   //! Get strain tensor
   void get_strain(libMesh::RealTensor& strain) const;

   //! Get strain as vector, in Voigt order
   void get_strain(std::vector<double>& strain) const;


  private:

   //! The polarization in C/m^2
   libMesh::RealVectorValue _polarization;

   //! Polarization might typically depend on strain
   Tensor2 _strain;

   //! We may take it from some other module
   SolutionProvider _polarization_source;

   //! \c true if polarization is given from input or from other module
   bool _fixed_or_external;

};


inline
void 
PolarizationModel::set_strain(const Tensor2& strain)
{
  _strain = strain;

}


inline
Tensor2&
PolarizationModel::get_strain(void)
{
  return _strain;
}


inline
void
PolarizationModel::get_strain(libMesh::RealTensor& strain) const
{
  strain(0,0) = _strain(1,1);
  strain(1,1) = _strain(2,2);
  strain(2,2) = _strain(3,3);
  strain(1,2) = strain(2,1) = _strain(3,2);
  strain(0,2) = strain(2,0) = _strain(3,1);
  strain(0,1) = strain(1,0) = _strain(2,1);
}

inline
void
PolarizationModel::get_strain(std::vector<double>& strain) const
{
  strain.resize(6, 0.0);
  strain[0] =_strain(1,1); 
  strain[1] =_strain(2,2); 
  strain[2] =_strain(3,3); 
  strain[3] =_strain(3,2); 
  strain[4] =_strain(3,1); 
  strain[5] =_strain(2,1); 
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






#endif // _POLARIZATIONMODEL_H_
