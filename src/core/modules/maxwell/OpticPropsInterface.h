#ifndef OPTIC_PROPS_INTERFACE_H_
#define OPTIC_PROPS_INTERFACE_H_

/*
 * ExcitonLayer.h
 *
 *  Created on: Sep 26, 2011
 *      Author: paveryan
 */
#include "TypeDefs.h"
#include "PhysicalModel.h"
#include "PhysicalModel.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "InitFailedException.h"
#include <tensor_value.h>

class TBDLEXPORT OpticPropsInterface: public PhysicalModel

{
  public:
     virtual ~OpticPropsInterface(void);

     virtual libMesh::Complex get_dielectric_constant() const {
       return libMesh::Complex(1, 0);
     }

     virtual const TensorValue<Complex> get_optical_epsilon() const {
       return TensorValue<Complex>(1);
     }

     virtual double get_permeability_constant() const {
       return 1.0;
     }

     virtual double get_spml() const {
       return -1.0;
     }

     static OpticPropsInterface* create(const std::string& name,
         const PhysicalObject* owner, const ModelOptions& options = ModelOptions());


   protected:

     //! \copydoc DriftDiffusionProperties::DriftDiffusionProperties()
     OpticPropsInterface(const ModelOptions& options);

};

inline
OpticPropsInterface::OpticPropsInterface(const ModelOptions& options)
 : PhysicalModel(options)
{
}

inline
OpticPropsInterface::~OpticPropsInterface(void)
{
}


inline
OpticPropsInterface*
OpticPropsInterface::create(const std::string& name,
    const PhysicalObject* owner, const ModelOptions& options)
{
  return dynamic_cast<OpticPropsInterface*>(
      PhysicalModel::create("opticmodel_" + name, owner, options));
}
#endif
