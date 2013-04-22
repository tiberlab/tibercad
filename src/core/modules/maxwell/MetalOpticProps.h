
/*
 * ExcitonLayer.h
 *
 *  Created on: Sep 26, 2011
 *      Author: paveryan
*/

#ifndef METAL_OPTIC_PROPS_H_
#define METAL_OPTIC_PROPS_H_

#include "TypeDefs.h"
#include "PhysicalModelInterface.h"
#include "PhysicalModel.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "OpticPropsInterface.h"
#include "Database.h"
#include "MaxwellEquations.h"
#include "OpticPropsModel.h"


class MetalOpticProps: public OpticPropsModel, public ICubic

{
  public:
    MetalOpticProps(const ModelOptions& options);

     static MetalOpticProps* create(const ModelOptions& options);

     PhysicalModelInterface* create_new() const;

     virtual bool addCData(CubicEigenSystem& system);
  protected:
     Complex alpha;

     virtual void read_database();

     virtual void do_init(void);
};

inline MetalOpticProps* MetalOpticProps::create(const ModelOptions& options)
{
  return new MetalOpticProps(options);
}

inline
MetalOpticProps::MetalOpticProps(const ModelOptions& options) : OpticPropsModel(options) {
}

inline
PhysicalModelInterface*
MetalOpticProps::create_new(void) const
{
  return new MetalOpticProps(get_options());
}

#endif

