
/*
 * ExcitonLayer.h
 *
 *  Created on: Sep 26, 2011
 *      Author: paveryan
*/

#ifndef OPTIC_PROPS_W_DEP_H_
#define OPTIC_PROPS_W_DEP_H_

#include "TypeDefs.h"
#include "PhysicalModelInterface.h"
#include "PhysicalModel.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "OpticPropsInterface.h"
#include "Database.h"


class OpticPropsWDependant: public OpticPropsInterface

{
  public:
    OpticPropsWDependant(const ModelOptions& options);

     static OpticPropsWDependant* create(const ModelOptions& options);

     PhysicalModelInterface* create_new() const;

     virtual Complex get_dielectric_constant() const;

     virtual const TensorValue<Complex> get_optical_epsilon() const {
       Complex t = get_dielectric_constant();
       return TensorValue<Complex>(t, 0, 0, 0, t, 0, 0, 0, t);
     }

     virtual double get_permeability_constant() const {
       return 1;
     }

     virtual double get_spml() const {
       return -1;
     }

     virtual void do_reinit(void);
     //virtual Complex get_dielectric_constant(double param) const;

  protected:
     std::vector<double> params;
     std::vector<Complex> epsilon;

     int currentIndex;
     Complex currentEpsilon;

     std::string dataFile;

     virtual void read_database();

     virtual void do_init(void);
};

inline OpticPropsWDependant* OpticPropsWDependant::create(const ModelOptions& options)
{
  return new OpticPropsWDependant(options);
}

inline
OpticPropsWDependant::OpticPropsWDependant(const ModelOptions& options) : OpticPropsInterface(options) {
}

inline
PhysicalModelInterface*
OpticPropsWDependant::create_new(void) const
{
  return new OpticPropsWDependant(get_options());
}

#endif

