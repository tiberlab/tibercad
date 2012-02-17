#ifndef __ISCALAR_FEBASE_H__
#define __ISCALAR_FEBASE_H__

#include "IFEBase.h"
#include "ScalarFunction.h"

class IScalarFEBase : public IFEBase<ScalarFunction> {
  public:
    IScalarFEBase(const unsigned int dim, double scaling) : IFEBase<ScalarFunction>(dim, scaling){
    }

    //TODO symmetry
    void applyScaling() {
      double x0 = length_scaling;

      for (unsigned int j = 0; j < functions.size(); j++) {
        for (unsigned int i = 0; i < functions[j].phi.size(); i++) {
          functions[j].grad[i] *= x0;
        }
      }
    }
};

#endif
