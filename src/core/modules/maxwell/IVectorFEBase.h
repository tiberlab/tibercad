#ifndef __IVECTOR_FEBASE_H__
#define __IVECTOR_FEBASE_H__

#include "IFEBase.h"
#include "VectorFunction.h"

class IVectorFEBase : public IFEBase<VectorFunction> {
  public:
    IVectorFEBase(const unsigned int dim, double scaling) : IFEBase<VectorFunction>(dim, scaling){
    }

    //TODO symmetry
    virtual void applyScaling() {
      double x0 = length_scaling;

      for (unsigned int j = 0; j < functions.size(); j++) {
        for (unsigned int i = 0; i < functions[j].phi.size(); i++) {
          functions[j].phi[i] *= x0;
          functions[j].dphidx[i] *= x0*x0;
          functions[j].dphidy[i] *= x0*x0;
          functions[j].dphidz[i] *= x0*x0;
        }
      }
    }
};

#endif
