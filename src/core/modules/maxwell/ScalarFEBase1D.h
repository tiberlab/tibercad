#ifndef __SCALAR_FEBASE_1D_H__
#define __SCALAR_FEBASE_1D_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "VectorFunction.h"
#include "FunctionInfo.h"
#include "IFEBase.h"
#include "IScalarFEBase.h"

class ScalarFEBase1D : public IScalarFEBase {
  public:
    static const unsigned int DIM = 1;

    ScalarFEBase1D(double scaling) : IScalarFEBase(DIM, scaling) {
    }

  protected:

    virtual void applyScaling() {
      double x0 = length_scaling;

      for (unsigned int j = 0; j < functions.size(); j++) {
        for (unsigned int i = 0; i < functions[j].phi.size(); i++) {
          functions[j].grad[i] *= x0;
        }
      }
    }

    virtual void addFunctions(const Elem *elem, const std::vector<Point>& pts, unsigned int order) {
      addFunction(phi_functions[0], FunctionInfo::VERTEX, 0, 0);
      addFunction(phi_functions[1], FunctionInfo::VERTEX, 1, 0);

      std::vector<ScalarFunction> legendrePolynoms;
      ScalarFunction::calculateScaledIntegrateLegendrePlynomials(phi_functions[0] - phi_functions[1], phi_functions[0] + phi_functions[1], legendrePolynoms, order + 1);
      for (unsigned int t_order = 0; t_order < order; t_order++) {
        addFunction(legendrePolynoms[t_order + 2], FunctionInfo::EDGE, 0, 0);// Interior functions. We should not care about direction etc.
      }
    }
};

#endif
