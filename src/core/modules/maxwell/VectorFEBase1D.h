#ifndef __VECTOR_FEBASE_1D_H__
#define __VECTOR_FEBASE_1D_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "VectorFunction.h"
#include "FunctionInfo.h"
#include "IFEBase.h"
#include "IVectorFEBase.h"

class VectorFEBase1D : public IVectorFEBase {
  public:
    static const unsigned int DIM = 1;

    VectorFEBase1D(double scaling) : IVectorFEBase(DIM, scaling) {
    }

  protected:
    virtual void addFunctions(const Elem *elem, const std::vector<Point>& pts, unsigned int order) {
      addSFunction(phi_functions[0], FunctionInfo::VERTEX, 0, 0);
      addSFunction(phi_functions[1], FunctionInfo::VERTEX, 1, 0);

      std::vector<ScalarFunction> legendrePolynoms;
      ScalarFunction::calculateScaledIntegrateLegendrePlynomials(phi_functions[0] - phi_functions[1], phi_functions[0] + phi_functions[1], legendrePolynoms, order + 1);
      for (unsigned int t_order = 0; t_order < order; t_order++) {
        addSFunction(legendrePolynoms[t_order + 2], FunctionInfo::EDGE, 0, 0);// Interior functions. We should not care about direction etc.
      }
    }

    virtual void addSFunction(const ScalarFunction& function, FunctionInfo::ItemType itemType, int itemId, int order) {

      VectorFunction fv(Point(0, 0, 0), function.phi.size());
      for (unsigned int qp = 0; qp < function.phi.size(); qp++) {
        fv.phi[qp] = Point(0, function.phi[qp], 0);
        fv.dphidx[qp] = Point(0, function.grad[qp](0), 0);
        fv.dphidy[qp] = 0;
        fv.dphidz[qp] = 0;
      }
      addFunction(fv, itemType, itemId, order);
    }
};

#endif
