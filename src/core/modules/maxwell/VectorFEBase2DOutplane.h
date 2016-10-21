#ifndef __VECTOR_FEBASE_OUTPLANE_2D_H__
#define __VECTOR_FEBASE_OUTPLANE_2D_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "VectorFunction.h"
#include "FunctionInfo.h"
#include "IFEBase.h"
#include "IVectorFEBase.h"

//TODO rewrite me mb

class VectorFEBase2DOutplane : public IVectorFEBase {
  private:
    ScalarFEBase2D fe;

  public:
    VectorFEBase2DOutplane(double scaling) : IVectorFEBase(2, scaling), fe(scaling) {
    }

    virtual void attach_quadrature_rule(libMesh::QBase* q) {
      fe.attach_quadrature_rule(q);
    }

    virtual void reinit(const Elem *elem, unsigned int order, const std::vector<Point> *const pts=NULL) {
      fe.length_scaling = length_scaling;
      fe.reinit(elem, order, pts);
      copyFunctions();
    }

    virtual void reinit(const Elem *elem, unsigned int order, const unsigned int side, const Real tolerance = libMesh::TOLERANCE) {
      fe.length_scaling = length_scaling;
      fe.reinit(elem, order, side, tolerance);
      copyFunctions();
    }

    virtual const std::vector<Real>& get_JxW() const {
      return fe.get_JxW();
    }

    virtual const std::vector<Point>& get_xyz() const {
      return fe.get_xyz();
    }

    virtual void copyFunctions() {
      functions_info.resize(0);
      functions.resize(0);

      functions_info = fe.functions_info;

      for (int i = 0; i < fe.functions.size(); i++) {
        ScalarFunction sf = fe.functions[i];
        VectorFunction function(Point(0), sf.phi.size());
        for (unsigned int qp = 0; qp < sf.phi.size(); qp++) {
          function.phi[qp] = Point(0, 0, sf.phi[qp]);
          function.dphidx[qp] = Point(0, 0, sf.grad[qp](0));
          function.dphidy[qp] = Point(0, 0, sf.grad[qp](1));;
          function.dphidz[qp] = 0;
        }
        functions.push_back(function);
      }
    }
};

#endif
