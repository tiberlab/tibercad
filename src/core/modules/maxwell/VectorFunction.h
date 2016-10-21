#ifndef __vector_function_h__
#define __vector_function_h__

#include "ScalarFunction.h"
#include <vector>
#include "point.h"
#include "vector_value.h"
#include "Messages.h"

class VectorFunction {
  public:
    std::vector<libMesh::Point> phi;
    std::vector<libMesh::Point> dphidx;
    std::vector<libMesh::Point> dphidy;
    std::vector<libMesh::Point> dphidz;

    VectorFunction() {
    }

    VectorFunction(libMesh::Point value, unsigned int size) {
      for (unsigned int i = 0; i < size; i++) {
        phi.push_back(value);
        dphidx.push_back(0);
        dphidy.push_back(0);
        dphidz.push_back(0);
      }
    }

    VectorFunction(const std::vector<libMesh::Point>& values, const std::vector<libMesh::Point>& d_dx, const std::vector<libMesh::Point>& d_dy, const std::vector<libMesh::Point>& d_dz, double multiplier) {
      for (unsigned int i = 0; i < values.size(); i++) {
        phi.push_back(multiplier * values[i]);
        dphidx.push_back(multiplier * d_dx[i]);
        dphidy.push_back(multiplier * d_dy[i]);
        dphidz.push_back(multiplier * d_dz[i]);
      }
    }

    VectorFunction(const VectorFunction* f1p, const VectorFunction& f2, double a1, double a2) {
      const VectorFunction f1 = *f1p;
      for (unsigned int i = 0; i < f1.phi.size(); i++) {
        phi.push_back(a1 * f1.phi[i] + a2*f2.phi[i]);
        dphidx.push_back(a1 * f1.dphidx[i] + a2*f2.dphidx[i]);
        dphidy.push_back(a1 * f1.dphidy[i] + a2*f2.dphidy[i]);
        dphidz.push_back(a1 * f1.dphidz[i] + a2*f2.dphidz[i]);
      }
    }

    VectorFunction(const VectorFunction* f1p, double a) {
      const VectorFunction f1 = *f1p;
      for (unsigned int i = 0; i < f1.phi.size(); i++) {
        phi.push_back(a * f1.phi[i]);
        dphidx.push_back(a * f1.dphidx[i]);
        dphidy.push_back(a * f1.dphidy[i]);
        dphidz.push_back(a * f1.dphidz[i]);
      }
    }

    VectorFunction(const VectorFunction* f1p, const ScalarFunction& f2) {
      const VectorFunction f1 = *f1p;
      for (int i = 0; i < f1.phi.size(); i++) {
        phi.push_back(f1.phi[i] * f2.phi[i]);
        dphidx.push_back(f1.dphidx[i] * f2.phi[i] + f1.phi[i] * f2.grad[i](0));
        dphidy.push_back(f1.dphidy[i] * f2.phi[i] + f1.phi[i] * f2.grad[i](1));
        dphidz.push_back(f1.dphidz[i] * f2.phi[i] + f1.phi[i] * f2.grad[i](2));

        /*
        div.push_back(f1.div[i] * f2.phi[i] + f1.phi[i] * f2.grad[i]);
        Point point(f2.grad[i](1) * f1.phi[i](2) - f2.grad[i](2) * f1.phi[i](1),
                    f2.grad[i](2) * f1.phi[i](0) - f2.grad[i](0) * f1.phi[i](2),
                    f2.grad[i](0) * f1.phi[i](1) - f2.grad[i](1) * f1.phi[i](0));
        curl.push_back(f1.curl[i] * f2.phi[i] + point);
*/
      }
    }


    VectorFunction operator+(const VectorFunction& f) {
      return VectorFunction(this, f, 1.0, 1.0);
    }

    VectorFunction operator-(const VectorFunction& f) {
      return VectorFunction(this, f, 1.0, -1.0);
    }

    VectorFunction operator*(double a) {
      return VectorFunction(this, a);
    }

    VectorFunction operator*(const ScalarFunction& f) {
      return VectorFunction(this, f);
    }

    static VectorFunction gradient(const ScalarFunction& scalarFunction) {
      if (!scalarFunction.secondOrderDerivatives) {
        Messages::error("Can not get gradient because second order derivatives is undefined.");
        return VectorFunction();
      }
      VectorFunction result(libMesh::Point(0, 0), scalarFunction.phi.size());
      for (int i = 0; i < scalarFunction.phi.size(); i++) {
        result.phi[i] = scalarFunction.grad[i];
        result.dphidx[i] = scalarFunction.dgraddx[i];
        result.dphidy[i] = scalarFunction.dgraddy[i];
        result.dphidz[i] = scalarFunction.dgraddz[i];
      }

      return result;
    }


    VectorValue<Complex> curl(unsigned int qp) const {
      return  VectorValue<Complex>(
          dphidy[qp](2) - dphidz[qp](1),
          dphidz[qp](0) - dphidx[qp](2),
          dphidx[qp](1) - dphidy[qp](0)
      );
    }
};

#endif
