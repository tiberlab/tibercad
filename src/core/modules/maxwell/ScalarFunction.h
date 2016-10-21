#ifndef __scalar_function_h__
#define __scalar_function_h__

#include <vector>
#include "vector_value.h"

class ScalarFunction {
  public:
    bool secondOrderDerivatives;

    std::vector<double> phi;
    std::vector<libMesh::Point> grad;

    std::vector<libMesh::Point> dgraddx;
    std::vector<libMesh::Point> dgraddy;
    std::vector<libMesh::Point> dgraddz;

    bool getSecondOrderDerivatives() const {
        return secondOrderDerivatives;
    }

    ScalarFunction(bool sOrderDerivatives = false) : secondOrderDerivatives(sOrderDerivatives)  {
    }

    ScalarFunction(double value, unsigned int size, bool sOrderDerivatives = false) : secondOrderDerivatives(sOrderDerivatives) {
      for (unsigned int i = 0; i < size; i++) {
        phi.push_back(value);
        grad.push_back(0);
        if (secondOrderDerivatives) {
          dgraddx.push_back(0);
          dgraddy.push_back(0);
          dgraddz.push_back(0);
        }
      }
    }

    ScalarFunction(const std::vector<double>& values, const std::vector<double>& d_dx, const std::vector<double>& d_dy, const std::vector<double>& d_dz, double multiplier, bool sOrderDerivatives = false) : secondOrderDerivatives(sOrderDerivatives) {
      for (unsigned int i = 0; i < values.size(); i++) {
        phi.push_back(multiplier * values[i]);
        grad.push_back(multiplier * libMesh::Point(d_dx[i], d_dy[i], d_dz[i]));
        if (secondOrderDerivatives) {
          dgraddx.push_back(0);
          dgraddy.push_back(0);
          dgraddz.push_back(0);
        }
      }
    }

    ScalarFunction(const ScalarFunction* f1p, const ScalarFunction& f2, double a1, double a2) : secondOrderDerivatives(f1p->secondOrderDerivatives) {
      const ScalarFunction f1 = *f1p;
      for (unsigned int i = 0; i < f1.phi.size(); i++) {
        phi.push_back(a1 * f1.phi[i] + a2*f2.phi[i]);
        grad.push_back(a1 * f1.grad[i] + a2*f2.grad[i]);
        if (secondOrderDerivatives) {
          dgraddx.push_back(a1 * f1.dgraddx[i] + a2*f2.dgraddx[i]);
          dgraddy.push_back(a1 * f1.dgraddy[i] + a2*f2.dgraddy[i]);
          dgraddz.push_back(a1 * f1.dgraddz[i] + a2*f2.dgraddz[i]);
        }
      }
    }

    ScalarFunction(const ScalarFunction* f1p, double a1) : secondOrderDerivatives(f1p->secondOrderDerivatives) {
      const ScalarFunction f1 = *f1p;
      for (unsigned int i = 0; i < f1.phi.size(); i++) {
        phi.push_back(a1 * f1.phi[i]);
        grad.push_back(a1 * f1.grad[i]);
        if (secondOrderDerivatives) {
          dgraddx.push_back(a1 * f1.dgraddx[i]);
          dgraddy.push_back(a1 * f1.dgraddy[i]);
          dgraddz.push_back(a1 * f1.dgraddz[i]);
        }
      }
    }

    void setSecondOrder(bool flag) {
      int newSize = flag ? phi.size() : 0;
      dgraddx.resize(newSize);
      dgraddy.resize(newSize);
      dgraddz.resize(newSize);
      secondOrderDerivatives = flag;
    }

    ScalarFunction(const ScalarFunction* f1p, const ScalarFunction& f2) : secondOrderDerivatives(f1p->secondOrderDerivatives) {
      const ScalarFunction f1 = *f1p;
      for (unsigned int i = 0; i < f1.phi.size(); i++) {
        phi.push_back(f1.phi[i] * f2.phi[i]);
        grad.push_back(f1.grad[i] * f2.phi[i] + f1.phi[i] * f2.grad[i]);
        if (secondOrderDerivatives) {
          //TT
          dgraddx.push_back(f1.dgraddx[i] * f2.phi[i] + f2.dgraddx[i] * f1.phi[i] + f1.grad[i] * f2.grad[i](0) + f1.grad[i](0) * f2.grad[i]);
          dgraddy.push_back(f1.dgraddy[i] * f2.phi[i] + f2.dgraddy[i] * f1.phi[i] + f1.grad[i] * f2.grad[i](1) + f1.grad[i](1) * f2.grad[i]);
          dgraddz.push_back(f1.dgraddz[i] * f2.phi[i] + f2.dgraddz[i] * f1.phi[i] + f1.grad[i] * f2.grad[i](2) + f1.grad[i](2) * f2.grad[i]);
        }
      }
    }

    ScalarFunction(const std::vector<libMesh::Point>& values, int coordNum, bool sOrderDerivatives = false) : secondOrderDerivatives(sOrderDerivatives) {
      for (unsigned int i = 0; i < values.size(); i++) {
        phi.push_back(values[i](coordNum));
        if (coordNum == 0) {
          grad.push_back(libMesh::Point(1, 0, 0));
        } else if (coordNum == 1) {
          grad.push_back(libMesh::Point(0, 1, 0));
        } else {
          grad.push_back(libMesh::Point(0, 0, 1));
        }

        if (secondOrderDerivatives) {
          dgraddx.push_back(0);
          dgraddy.push_back(0);
          dgraddz.push_back(0);
        }
      }
    }

    ScalarFunction operator+(const ScalarFunction& f) const {
      return ScalarFunction(this, f, 1.0, 1.0);
    }

    ScalarFunction operator-(const ScalarFunction& f) const {
      return ScalarFunction(this, f, 1.0, -1.0);
    }

    ScalarFunction operator*(const double a) const {
      return ScalarFunction(this, a);
    }

    ScalarFunction operator*(const ScalarFunction& f) const {
      return ScalarFunction(this, f);
    }

    // n is int, because sometimes it is suitable to pass negative values here
    // which has no meaning... but if is unsigned: unsigned = 0 - 1 >> 1 -)

    static void calculateLegendrePlynomials(const ScalarFunction& base, std::vector<ScalarFunction>& result, int n) {
       result.resize(0);
       result.push_back(ScalarFunction(1.0, base.phi.size(), base.secondOrderDerivatives)); //1
       result.push_back(base * 1.0); //x

       for (int p = 1; p < n; p++) {
         ScalarFunction tmp = (result[p - 1]* p);
         ScalarFunction res = (result[p] * base * (2*p+1)) - tmp;

         res = res * (1.0 / (p + 1.0));
         result.push_back(res);
       }
    }

    static void calculateIntegrateLegendrePlynomials(const ScalarFunction& base, std::vector<ScalarFunction>& result, int n) {
       result.resize(0);
       result.push_back(ScalarFunction(0.0, base.phi.size(), base.secondOrderDerivatives)); //Just smth :)
       result.push_back(base * 1.0); //x

       ScalarFunction tmp(0.5, base.phi.size(), base.secondOrderDerivatives);
       result.push_back(base * base * 0.5 - tmp); //1/2(x*x - 1)

       for (int p = 2; p < n; p++) {
         ScalarFunction tmp = result[p - 1]* (p - 2);
         ScalarFunction res = (result[p] * base * (2*p-1)) - tmp;

         res = res * (1.0 / (p + 1.0));
         result.push_back(res);
       }
    }

    static void calculateScaledIntegrateLegendrePlynomials(const ScalarFunction& base1, const ScalarFunction& base2, std::vector<ScalarFunction>& result, int n) {
      result.resize(0);
      result.push_back(ScalarFunction(0.0, base1.phi.size(), base1.secondOrderDerivatives)); //Just smth :)
      result.push_back(base1 * 1.0); //x

      ScalarFunction tmp = base2 * base2 * 0.5;
      result.push_back(base1 * base1 * 0.5 - tmp); //1/2(x*x - t*t)

      for (int p = 2; p < n; p++) {
        ScalarFunction tmp2 = base2 * base2 * result[p - 1] * (p - 2);
        ScalarFunction res = (result[p] * base1 * (2*p-1)) - tmp2;

        res = res * (1.0 / (p + 1.0));
        result.push_back(res);
      }
    }

    static void calculateScaledLegendrePlynomials(const ScalarFunction& base1, const ScalarFunction& base2, std::vector<ScalarFunction>& result, int n) {
      result.resize(0);
      result.push_back(ScalarFunction(1.0, base1.phi.size(), base1.secondOrderDerivatives));
      result.push_back(base1 * 1.0); //x

      for (int p = 1; p < n; p++) {
        ScalarFunction res = (result[p] * base1 * (2*p+1)) - (result[p - 1]* p) * base2 * base2;

        res = res * (1.0 / (p + 1.0));
        result.push_back(res);
      }
   }

    libMesh::VectorValue<libMesh::Complex> grads(unsigned int qp,  libMesh::VectorValue<libMesh::Complex> sVector) const {
     return  libMesh::VectorValue<libMesh::Complex>(grad[qp](0)*sVector(0), grad[qp](1)*sVector(1), grad[qp](2)*sVector(2));
   }
};

#endif
