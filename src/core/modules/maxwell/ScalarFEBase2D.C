#include "ScalarFEBase2D.h"
#include "fe_type.h"
#include "VectorFunction.h"
#include "ElementUtils.h"
using namespace libMesh;

ScalarFEBase2D::ScalarFEBase2D(double scaling) : IScalarFEBase(DIM, scaling) {
}

void ScalarFEBase2D::addTriFunctions(const Elem *elem, const std::vector<Point>& pts, int order) {
  std::vector<int> verticesIds;
  ElementUtils::getVertecesIds(elem, verticesIds);

  if (order >= 0) {
    for (int i = 0; i < 3; i++) {
      addFunction(phi_functions[i], FunctionInfo::VERTEX, i, 0);
    }
  }

  for (int i = 0; i < 3; i++) {
    int i1 = ElementUtils::getNext(i, 3);
    int i2 = ElementUtils::getNext(i1, 3);

    std::vector<ScalarFunction> LP;
    ScalarFunction test1 = phi_functions[i] + phi_functions[i1];
    ScalarFunction test2 = (phi_functions[i1] - phi_functions[i]) * ElementUtils::getDirection(elem, i);

    ScalarFunction::calculateScaledIntegrateLegendrePlynomials(test2, test1, LP, order + 1);

    for (int t_order = 2; t_order <= order + 1; t_order++) {
      ScalarFunction result = LP[t_order];
      addFunction(result, FunctionInfo::EDGE, i, t_order - 1);
    }

    if (i == 0) {
      std::vector<ScalarFunction> integrated3Polynoms;
      ScalarFunction test3(-0.5, pts.size());
      test3 = (test3 + phi_functions[i2]) * 2;
      ScalarFunction::calculateLegendrePlynomials(test3, integrated3Polynoms, order - 2);

      for (int t_i = 0; t_i <= order - 2; t_i++) {
        for (int t_j = 0; t_j + t_i <= order - 2; t_j++) {
          ScalarFunction result = LP[t_i + 2] * phi_functions[i2] * integrated3Polynoms[t_j];
          addFunction(result, FunctionInfo::SURFACE, -1, 0);
        }
      }
    }
  }
}

void ScalarFEBase2D::addQuadFunctions(const Elem *elem, const std::vector<Point>& pts, int order) {
  ScalarFunction function_x(pts, 0);
  ScalarFunction function_y(pts, 1);
  ScalarFunction function_one(1.0, pts.size());
  ScalarFunction function_half(0.5, pts.size());

  if (order >= 0) {
    addFunction((function_one - function_x)*(function_one - function_y) * 0.25, FunctionInfo::VERTEX, 0, 0);
    addFunction((function_one + function_x)*(function_one - function_y) * 0.25, FunctionInfo::VERTEX, 1, 0);
    addFunction((function_one + function_x)*(function_one + function_y) * 0.25, FunctionInfo::VERTEX, 2, 0);
    addFunction((function_one - function_x)*(function_one + function_y) * 0.25, FunctionInfo::VERTEX, 3, 0);
  }

  std::vector<ScalarFunction> legendrePolynoms_1, legendrePolynoms_2, legendrePolynoms_3, legendrePolynoms_4;

  ScalarFunction::calculateIntegrateLegendrePlynomials(function_x * ElementUtils::getDirection(elem, 0), legendrePolynoms_1, order + 1);
  ScalarFunction::calculateIntegrateLegendrePlynomials(function_y * ElementUtils::getDirection(elem, 1), legendrePolynoms_2, order + 1);
  ScalarFunction::calculateIntegrateLegendrePlynomials(function_x * ElementUtils::getDirection(elem, 2) * (-1.0), legendrePolynoms_3, order + 1);
  ScalarFunction::calculateIntegrateLegendrePlynomials(function_y * ElementUtils::getDirection(elem, 3) * (-1.0), legendrePolynoms_4, order + 1);

  for (int i = 0; i < order; i++) {
    addFunction((function_one - function_y) * legendrePolynoms_1[i + 2] * 0.5, FunctionInfo::EDGE, 0, i + 1);
    addFunction((function_one + function_x) * legendrePolynoms_2[i + 2] * 0.5, FunctionInfo::EDGE, 1, i + 1);
    addFunction((function_one + function_y) * legendrePolynoms_3[i + 2] * 0.5, FunctionInfo::EDGE, 2, i + 1);
    addFunction((function_one - function_x) * legendrePolynoms_4[i + 2] * 0.5, FunctionInfo::EDGE, 3, i + 1);
  }

  for (int i = 0; i < order; i++) {
    for (int j = 0; j < order; j++) {
      addFunction(legendrePolynoms_1[i + 2] * legendrePolynoms_2[j + 2], FunctionInfo::SURFACE, -1, 0);//Order for interior functions may be any
    }
  }

  change_coordinate_system();
}

void ScalarFEBase2D::change_coordinate_system() {
  unsigned int ptsCount = functions.size() > 0 ? functions[0].phi.size() : 0;
  if (ptsCount > 0) {
    const std::vector<RealGradient>& dxyzdxi = scalarFe->get_dxyzdxi();
    const std::vector<RealGradient>& dxyzdeta = scalarFe->get_dxyzdeta();

    //For each point
    for (unsigned int qp = 0; qp < ptsCount; qp++) {
      double DF_T[2][2];

      double detDF_T = dxyzdeta[qp](1) * dxyzdxi[qp](0) - dxyzdxi[qp](1) * dxyzdeta[qp](0);

      DF_T[0][0] = dxyzdeta[qp](1) / detDF_T;
      DF_T[0][1] = - dxyzdxi[qp](1) / detDF_T;
      DF_T[1][0] = - dxyzdeta[qp](0) / detDF_T;
      DF_T[1][1] = dxyzdxi[qp](0) / detDF_T;

      //For each function
      for (unsigned int j = 0; j < functions.size(); j++) {
        Point newPoint = changePoint(DF_T, Point(functions[j].grad[qp](0), functions[j].grad[qp](1)), 1);

        functions[j].grad[qp] = newPoint;
      }
    }
  }
}
