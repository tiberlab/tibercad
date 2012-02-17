#include "VectorFEBase2D.h"
#include "fe_type.h"
#include <assert.h>
#include "VectorFunction.h"
#include "ElementUtils.h"

VectorFEBase2D::VectorFEBase2D(double scaling) : IVectorFEBase(DIM, scaling) {
}

void VectorFEBase2D::addTriFunctions(const Elem *elem, const std::vector<Point>& pts, int order) {
  /* Lower order (order == 0) */
  if (order >= 0) {
    VectorFunction nablal[THREE];
    if (pts.size() > 0) {
      // We use here the fact that nabla li is const for triangles.
      nablal[0] = VectorFunction(scalarFe->get_dphi()[0][0], pts.size());
      nablal[1] = VectorFunction(scalarFe->get_dphi()[1][0], pts.size());
      nablal[2] = VectorFunction(scalarFe->get_dphi()[2][0], pts.size());
    }

    for (unsigned int i = 0; i < THREE; i++) {
      int i1 = ElementUtils::getNext(i, THREE);
      int multiplier = ElementUtils::getDirection(elem, i);

      VectorFunction tmp = nablal[i] * phi_functions[i1];
      VectorFunction tmp2 = nablal[i1] * phi_functions[i] - tmp;

      addFunction(tmp2 * multiplier * 2 /* 2 is for compatible with quads*/, FunctionInfo::EDGE, i, 0);
    }
  }

  // From here we will need second order derivatives
  phi_functions[0].setSecondOrder(true);
  phi_functions[1].setSecondOrder(true);
  phi_functions[2].setSecondOrder(true);

  std::vector<std::vector<ScalarFunction>> scaledIntegrateLegendrePlynomials;
  scaledIntegrateLegendrePlynomials.resize(3);

  for (int i = 0; i < THREE; i++) {
    int i1 = ElementUtils::getNext(i, THREE);
    ScalarFunction::calculateScaledIntegrateLegendrePlynomials((phi_functions[i1] - phi_functions[i]) * ElementUtils::getDirection(elem, i),
        phi_functions[i] + phi_functions[i1],
        scaledIntegrateLegendrePlynomials[i], order + 1);
  }

  /* Edges (order == 1+) */
  for (int t_order = 2; t_order <= order + 1; t_order++) {
    for (int i = 0; i < THREE; i++) {
      addFunction(VectorFunction::gradient(scaledIntegrateLegendrePlynomials[i][t_order]), FunctionInfo::EDGE, i, t_order - 1);
    }
  }

  /*Bubbles */
  std::vector<ScalarFunction> legendrePlynomials;
  ScalarFunction tmp(-0.5, pts.size(), true);
  ScalarFunction::calculateLegendrePlynomials((tmp + phi_functions[2]) * 2, legendrePlynomials, order - 2);

  /* Bubbles type 3 */
  for (int t_order = 3; t_order <= order + 1; t_order++) {
    addFunction(functions[0] * phi_functions[2] * legendrePlynomials[t_order - 3], FunctionInfo::SURFACE, -1, t_order);
  }

  /* Bubbles type 1 and 2 */
  for (int t_i = 0; t_i <= order - 2; t_i++) {
    for (int t_j = 0; t_j <= (order - 2 - t_i); t_j++) {
      ScalarFunction vj = legendrePlynomials[t_j] * phi_functions[2];
      VectorFunction nablaUi = VectorFunction::gradient(scaledIntegrateLegendrePlynomials[0][t_i + 2]);
      VectorFunction uiNablaVj = VectorFunction::gradient(vj) * scaledIntegrateLegendrePlynomials[0][t_i + 2];
      addFunction(nablaUi * vj + uiNablaVj, FunctionInfo::SURFACE, -1, t_i + t_j + 3);
      addFunction(nablaUi * vj - uiNablaVj, FunctionInfo::SURFACE, -1, t_i + t_j + 3 );
    }
  }

  // Clean up)
  phi_functions[0].setSecondOrder(false);
  phi_functions[1].setSecondOrder(false);
  phi_functions[2].setSecondOrder(false);
}

void VectorFEBase2D::addQuadFunctions(const Elem *elem, const std::vector<Point>& pts, int order) {
  ScalarFunction function_x(pts, 0);
  ScalarFunction function_y(pts, 1);
  ScalarFunction function_one(1.0, pts.size());
  VectorFunction ex(Point(1, 0), pts.size());
  VectorFunction ey(Point(0, 1), pts.size());

  ////////////////////////////////////////////////////

  std::vector<ScalarFunction> legendrePolynoms_x, legendrePolynoms_y, legendrePolynoms_mx, legendrePolynoms_my;

  //TODO mb use 2 polynoms to save time...same for scalar
  ScalarFunction::calculateLegendrePlynomials(function_x * ElementUtils::getDirection(elem, 0), legendrePolynoms_x, order);
  ScalarFunction::calculateLegendrePlynomials(function_y * ElementUtils::getDirection(elem, 1), legendrePolynoms_y, order);
  ScalarFunction::calculateLegendrePlynomials(function_x * (-ElementUtils::getDirection(elem, 2)), legendrePolynoms_mx, order);
  ScalarFunction::calculateLegendrePlynomials(function_y * (-ElementUtils::getDirection(elem, 3)), legendrePolynoms_my, order);

  for (int t_order = 0; t_order <= order; t_order++) {
    double half = 0.5;

    addFunction(ex * legendrePolynoms_x[t_order] * (function_one - function_y) * (half * ElementUtils::getDirection(elem, 0)), FunctionInfo::EDGE, 0, t_order);
    addFunction(ey * legendrePolynoms_y[t_order] * (function_one + function_x) * (half * ElementUtils::getDirection(elem, 1)), FunctionInfo::EDGE, 1, t_order);
    addFunction(ex * legendrePolynoms_mx[t_order] * (function_one + function_y) * (- half * ElementUtils::getDirection(elem, 2)), FunctionInfo::EDGE, 2, t_order);
    addFunction(ey * legendrePolynoms_my[t_order] * (function_one - function_x) * (- half * ElementUtils::getDirection(elem, 3)), FunctionInfo::EDGE, 3, t_order);
  }

  std::vector<ScalarFunction> intergratedLegendrePolynoms_x, intergratedLegendrePolynoms_y;

  ScalarFunction::calculateIntegrateLegendrePlynomials(function_x, intergratedLegendrePolynoms_x, order + 1);
  ScalarFunction::calculateIntegrateLegendrePlynomials(function_y, intergratedLegendrePolynoms_y, order + 1);

  for (int j = 0; j < order + 1; j++) {
    for (int k = 2; k < order + 2; k++) {
      addFunction(ex * legendrePolynoms_x[j] * intergratedLegendrePolynoms_y[k], FunctionInfo::SURFACE, -1, k);
      addFunction(ey * legendrePolynoms_y[j] * intergratedLegendrePolynoms_x[k], FunctionInfo::SURFACE, -1, k);
    }
  }

  change_coordinate_system();
}

void VectorFEBase2D::change_coordinate_system() {
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
        functions[j].phi[qp] = changePoint(DF_T, functions[j].phi[qp], 1);

        // TODO: curl transforms in this way, div may be in another

        functions[j].dphidx[qp] = functions[j].dphidx[qp] / detDF_T;
        functions[j].dphidy[qp] = functions[j].dphidy[qp] / detDF_T;
        functions[j].dphidz[qp] = functions[j].dphidz[qp] / detDF_T;
      }
    }
  }
}
