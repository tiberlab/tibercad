#include "ScalarFEBase3D.h"
#include "fe_type.h"
#include "VectorFunction.h"
#include "ElementUtils.h"
using namespace libMesh;

ScalarFEBase3D::ScalarFEBase3D(double scaling) : IScalarFEBase(DIM, scaling) {
}

void ScalarFEBase3D::applyScaling() {
  double x0 = length_scaling;

  for (unsigned int j = 0; j < functions.size(); j++) {
    for (unsigned int i = 0; i < functions[j].phi.size(); i++) {
      functions[j].grad[i] *= x0;
    }
  }

}

//TODO do not compute polynoms 2 times...
void ScalarFEBase3D::addTetraFunctions(const Elem *elem, const std::vector<Point>& pts, int order) {
  for (unsigned int i = 0; i < phi_functions.size(); i++) {
    addFunction(phi_functions[i], FunctionInfo::VERTEX, i, 0);
  }

  if (elem->n_nodes() != 4) {
    return;
  }

  for (int i = 0; i < elem->n_edges(); i++) {
    std::vector<std::pair<unsigned int, unsigned int>> vertices = ElementUtils::getSortedEdge(elem, i);

    int i0 = vertices[0].first;
    int i1 = vertices[1].first;

    std::vector<ScalarFunction> LP;
    ScalarFunction::calculateScaledIntegrateLegendrePlynomials(phi_functions[i0] - phi_functions[i1],
        phi_functions[i0] + phi_functions[i1], LP, order + 1);

    for (int t_order = 2; t_order <= order + 1; t_order++) {
      addFunction(LP[t_order], FunctionInfo::EDGE, i, t_order - 1);
    }
  }

  // Face functions
  for (int side = 0; side < elem->n_sides(); side++) {
    std::vector<std::pair<unsigned int, unsigned int>> vertices = ElementUtils::getSortedSide(elem, side);

    std::vector<ScalarFunction> LP;

    ScalarFunction::calculateScaledIntegrateLegendrePlynomials(phi_functions[vertices[0].first] - phi_functions[vertices[1].first],
        phi_functions[vertices[0].first] + phi_functions[vertices[1].first], LP, order);

    std::vector<ScalarFunction> lsj;

    ScalarFunction::calculateScaledLegendrePlynomials(phi_functions[vertices[2].first] - phi_functions[vertices[0].first] - phi_functions[vertices[1].first],
        phi_functions[vertices[2].first] + phi_functions[vertices[0].first] + phi_functions[vertices[1].first], lsj, order - 2);


    for (int t_j = 0; t_j <= order - 2; t_j++) {
      ScalarFunction vj = lsj[t_j] * phi_functions[vertices[2].first];
      for (int t_i = 0; t_i + t_j <= order - 2; t_i++) {
        //addFunction(LP[t_i + 2] * vj, FunctionInfo::SURFACE, side, t_i * FunctionInfo::MAX_ORDER + t_j);
        addFunction(phi_functions[vertices[0].first] * phi_functions[vertices[1].first] * phi_functions[vertices[2].first], FunctionInfo::SURFACE, side, t_i * FunctionInfo::MAX_ORDER + t_j);
      }
    }
  }

  // Cell functions
  std::vector<ScalarFunction> LP;

  ScalarFunction::calculateScaledIntegrateLegendrePlynomials(phi_functions[0] - phi_functions[1],
      phi_functions[0] + phi_functions[1], LP, order - 1);

  std::vector<ScalarFunction> lsj;

  ScalarFunction::calculateScaledLegendrePlynomials(phi_functions[2] - phi_functions[1] - phi_functions[0],
      phi_functions[2] + phi_functions[1] + phi_functions[0], lsj, order - 3);

  std::vector<ScalarFunction> lk;

  ScalarFunction::calculateLegendrePlynomials(phi_functions[3] * 2 - ScalarFunction(1.0, pts.size()), lk, order - 3);

  for (int t_i = 0; t_i <= order - 3; t_i++) {
    for (int t_j = 0; t_j <= order - 3; t_j++) {
      for (int t_k = 0; t_k + t_j + t_i <= order - 3; t_k++) {
        addFunction(phi_functions[2] * phi_functions[3] * LP[t_i + 2] * lsj[t_j] * lk[t_k], FunctionInfo::VOLUME, 0, 0);
      }
    }
  }
}
