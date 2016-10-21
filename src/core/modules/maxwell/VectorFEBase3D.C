#include "VectorFEBase3D.h"
#include "fe_type.h"
#include "VectorFunction.h"
#include "ElementUtils.h"


VectorFEBase3D::VectorFEBase3D(double scaling) : IVectorFEBase(DIM, scaling) {
}

void VectorFEBase3D::addTetraFunctions(const libMesh::Elem *elem, const std::vector<libMesh::Point>& pts, int order) {

  //0,1 1,2 2,0 0,3 1,3 2,3
  VectorFunction nablal[elem->n_nodes()];

  for (int ii = 0; ii < elem->n_nodes(); ii++) {
    phi_functions[ii].setSecondOrder(true);
    // We use here the fact that nabla li is const fro tetra.
    if (pts.size() > 0) {
      nablal[ii] = VectorFunction(scalarFe->get_dphi()[ii][0], pts.size());
    }
  }

  // Edge based.

  for (int i = 0; i < elem->n_edges(); i++) {
    std::vector<std::pair<unsigned int, unsigned int>> vertices = ElementUtils::getSortedEdge(elem, i);

    int i0 = vertices[0].first;
    int i1 = vertices[1].first;

    addFunction(nablal[i0] * phi_functions[i1] - nablal[i1] * phi_functions[i0], FunctionInfo::EDGE, i, 0);

    std::vector<ScalarFunction> LP;
    ScalarFunction::calculateScaledIntegrateLegendrePlynomials(phi_functions[i0] - phi_functions[i1],
        phi_functions[i0] + phi_functions[i1], LP, order + 1);

    for (int t_order = 2; t_order <= order + 1; t_order++) {
      addFunction(VectorFunction::gradient(LP[t_order]), FunctionInfo::EDGE, i, t_order - 1);
    }
  }

  if (elem->n_nodes() != 4) {
    return; //TODO
  }

  //Face based.
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
        // order must be unique!

        addFunction(VectorFunction::gradient(LP[t_i + 2] * vj), FunctionInfo::SURFACE, side, t_i * FunctionInfo::MAX_ORDER + t_j);
        //addFunction(nablal[vertices[2].first] * phi_functions[vertices[0].first] * phi_functions[vertices[1].first], FunctionInfo::SURFACE, side, t_i * FunctionInfo::MAX_ORDER + t_j);

        addFunction(VectorFunction::gradient(LP[t_i + 2]) * vj - VectorFunction::gradient(vj) * LP[t_i + 2], FunctionInfo::SURFACE, side, FunctionInfo::MAX_ORDER * FunctionInfo::MAX_ORDER + t_i * FunctionInfo::MAX_ORDER + t_j);
        //addFunction(nablal[vertices[0].first] * phi_functions[vertices[1].first] * phi_functions[vertices[2].first], FunctionInfo::SURFACE, side,  FunctionInfo::MAX_ORDER * FunctionInfo::MAX_ORDER + t_i * FunctionInfo::MAX_ORDER + t_j);
      }
      //addFunction(nablal[vertices[1].first] * phi_functions[vertices[0].first] * phi_functions[vertices[2].first], FunctionInfo::SURFACE, side, 2 *  FunctionInfo::MAX_ORDER * FunctionInfo::MAX_ORDER + t_j);
      //addFunction((nablal[vertices[0].first] * phi_functions[vertices[1].first] - nablal[vertices[1].first] * phi_functions[vertices[0].first]) * vj, FunctionInfo::SURFACE, side, 2 * FunctionInfo::MAX_ORDER * FunctionInfo::MAX_ORDER + t_j);
      addFunction((VectorFunction::gradient(phi_functions[vertices[0].first]) * phi_functions[vertices[1].first] - VectorFunction::gradient(phi_functions[vertices[1].first]) * phi_functions[vertices[0].first]) * vj, FunctionInfo::SURFACE, side, 2 * FunctionInfo::MAX_ORDER * FunctionInfo::MAX_ORDER + t_j);
    }
  }


  //Cell based
  std::vector<ScalarFunction> LP;

  ScalarFunction::calculateScaledIntegrateLegendrePlynomials(phi_functions[0] - phi_functions[1],
      phi_functions[0] + phi_functions[1], LP, order - 1);

  std::vector<ScalarFunction> lsj;

  ScalarFunction::calculateScaledLegendrePlynomials(phi_functions[2] - phi_functions[1] - phi_functions[0],
      phi_functions[2] + phi_functions[1] + phi_functions[0], lsj, order - 3);

  std::vector<ScalarFunction> lk;

  ScalarFunction::calculateLegendrePlynomials(phi_functions[3] * 2 - ScalarFunction(1.0, pts.size(), true), lk, order - 3);

  for (int i = 0; i <= order - 3; i++) {
    for (int j = 0; j <= order - 3; j++) {
      ScalarFunction vj = lsj[j] * phi_functions[2];

      for (int k = 0; k + j + i <= order - 3; k++) {
        ScalarFunction wk = lk[k] * phi_functions[3];

        addFunction(VectorFunction::gradient(LP[i + 2] * vj * wk), FunctionInfo::VOLUME, 0, 0);
        addFunction(VectorFunction::gradient(LP[i + 2]) * vj * wk - VectorFunction::gradient(vj) * LP[i + 2] * wk + VectorFunction::gradient(wk) * LP[i + 2] * vj, FunctionInfo::VOLUME, 0, 0);
        addFunction(VectorFunction::gradient(LP[i + 2]) * vj * wk + VectorFunction::gradient(vj) * LP[i + 2] * wk - VectorFunction::gradient(wk) * LP[i + 2] * vj, FunctionInfo::VOLUME, 0, 0);
        if (i == 0) {
          addFunction((nablal[0] * phi_functions[1] - nablal[1] * phi_functions[0]) * vj * wk, FunctionInfo::VOLUME, 0, 0);
        }
      }
    }
  }

  if (pts.size() > 0) {
    for (int i = 0; i < elem->n_nodes(); i++) {
      phi_functions[i].setSecondOrder(false);
    }
  }

}
