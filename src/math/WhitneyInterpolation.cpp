/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file WhitneyInterpolation.cpp
 * \brief tiberCAD API implementation.
 */


#include "tibercad/math/WhitneyInterpolation.h"

#include <numeric>

using namespace std;

void
WhitneyInterpolation::reinit(const std::vector<libMesh::Point>& nodes,
    const std::vector<libMesh::Point>& points)
{
  _nodes = nodes;

  setup_1cells();
  compute_gradients();

  unsigned int nn = _nodes.size();
  unsigned int ne = _cells.size();
  unsigned int np = points.size();

  _w0.assign(nn, vector<double>(np));
  _w1.assign(ne, vector<libMesh::RealGradient>(np));

  // The barycentric coordinate functions (0-forms) are affine, and vanish
  // at every node but their own (where they take the value 1), so they
  // can be evaluated from their (constant) gradient alone.
  for (unsigned int p = 0; p < np; ++p)
  {
    libMesh::Point dx = points[p] - _nodes[0];

    double n0 = 1.0;
    for (unsigned int i = 1; i < nn; ++i)
    {
      double ni = _grad[i] * dx;
      _w0[i][p] = ni;
      n0 -= ni;
    }
    _w0[0][p] = n0;
  }

  // Standard Whitney 1-forms: lambda_ij = N_i dN_j - N_j dN_i
  for (unsigned int e = 0; e < ne; ++e)
  {
    unsigned int ni = _cells[e].first;
    unsigned int nj = _cells[e].second;
    for (unsigned int p = 0; p < np; ++p)
      _w1[e][p] = _w0[ni][p] * _grad[nj] - _w0[nj][p] * _grad[ni];
  }
}


void
WhitneyInterpolation::setup_1cells(void)
{
  _cells.clear();

  unsigned int nn = _nodes.size();
  _cells.reserve(nn * (nn - 1) / 2);

  for (unsigned int i = 0; i < nn; ++i)
    for (unsigned int j = i + 1; j < nn; ++j)
      _cells.push_back(make_pair(i, j));
}


void
WhitneyInterpolation::compute_gradients(void)
{
  unsigned int nn = _nodes.size();
  unsigned int dim = nn - 1;

  _grad.assign(nn, libMesh::RealGradient(0.0));

  // Edge vectors from node 0 to the other nodes.
  vector<libMesh::Point> e(dim);
  for (unsigned int i = 0; i < dim; ++i)
    e[i] = _nodes[i + 1] - _nodes[0];

  // The gradient of the barycentric coordinate N_{i+1} is obtained from
  // the (generalized, for a simplex not spanning the full ambient space)
  // inverse of the Gram matrix G_ab = e_a . e_b of the edge vectors, as
  // grad(N_{i+1}) = sum_j Ginv(i,j) * e_j. This reduces to the exact
  // inverse for a tetrahedron in 3D space, and to a least-squares
  // (Moore-Penrose) solution for an edge or a triangle embedded in 3D.
  if (dim == 1)
  {
    double g00 = e[0] * e[0];
    _grad[1] = e[0] / g00;
  }
  else if (dim == 2)
  {
    double g00 = e[0] * e[0];
    double g01 = e[0] * e[1];
    double g11 = e[1] * e[1];
    double det = g00 * g11 - g01 * g01;

    _grad[1] = (g11 * e[0] - g01 * e[1]) / det;
    _grad[2] = (g00 * e[1] - g01 * e[0]) / det;
  }
  else // dim == 3
  {
    double g00 = e[0] * e[0], g01 = e[0] * e[1], g02 = e[0] * e[2];
    double g11 = e[1] * e[1], g12 = e[1] * e[2], g22 = e[2] * e[2];

    // Cofactors of the (symmetric) Gram matrix.
    double c00 = g11 * g22 - g12 * g12;
    double c01 = g02 * g12 - g01 * g22;
    double c02 = g01 * g12 - g02 * g11;
    double c11 = g00 * g22 - g02 * g02;
    double c12 = g01 * g02 - g00 * g12;
    double c22 = g00 * g11 - g01 * g01;

    double det = g00 * c00 + g01 * c01 + g02 * c02;

    _grad[1] = (c00 * e[0] + c01 * e[1] + c02 * e[2]) / det;
    _grad[2] = (c01 * e[0] + c11 * e[1] + c12 * e[2]) / det;
    _grad[3] = (c02 * e[0] + c12 * e[1] + c22 * e[2]) / det;
  }

  // N_0 = 1 - sum_{i>0} N_i, so grad(N_0) = -sum_{i>0} grad(N_i)
  libMesh::RealGradient sum(0.0);
  for (unsigned int i = 1; i < nn; ++i)
    sum += _grad[i];
  _grad[0] = -sum;
}
