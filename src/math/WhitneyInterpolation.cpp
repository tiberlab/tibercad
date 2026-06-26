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


#include "libmesh/elem.h"
#include "libmesh/fe_map.h"
#include "tibercad/base/libMeshDefs.h"

using namespace std;

void
WhitneyInterpolation::reinit(const libMesh::Elem& elem,
    const std::vector<libMesh::Point>& points)
{

  unsigned int dim = elem.dim();

  unsigned int nn = elem.n_nodes();
  unsigned int ne = elem.n_edges();
  if (dim == 1)
    ne = 1;

  unsigned int np = points.size();

  _w0.resize(nn);
  _w1.resize(ne);

  unique_ptr<libMesh::FEBase> fe = libMesh::FEBase::build(dim, libMesh::FEType(1, libMesh::LAGRANGE));
  const vector<vector<libMesh::Real>> &phi = fe->get_phi();
  const vector<vector<libMesh::RealGradient>> &dphi = fe->get_dphi();

  fe->reinit(&elem, &points);

  _w0 = phi;

  if (dim == 1)
  {
    // for Lagrange elements of order 1, the 1-form is constant along the edge, and its value
    // is given by the derivative of the basis function associated with the second node of the
    // edge (the first node's basis function derivative is negative of that).
    for (unsigned int p = 0; p < np; ++p)
      _w1[0][p] = dphi[1][0];
  }
  else
  {
    for (unsigned int i = 0; i < ne; ++i)
    {
      unsigned int ni = elem.local_edge_node(i, 0);
      unsigned int nj = elem.local_edge_node(i, 1);
      for (unsigned int p = 0; p < np; ++p)
      {
        _w1[i][p] = _w0[ni][p] * dphi[nj][p] - _w0[nj][p] * dphi[ni][p];
      }
    }
  }
}
