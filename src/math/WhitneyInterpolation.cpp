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

#include "libmesh/fe_interface.h"

using namespace std;

void
WhitneyInterpolation::reinit(const libMesh::Elem& elem,
    const std::vector<libMesh::Point>& points, bool local_coordinates)
{

  // check if element has changed. If so,
  // we need to recalculate the subdivison.
  if (_elem != &elem)
  {
    setup_1cells(elem, _primal_1cells);
    _elem = &elem;
  }

  unsigned int dim = elem.dim();

  unsigned int nn = elem.n_nodes();
  unsigned int ne = _primal_1cells.size();

  unsigned int np = points.size();

  _w1.resize(ne);

  unique_ptr<libMesh::FEBase> fe = libMesh::FEBase::build(dim, libMesh::FEType(1, libMesh::LAGRANGE));
  const vector<vector<libMesh::Real>> &phi = fe->get_phi();
  const vector<vector<libMesh::RealGradient>> &dphi = fe->get_dphi();

  vector<libMesh::Point> ref_points(np);
  if (local_coordinates)
  {
    ref_points = points;
    _xyz.resize(np);
    for (unsigned int p = 0; p < np; ++p)
      _xyz[p] = libMesh::FEMap::map(dim, &elem, ref_points[p]);
  }
  else
  {
    libMesh::FEMap::inverse_map(dim, &elem, points, ref_points);
    //libMesh::FEInterface::inverse_map(dim, libMesh::FEType(), &elem, points, ref_points);
    _xyz = points;
  }

  if (np > 0)
    fe->reinit(&elem, &ref_points);

  _w0 = phi;

  if (dim == 1)
  {
    // for Lagrange elements of order 1, the 1-form is constant along the edge, and its value
    // is given by the derivative of the basis function associated with the second node of the
    // edge (the first node's basis function derivative is negative of that).
    _w1[0].resize(np);
    for (unsigned int p = 0; p < np; ++p)
      _w1[0][p] = dphi[1][p];
  }
  else
  {
    // For simplices we can use the standard Whitney 1-forms.
    // For quands, we use a subdivison into two triangles.
    // For other non-simplices we use Whitney 1-forms for now, until we implement
    // a subdivision into simplices for these elements as well.
    if (elem.type() == libMesh::QUAD4)
    {
      // For quads, we use a subdivision into two triangles, and then use
      // the standard Whitney 1-forms on each triangle. The subdivison has
      // already been done in setup_1cells(), which added a virtual edge connecting
      // the two diagonally opposite nodes with larger sum of angles. The first
      // four primal 1-cells are the real edges of the quad, and the fifth
      // primal 1-cell is the virtual edge. We use the standard Whitney 1-forms
      // on each triangle.

      // Normal from the diagonals (unbiased w.r.t. any single vertex)
      const libMesh::Point N = (elem.point(2) - elem.point(0)).cross(elem.point(3) - elem.point(1));

      libMesh::Point p0 = elem.point(_primal_1cells[4].first);
      libMesh::Point p2 = elem.point(_primal_1cells[4].second);
      libMesh::Point v0 = p2 - p0;

      // First, group points according to which triangle they belong to. We use the diagonal as the dividing line.
      std::vector<unsigned int> left_points, right_points;
      for (unsigned int p = 0; p < np; ++p)
      {
        // Check if point is on the left or right of the diagonal. We use the cross product to determine this.
        libMesh::Point vp = _xyz[p] - p0;
        double cross = vp.cross(v0) * N;
        if (cross < 0)
          left_points.push_back(p);
        else
          right_points.push_back(p);
      }

      // Size the 1-forms for the 5 primal 1-cells (4 real edges + 1 virtual) to hold
      // all np points, since each triangle call below only fills in the subset of
      // points (and edges) belonging to that triangle.
      _w1.assign(5, std::vector<libMesh::RealGradient>(np));

      // For each triangle, we compute Whitney forms on the points that belong to that triangle.
      if (!left_points.empty())
      {
        unsigned int n1 = _primal_1cells[4].first;
        unsigned int n2 = _primal_1cells[4].second;
        unsigned int n3 = (n2 + 1) % 4; // the third node of the triangle is the next node in the quad

        calculate_subtriangle_whitney_forms(n1, n2, n3, left_points, ref_points);
      }

      if (!right_points.empty())
      {
        unsigned int n1 = _primal_1cells[4].second;
        unsigned int n2 = _primal_1cells[4].first;
        unsigned int n3 = (n2 + 1) % 4; // the third node of the triangle is the next node in the quad

        calculate_subtriangle_whitney_forms(n1, n2, n3, right_points, ref_points);
      }
    }
    else
    {
      // For simplices and other non-simplices, we use the standard Whitney 1-forms.
      // The Whitney 1-forms are defined as:
      // lambda_ij = N_i dN_j - N_j dN_i
      // where N_i is the basis function associated with node i, and dN_i is its gradient.
      for (unsigned int i = 0; i < ne; ++i)
      {
        _w1[i].resize(np);
        unsigned int ni = elem.local_edge_node(i, 0);
        unsigned int nj = elem.local_edge_node(i, 1);
        for (unsigned int p = 0; p < np; ++p)
        {
          _w1[i][p] = _w0[ni][p] * dphi[nj][p] - _w0[nj][p] * dphi[ni][p];
        }
      }
    }
  }
}


void WhitneyInterpolation::calculate_subtriangle_whitney_forms(unsigned int n1, unsigned int n2, unsigned int n3,
        const std::vector<unsigned int>& points, const std::vector<libMesh::Point>& ref_points)
{
  // parameters for the 0-forms, as a + bx + cy
  // area of the triangle in quad reference coordinates
  // we also use N3 = 1 - N1 - N2
  // NOTE: we use the hardcoded knowledge of the libMesh reference quad.
  const double area = 2.0;
  double a1 = 0.5 * ((_elem->master_point(n3)(1) - _elem->master_point(n2)(1)) * _elem->master_point(n3)(0) -
                     (_elem->master_point(n3)(0) - _elem->master_point(n2)(0)) * _elem->master_point(n3)(1)) / area;
  double b1 = 0.5 * (_elem->master_point(n2)(1) - _elem->master_point(n3)(1)) / area;
  double c1 = 0.5 * (_elem->master_point(n3)(0) - _elem->master_point(n2)(0)) / area;

  double a2 = 0.5 * ((_elem->master_point(n1)(1) - _elem->master_point(n3)(1)) * _elem->master_point(n3)(0) -
                     (_elem->master_point(n1)(0) - _elem->master_point(n3)(0)) * _elem->master_point(n3)(1)) / area;
  double b2 = 0.5 * (_elem->master_point(n3)(1) - _elem->master_point(n1)(1)) / area;
  double c2 = 0.5 * (_elem->master_point(n1)(0) - _elem->master_point(n3)(0)) / area;

  // the gradients are constant, we use the TRI reference coordinates to compute them
  // The reference element is (0,0), (1,0), (0,1) for the triangle. The gradients of the basis functions are:
  // grad(N1) = (-1, -1)
  // grad(N2) = (1, 0)
  // grad(N3) = (0, 1)

  // now we need the Jacobian of the mapping from the triangle reference coordinates to the real coordinates.
  libMesh::Point x1 = _elem->point(n1);
  libMesh::Point x2 = _elem->point(n2);
  libMesh::Point x3 = _elem->point(n3);

  const libMesh::Point e1 = x2 - x1;
  const libMesh::Point e2 = x3 - x1;

  const Real e1e1 = e1 * e1;
  const Real e2e2 = e2 * e2;
  const Real e1e2 = e1 * e2;

  const Real detG = e1e1 * e2e2 - e1e2 * e1e2; // = |e1 x e2|^2, i.e. (2*Area)^2

  const libMesh::Point col_xi = (e2e2 * e1 - e1e2 * e2) / detG;
  const libMesh::Point col_eta = (e1e1 * e2 - e1e2 * e1) / detG;

  libMesh::RealGradient grad_N1 = -(col_xi + col_eta);
  libMesh::RealGradient grad_N2 = col_xi;
  libMesh::RealGradient grad_N3 = col_eta;

  double sign = 1.0;
  if (n1 > n2)
    sign = -1.0;

  // NOTE: _w1 is sized to hold all points by the caller (reinit()); this function
  // only fills in the entries for the points (and edges) belonging to this triangle.
  for (unsigned int p : points)
  {
    double N1 = a1 + b1 * ref_points[p](0) + c1 * ref_points[p](1);
    double N2 = a2 + b2 * ref_points[p](0) + c2 * ref_points[p](1);
    double N3 = 1.0 - N1 - N2;
    _w1[4][p] = sign * (N1 * grad_N2 - N2 * grad_N1);
    _w1[n2][p] = N2 * grad_N3 - N3 * grad_N2;
    _w1[n3][p] = N3 * grad_N1 - N1 * grad_N3;
  }
}



void
WhitneyInterpolation::setup_1cells(const libMesh::Elem& elem,
    std::vector<std::pair<unsigned int, unsigned int>>& primal_1cells)
{
  primal_1cells.clear();
  unsigned int dim = elem.dim();
  unsigned int ne = (dim == 1) ? 1 : elem.n_edges();

  primal_1cells.reserve(ne);

  if (dim == 1)
  {
    primal_1cells.push_back(make_pair(0, 1));
  }
  else
  {
    for (unsigned int e = 0; e < ne; ++e)
    {
      unsigned int ni = elem.local_edge_node(e, 0);
      unsigned int nj = elem.local_edge_node(e, 1);
      primal_1cells.push_back(make_pair(ni, nj));
    }

    if (elem.type() == libMesh::QUAD4)
    {
      // for quadrilaterals, we add one virtual edge connecting
      // the two diagonally opposite nodes with larger sum of angles.
      
      primal_1cells.push_back(larger_angle_pair(elem.point(0), elem.point(1),
                                                elem.point(2), elem.point(3)));
    }
    else if (elem.type() == libMesh::HEX8)
    {
      // for hexahedra, we add the virtual edges of the six tetrahedra
      primal_1cells.push_back(make_pair(0, 6));
      primal_1cells.push_back(make_pair(1, 7));
      primal_1cells.push_back(make_pair(2, 4));
      primal_1cells.push_back(make_pair(3, 5));
      primal_1cells.push_back(make_pair(0, 5));
      primal_1cells.push_back(make_pair(1, 4));
    }
  }
}


std::pair<unsigned int, unsigned int>
WhitneyInterpolation::larger_angle_pair(const libMesh::Point & x1,
                                        const libMesh::Point & x2,
                                        const libMesh::Point & x3,
                                        const libMesh::Point & x4) const
{
  // Normal from the diagonals (unbiased w.r.t. any single vertex)
  const libMesh::Point N = (x3 - x1).cross(x4 - x2);

  // Vectors at x2 and x4
  const libMesh::Point a = x1 - x2, b = x3 - x2;
  const libMesh::Point c = x1 - x4, d = x3 - x4;

  const double P = a.cross(b) * N;   // ∝ |a||b| sin(angle2)
  const double Q = a * b;            // ∝ |a||b| cos(angle2)
  const double R = c.cross(d) * N;   // ∝ |c||d| sin(angle4)
  const double S = c * d;            // ∝ |c||d| cos(angle4)

  const double test = P * S + Q * R; // ∝ sin(angle2 + angle4)

  return (test > 0 ? make_pair(0, 2) : make_pair(1, 3)); 
}