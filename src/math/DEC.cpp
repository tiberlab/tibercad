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
 * \file DEC.cpp
 * \brief tiberCAD API implementation.
 */


#include "tibercad/math/DEC.h"
#include "tibercad/io/Messages.h"


#include "libmesh/elem.h"
#include "libmesh/dense_vector.h"
#include "libmesh/fe_map.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"

#include <fstream>
#include <cassert>
#include <set>
#include <numeric>
#include <limits>

using namespace std;
using namespace libMesh;


DEC::DEC(const libMesh::Elem& elem,
         DEC::DualConstruction dual_constr,
         DEC::HodgeConstruction hodge_constr)
: _elem(&elem),
  _dual_constr(dual_constr),
  _hodge_constr(hodge_constr) 
{
}

void
DEC::reinit(const libMesh::Elem& elem)
{
  _elem = &elem;

  init();
}

void
DEC::init(void)
{
  const libMesh::Elem& elem = *_elem;

  _center = get_center(elem);

  setup_subdivision();

  // Keyed by sorted node pair, since _primal_1cells[e] may store either
  // orientation of the edge (see setup_subdivision()).
  _cell_index.clear();
  for (unsigned int e = 0; e < _primal_1cells.size(); ++e)
  {
    std::pair<unsigned int, unsigned int> key = _primal_1cells[e];
    if (key.first > key.second)
      std::swap(key.first, key.second);
    _cell_index[key] = e;
  }

  unsigned int dim = elem.dim();
  unsigned int ne = (dim == 1) ? 1 : elem.n_edges();
  unsigned int nn = elem.n_nodes();

  // For non-simplex elements, Whitney interpolation for 1-forms
  // might require additional virtual 1-cells, which are not edges of the element.
  // If the Hodge is constructed using interpolation,
  // we need to augment the Hodge dimension accordingly. If we use MFD
  // construction, only real element edges are considered.
  if (_hodge_constr == INTERPOLATION)
  {
    ne = _primal_1cells.size();
  }

  _primal.resize(ne);
  _midpoints.resize(ne);
  _incidence.resize(ne, nn);
  _dual_volumes.resize(nn, 0.0);

  if (dim == 0)
  {
    _dual_volumes[0] = 1.0;
  }
  else if (dim == 1)
  {
    _primal[0] = elem.point(1) - elem.point(0);
    _midpoints[0] = 0.5 * (elem.point(0) + elem.point(1));
    _incidence(0, 0) = -1;
    _incidence(0, 1) =  1;

    double vol = 0.5 * elem.volume();
    _dual_volumes[0] = vol;
    _dual_volumes[1] = vol;
  }
  else
  {
    for (unsigned int e = 0; e < ne; ++e)
    {
      unsigned int ni = _primal_1cells[e].first;
      unsigned int nj = _primal_1cells[e].second;
      _primal[e] = (elem.point(nj) - elem.point(ni));
      _midpoints[e] = 0.5 * (elem.point(ni) + elem.point(nj));

      _incidence(e, ni) = -1;
      _incidence(e, nj) =  1;
    }

    // We need only the real edges of the element for dual volume
    // calculation, even if we use interpolation for Hodge construction.
    // This is due to the fact that we use a single center point, i.e.
    // subdivision is only for interpolation.
    for (unsigned int e = 0; e < elem.n_edges(); ++e)
    {
      unsigned int ni = elem.local_edge_node(e, 0);
      unsigned int nj = elem.local_edge_node(e, 1);
      Point a = _midpoints[e] - elem.point(ni);
      Point b = _center - elem.point(ni);
      Point axb = a.cross(b);

      double vol = 0.0;

      if (dim == 2)
      {
        // volume is the area of the quadilateral formed by two
        // half-edges and the center point, which can be divided
        // into two triangles. This is the contribution
        // of the triangle formed with edge e.
        vol = 0.5 * axb.norm();
      }
      else
      {
        // volume is the sum of cones, each cone has volume
        // 1/6 a * (b x c) where a is the vector from the
        // edge midpoint to the ceenter, b is the vector from
        // the edge midpoint to the center of the first side,
        // and c is the half-edge.

        // now we need the two sides containing the edge, and their centers
        unsigned int ns = elem.n_sides();
        for (unsigned int s = 0; s < ns; ++s)
        {
          if (elem.is_edge_on_side(e, s))
          {
            Point c = elem.side_ptr(s)->vertex_average() - elem.point(ni);
            double part_vol = 1.0 / 6.0 * a * (b.cross(c));
            vol += std::abs(part_vol);
          }
        }
      }

      _dual_volumes[ni] += vol;
      _dual_volumes[nj] += vol;
    }
  }
}


void
DEC::reinit_forms(const std::vector<libMesh::Point>& points, bool reference_coord)
{
  const libMesh::Elem& elem = *_elem;

  unsigned int dim = elem.dim();
  unsigned int nn = elem.n_nodes();
  unsigned int nc = _primal_1cells.size();
  unsigned int np = points.size();

  if (reference_coord)
  {
    _xyz.resize(np);
    for (unsigned int p = 0; p < np; ++p)
      _xyz[p] = libMesh::FEMap::map(dim, &elem, points[p]);
  }
  else
  {
    _xyz = points;
  }

  _w0.assign(nn, std::vector<double>(np, 0.0));
  _w1.assign(nc, std::vector<RealGradient>(np));

  bool is_simplex = (nn == dim + 1);

  if (!is_simplex && np > 0)
  {
    // For a logically subdivided (non-simplicial) element, the 0-forms
    // use the element's own native (e.g. bilinear/trilinear) Lagrange
    // basis functions, which are smooth over the whole element -- unlike
    // the 1-forms, which are necessarily piecewise over the logical
    // simplicial subdivision (see setup_subdivision()), since Whitney
    // forms are only defined on simplices.
    std::vector<libMesh::Point> ref_pts;
    libMesh::FEMap::inverse_map(dim, &elem, _xyz, ref_pts);

    libMesh::FEType fe_type(1, libMesh::LAGRANGE);
    for (unsigned int p = 0; p < np; ++p)
      for (unsigned int i = 0; i < nn; ++i)
        _w0[i][p] = libMesh::FEInterface::shape(fe_type, &elem, i, ref_pts[p]);
  }

  for (unsigned int p = 0; p < np; ++p)
  {
    // Find the sub-simplex containing the point: for each candidate, the
    // point belongs to it if all of its barycentric 0-forms are
    // non-negative. In case of round-off at a sub-simplex boundary, we
    // keep the candidate with the smallest (least negative) violation.
    // The Whitney interpolation for the best candidate found so far is
    // kept around, so it does not need to be recomputed once the search
    // is over (in particular, a plain simplex element never subdivides,
    // so the search loop below runs (and reinits _whip) only once).
    double best_violation = std::numeric_limits<double>::max();

    for (unsigned int s = 0; s < _sub_simplices.size(); ++s)
    {
      const std::vector<unsigned int>& sub = _sub_simplices[s];

      std::vector<Point> sub_nodes(sub.size());
      for (unsigned int a = 0; a < sub.size(); ++a)
        sub_nodes[a] = elem.point(sub[a]);

      _whip.reinit(sub_nodes, {_xyz[p]});

      const auto& w0 = _whip.get_0forms();

      double violation = 0.0;
      for (unsigned int a = 0; a < sub.size(); ++a)
        violation = std::max(violation, -w0[a][0]);

      if (violation < best_violation)
      {
        best_violation = violation;

        // Record the results for the current best candidate; overwritten
        // if a later candidate turns out to be a better fit.
        const auto& w1 = _whip.get_1forms();
        const auto& cells = _whip.get_1cells();

        // For a simplex, the sub-simplex *is* the whole element, so its
        // barycentric 0-forms already are the element's native basis
        // functions. For a non-simplex, _w0 was already filled above
        // using the native basis, which is what should be returned.
        if (is_simplex)
          for (unsigned int a = 0; a < sub.size(); ++a)
            _w0[sub[a]][p] = w0[a][0];

        for (unsigned int c = 0; c < cells.size(); ++c)
        {
          unsigned int gi = sub[cells[c].first];
          unsigned int gj = sub[cells[c].second];

          unsigned int si = gi, sj = gj;
          if (si > sj)
            std::swap(si, sj);

          unsigned int row = _cell_index.at(std::make_pair(si, sj));

          // _primal_1cells[row] may store either orientation of this
          // edge (see setup_subdivision()); flip the sign of the
          // Whitney 1-form if it disagrees with the (gi,gj) direction
          // evaluated here.
          double sign = (_primal_1cells[row].first == gi) ? 1.0 : -1.0;
          _w1[row][p] = sign * w1[c][0];
        }
      }

      if (violation <= 1e-10)
        break;
    }
  }
}


void
DEC::setup_subdivision(void)
{
  _sub_simplices.clear();
  _primal_1cells.clear();

  const libMesh::Elem& elem = *_elem;

  unsigned int dim = elem.dim();
  unsigned int nn = elem.n_nodes();

  if (dim == 0)
    return;

  if (dim == 1)
  {
    // libMesh doesn't currently assign 1D elements any edges, although
    // they logically have a single one.
    _sub_simplices.push_back({0, 1});
    _primal_1cells.push_back(std::make_pair(0u, 1u));
    return;
  }

  // Real edges of the element, in libMesh's own order *and orientation*.
  // The orientation must be kept as libMesh gives it (rather than
  // normalized to the lower node index first): it encodes the element's
  // boundary orientation (e.g. counterclockwise in 2D), which the
  // primal/dual pairing sign in get_hodge() relies on to get a
  // positive-definite Hodge star. Virtual edges introduced by the
  // subdivision (if any) are appended below, in an arbitrary but fixed
  // orientation.
  unsigned int n_real_edges = elem.n_edges();
  _primal_1cells.reserve(n_real_edges);

  // Tracks, for each undirected node pair already covered, that it has
  // been seen (regardless of which of the two orientations was stored
  // in _primal_1cells), so virtual edges are not added twice.
  std::set<std::pair<unsigned int, unsigned int>> seen;
  for (unsigned int e = 0; e < n_real_edges; ++e)
  {
    unsigned int ni = elem.local_edge_node(e, 0);
    unsigned int nj = elem.local_edge_node(e, 1);
    _primal_1cells.push_back(std::make_pair(ni, nj));

    if (ni > nj)
      std::swap(ni, nj);
    seen.insert(std::make_pair(ni, nj));
  }

  bool is_simplex = (nn == dim + 1);

  if (is_simplex)
  {
    _sub_simplices.emplace_back(nn);
    std::iota(_sub_simplices.back().begin(), _sub_simplices.back().end(), 0);
  }
  else if (elem.type() == libMesh::QUAD4)
  {
    // Split the quad into two triangles by the diagonal connecting the
    // pair of nodes with the larger sum of subtended angles.
    std::pair<unsigned int, unsigned int> diag = larger_angle_pair(
        elem.point(0), elem.point(1), elem.point(2), elem.point(3));
    unsigned int d0 = diag.first;
    unsigned int d2 = diag.second;

    _sub_simplices.push_back({d0, d2, (d2 + 1) % 4});
    _sub_simplices.push_back({d2, d0, (d0 + 1) % 4});
  }
  else if (elem.type() == libMesh::HEX8)
  {
    // Split the hexahedron into six tetrahedra sharing the main diagonal
    // between nodes 0 and 6, fanning around the hexagonal "belt" of
    // vertices connected to either endpoint by a real edge.
    static const unsigned int belt[6] = {1, 2, 3, 7, 4, 5};
    for (unsigned int k = 0; k < 6; ++k)
      _sub_simplices.push_back({0, 6, belt[k], belt[(k + 1) % 6]});
  }
  else
  {
    Messages::error("DEC: no simplicial subdivision implemented for this element type.");
    return;
  }

  // Append the virtual 1-cells required by the subdivision, i.e. the
  // sub-simplex edges that are not already real edges of the element.
  // Each virtual edge keeps the orientation it is first encountered in
  // (the choice is arbitrary, but must be fixed: reinit_forms() looks it
  // up again to tell whether a given sub-simplex evaluates it in the
  // same or the opposite direction).
  for (const std::vector<unsigned int>& sub : _sub_simplices)
    for (unsigned int a = 0; a < sub.size(); ++a)
      for (unsigned int b = a + 1; b < sub.size(); ++b)
      {
        unsigned int gi = sub[a];
        unsigned int gj = sub[b];

        unsigned int si = gi, sj = gj;
        if (si > sj)
          std::swap(si, sj);

        if (seen.insert(std::make_pair(si, sj)).second)
          _primal_1cells.push_back(std::make_pair(gi, gj));
      }
}


std::pair<unsigned int, unsigned int>
DEC::larger_angle_pair(const libMesh::Point& x1,
                        const libMesh::Point& x2,
                        const libMesh::Point& x3,
                        const libMesh::Point& x4) const
{
  // Normal from the diagonals (unbiased w.r.t. any single vertex)
  const Point N = (x3 - x1).cross(x4 - x2);

  // Vectors at x2 and x4
  const Point a = x1 - x2, b = x3 - x2;
  const Point c = x1 - x4, d = x3 - x4;

  const double P = a.cross(b) * N;   // proportional to |a||b| sin(angle2)
  const double Q = a * b;            // proportional to |a||b| cos(angle2)
  const double R = c.cross(d) * N;   // proportional to |c||d| sin(angle4)
  const double S = c * d;            // proportional to |c||d| cos(angle4)

  const double test = P * S + Q * R; // proportional to sin(angle2 + angle4)

  return (test > 0 ? std::make_pair(0u, 2u) : std::make_pair(1u, 3u));
}


void
DEC::get_hodge(libMesh::DenseMatrix<double>& hodge,
               const libMesh::RealTensor& metric)
{
  unsigned int dim = _elem->dim();
  unsigned int nn  = _elem->n_nodes();

  RealTensor R;
  R(0, 1) = -1.0;
  R(1, 0) =  1.0;
  R(2, 2) =  1.0;

  hodge.resize(_primal.size(), _primal.size());
  hodge.zero();

  bool is_simplex = (nn == (dim + 1));

  if (dim == 1)
  {
    hodge(0, 0) =  metric(0, 0) / _primal[0].norm();
  }
  else if (is_simplex || (_hodge_constr == INTERPOLATION))
  {
    if (dim == 2)
    {
      Point center(_center);

      // for QUAD, we change the center point for Hodge calculation
      // to be the midpoint of the virtual edge.
      if (_elem->type() == libMesh::QUAD4)
      {
        center = _midpoints[_primal.size() - 1];
      }


      // for a simplex, we can perform calculations in physical coordinates.
      // Also, there is no need to use a mimetic Hodge.
      for (unsigned int e = 0; e < _primal.size(); ++e)
      {
        // we use the midpoint of the dual edge segment as integration point
        Point q_point = 0.5 * (_midpoints[e] + center);
        RealGradient dual = center - _midpoints[e];

        Point tmp = _primal[e].cross(dual);

        // Evaluate the Whitney interpolation 1-forms at the integration point
        reinit_forms({q_point});

        auto &w1 = get_1forms();

        for (unsigned int i = 0; i < _primal.size(); ++i)
        {
          RealGradient w_a = w1[i][0];
          w_a = metric * w_a;

          w_a = w_a.cross(dual);

          hodge(e, i) = w_a(2);
        }
      }
    }
    else // dim == 3
    {
      // Dual faces are only defined for real edges of the element: the
      // virtual 1-cells introduced by a logical subdivision (see
      // setup_subdivision()) only serve as Whitney interpolation basis
      // functions (i.e. as columns below), not as independent Hodge rows,
      // just like they are excluded from the dual volume calculation in
      // init().
      for (unsigned int e = 0; e < _elem->n_edges(); ++e)
      {
        unsigned int ni = _elem->local_edge_node(e, 0);
        Point a = _midpoints[e] - _elem->point(ni);

        // first basis vector for surface patch
        Point v = _center - _midpoints[e];

        for (unsigned int s = 0; s < _elem->n_sides(); ++s)
        {
          if (_elem->is_edge_on_side(e, s))
          {
            Point c = _elem->side_ptr(s)->vertex_average();

            // second basis vector for surface patch
            Point w = c - _midpoints[e];

            // the cross product of the two basis vectors gives the normal vector to the surface patch
            Point n = v.cross(w);

            // check orientation of the normal vector
            if (n * a < 0)
              n *= -1;

            // we use the midpoint of the dual edge patches as integration points
            Point q_point = 1.0 / 3.0 * (_center + _midpoints[e] + c);

            // Evaluate the Whitney interpolation 1-forms at the integration point
            reinit_forms({q_point});

            auto &w1 = get_1forms();

            for (unsigned int i = 0; i < _primal.size(); ++i)
            {
              RealGradient w_a = w1[i][0];
              w_a = metric * w_a;

              hodge(e, i) += 0.5 * w_a * n;
            }
          }
        }
      }
    }
  }
  else
  {
    compute_hodge_mfd(*_elem, hodge, metric);
  }
}



void
DEC::get_incidence_pairs(std::vector<std::pair<unsigned int, unsigned int>>& inc) const
{
  inc.clear();
  inc.reserve(_incidence.m());

  for (unsigned int i = 0; i < _incidence.m(); ++i)
  {
    auto p = std::make_pair(0, 0);
    for (unsigned int j = 0; j < _incidence.n(); ++j)
    {
      if (_incidence(i, j) > 0)
        p.second = j;
      else if (_incidence(i, j) < 0)
        p.first = j;
    }
    inc.push_back(p);
  }
}


libMesh::Point
DEC::get_center(const libMesh::Elem& elem) const
{
  unsigned int dim = elem.dim();
  unsigned int nn = elem.n_nodes();

  Point center;

  if (dim == 0)
    center = elem.point(0);

  else if (dim == 1)
    center = 0.5 * (elem.point(0) + elem.point(1));

  else if (dim == 2)
  {
    if (nn == 3)
    {
      if (_dual_constr == BARYCENTRIC)
        center = elem.vertex_average();
      else
      {
        center = circumcenter(elem);
        if ((_dual_constr == MIXED) && !elem.contains_point(_center))
          center = elem.vertex_average();
      }
    }
    else if (nn == 4)
    {
      // For quadrilaterals, use intersection of diagonals
      center = diagonal_intersection(elem);
    }
  }

  else if (dim == 3)
  {
    if (nn == 4)
    {
      if (_dual_constr == BARYCENTRIC)
        center = elem.vertex_average();
      else
      {
        center = circumcenter(elem);
        if ((_dual_constr == MIXED) && !elem.contains_point(_center))
          center = elem.vertex_average();
      }
    }
    else // for now use barycenter for other 3D elements
      center = elem.vertex_average();
  }

  return center;
}




libMesh::Point
DEC::circumcenter(const libMesh::Elem& elem) const
{
  Point x_i(0.0);

  unsigned int dim = elem.dim();

  if ((dim == 2) && (elem.n_nodes() == 3))
  {
    Point a, b, c;

    a = elem.point(0);
    b = elem.point(1);
    c = elem.point(2);

    double d = 2 * (a(0) * (b(1) - c(1)) +
                    b(0) * (c(1) - a(1)) + c(0) * (a(1) - b(1)));

    x_i(0) = a.norm_sq() * (b(1) - c(1)) + b.norm_sq() * (c(1) - a(1)) +
             c.norm_sq() * (a(1) - b(1));
    x_i(1) = a.norm_sq() * (b(0) - c(0)) + b.norm_sq() * (c(0) - a(0)) +
             c.norm_sq() * (a(0) - b(0));
    x_i(1) *= -1;
    x_i /= d;

  }
  else if ((dim == 3) && (elem.n_nodes() == 4))
  {
    // tetrahedron
    // circumcenter
    Point u1(elem.point(1) - elem.point(0));
    Point u2(elem.point(2) - elem.point(0));
    Point u3(elem.point(3) - elem.point(0));

    double l1 = u1.norm_sq();
    double l2 = u2.norm_sq();
    double l3 = u3.norm_sq();

    x_i = u2.cross(u3);
    double den = 2 * u1 * x_i;

    x_i *= l1;
    x_i += l2 * u3.cross(u1) + l3 * u1.cross(u2);

    x_i /= den;

    x_i += elem.point(0);
  }
  else
    x_i = elem.vertex_average();

  return(x_i);
}


/*
 * Compute intersection of quad diagonals
 * Diagonal 1: p[0] -> p[2]
 * Diagonal 2: p[1] -> p[3]
 * Returns the intersection point xD
 */
Point
DEC::diagonal_intersection(const libMesh::Elem& elem) const
{
  /*
  const Point& p0 = elem.point(0);
  const Point& p1 = elem.point(1);
  const Point& p2 = elem.point(2);
  const Point& p3 = elem.point(3);

  // Parametrize:
  // Diagonal 1: p0 + t*(p2-p0)
  // Diagonal 2: p1 + s*(p3-p1)
  // Solve: p0 + t*(p2-p0) = p1 + s*(p3-p1)
  // => t*(p2-p0) - s*(p3-p1) = p1-p0

  Point d1 = p2 - p0;  // direction of diagonal 1
  Point d2 = p3 - p1;  // direction of diagonal 2
  Point d  = p1 - p0;  // rhs

  // 2x2 system:
  // d1.x * t - d2.x * s = d.x
  // d1.y * t - d2.y * s = d.y
  // Solve by Cramer's rule
  Real denom = d1(0)*(-d2(1)) - d1(1)*(-d2(0));
  //         = -d1(0)*d2(1) + d1(1)*d2(0)
  //         = -(d1(0)*d2(1) - d1(1)*d2(0))

  libmesh_assert_greater(std::abs(denom), 1e-14);

  Real t = (d(0)*(-d2(1)) - d(1)*(-d2(0))) / denom;
  //     = (-d(0)*d2(1) + d(1)*d2(0)) / denom

  // Intersection point
  Point xD = p0 + t * d1;

  // Verify with s (debug check)
  Real s = (d1(0)*d(1) - d1(1)*d(0)) / denom;
  Point xD_check = p1 + s * d2;
  libmesh_assert_less((xD - xD_check).norm(), 1e-10);

  return xD;
  */

  const Point& p0 = elem.point(0);
  const Point& p1 = elem.point(1);
  const Point& p2 = elem.point(2);
  const Point& p3 = elem.point(3);

  Point d1 = p2 - p0;  // direction of diagonal 1
  Point d2 = p3 - p1;  // direction of diagonal 2
  Point d  = p1 - p0;  // rhs

  // The system t*d1 - s*d2 = d has 3 equations, 2 unknowns.
  // Find the two equations with largest |determinant| for stability.
  // This corresponds to projecting onto the plane of the quad
  // by dropping the coordinate most aligned with the quad normal.

  // Quad normal (unnormalized)
  Point normal = d1.cross(d2);

  // Drop the coordinate with largest absolute normal component
  // to get the most stable 2x2 subsystem
  Real nx = std::abs(normal(0));
  Real ny = std::abs(normal(1));
  Real nz = std::abs(normal(2));

  int i, j; // indices of the two coordinates to use
  if (nx >= ny && nx >= nz)
  {
    // Normal most aligned with x: use y,z equations
    i = 1; j = 2;
  }
  else if (ny >= nx && ny >= nz)
  {
    // Normal most aligned with y: use x,z equations
    i = 0; j = 2;
  }
  else
  {
    // Normal most aligned with z: use x,y equations
    i = 0; j = 1;
  }

  // 2x2 system using coordinates i and j:
  // d1[i]*t - d2[i]*s = d[i]
  // d1[j]*t - d2[j]*s = d[j]
  //
  // Matrix A = | d1[i]  -d2[i] |
  //            | d1[j]  -d2[j] |
  Real denom = d1(i)*(-d2(j)) - (-d2(i))*d1(j);
  //         = -d1(i)*d2(j) + d2(i)*d1(j)

  if (std::abs(denom) < 1e-14)
  {
    libmesh_warning("Degenerate quad: diagonals are parallel");
    return 0.25*(p0+p1+p2+p3);
  }

  Real t = (d(i)*(-d2(j)) - (-d2(i))*d(j)) / denom;
  Real s = (d1(i)*d(j)   - d(i)*d1(j))    / denom;

  Point xD = p0 + t*d1;

  // Sanity checks
  libmesh_assert_greater(t, -1e-10);
  libmesh_assert_less(t,    1.0+1e-10);
  libmesh_assert_greater(s, -1e-10);
  libmesh_assert_less(s,    1.0+1e-10);

#ifdef DEBUG
  // Verify both parametrizations agree
  Point xD_check = p1 + s*d2;
  libmesh_assert_less((xD-xD_check).norm(), 1e-8);
#endif

  return xD;
}



void
DEC::compute_hodge_mfd(const libMesh::Elem& elem,
                       libMesh::DenseMatrix<libMesh::Real>& H,
                       const libMesh::RealTensor& metric) const
{
  const unsigned int dim    = elem.dim();
  const unsigned int n_e    = elem.n_edges();
  const unsigned int n_nodes= elem.n_nodes();
  const unsigned int n_test = dim; // linear test fields: u=x, u=y, u=z

  H.resize(n_e, n_e);
  H.zero();

  //----------------------------------------------------------------
  // Element center (barycenter)
  //----------------------------------------------------------------
  Point x_c;
  for (unsigned int i=0; i<n_nodes; ++i)
    x_c += elem.point(i);
  x_c /= n_nodes;

  //----------------------------------------------------------------
  // Edge midpoints and primal edge vectors
  // C[r,k] = (x_head - x_tail)[k] for edge r, coordinate k
  //----------------------------------------------------------------
  DenseMatrix<Real> C(n_e, n_test);
  C.zero();

  std::vector<Point> midpoints(n_e);
  for (unsigned int r=0; r<n_e; ++r)
    for (unsigned int k=0; k<n_test; ++k)
      C(r,k) = _primal[r](k);

  //----------------------------------------------------------------
  // Dual edge/face vectors and flux matrix R
  // In 2D: R[r,k] = dual_edge[r] rotated 90 degrees, component k
  // In 3D: R[r,k] = dual_face_normal[r] * area, component k
  //----------------------------------------------------------------
  DenseMatrix<Real> R(n_e, n_test);
  R.zero();

  if (dim == 2)
  {
    //--------------------------------------------------------------
    // 2D: dual edge from midpoint to center
    // flux of star(dx^k) through dual edge = (dual_edge)_perp[k]
    // star(dx) = dy => flux = (dual_edge)_y
    // star(dy) = -dx => flux = -(dual_edge)_x
    //--------------------------------------------------------------
    for (unsigned int r=0; r<n_e; ++r)
    {
      Point dual = x_c - _midpoints[r];
      // Apply metric: R[r,k] = star(e_k) . dual_edge
      // star(e_0=dx) = mu_xx*dy - mu_yx*dx (with metric)
      // More precisely: flux of mu*star(du) through dual edge
      // For test field u=x^k: du = e_k, star(e_k) is the dual
      // R[r,0]: flux of star(dx) = metric applied
      R(r,0) =  (metric(0,0)*dual(1) - metric(1,0)*dual(0));
      R(r,1) =  (metric(0,1)*dual(1) - metric(1,1)*dual(0));
    }
  }
  else // dim == 3
  {
    //--------------------------------------------------------------
    // 3D: dual face for each edge
    // Dual face = polygon connecting:
    //   midpoint m_r, adjacent face centers, element center x_c
    // Area vector = sum of triangle areas from x_c
    //--------------------------------------------------------------
    const unsigned int n_f = elem.n_faces();

    // Precompute face centers
    std::vector<Point> face_centers(n_f);
    for (unsigned int f=0; f<n_f; ++f)
    {
      auto face = elem.build_side_ptr(f);
      for (unsigned int i=0; i<face->n_nodes(); ++i)
        face_centers[f] += face->point(i);
      face_centers[f] /= face->n_nodes();
    }

    for (unsigned int r=0; r<n_e; ++r)
    {
      // Find faces adjacent to edge r
      // A face is adjacent to edge r if it contains both
      // endpoint nodes of edge r
      auto edge = elem.build_edge_ptr(r);
      dof_id_type n0 = elem.get_node_index(edge->node_ptr(0));
      dof_id_type n1 = elem.get_node_index(edge->node_ptr(1));

      std::vector<unsigned int> adj_faces;
      for (unsigned int f=0; f<n_f; ++f)
      {
        auto face = elem.build_side_ptr(f);
        bool has_n0=false, has_n1=false;
        for (unsigned int i=0; i<face->n_nodes(); ++i)
        {
          dof_id_type nf = elem.get_node_index(face->node_ptr(i));
          if (nf == n0) has_n0 = true;
          if (nf == n1) has_n1 = true;
        }
        if (has_n0 && has_n1) adj_faces.push_back(f);
      }

      // Order face centers consistently around edge direction
      if (adj_faces.size() > 2)
      {
        auto edge = elem.build_edge_ptr(r);
        Point d = (edge->point(1) - edge->point(0)).unit();

        // Reference perpendicular direction
        Point ref = face_centers[adj_faces[0]] - _midpoints[r];
        ref = ref - (ref * d) * d; // project out edge component
        Real ref_norm = ref.norm();

        if (ref_norm > 1e-14)
        {
          ref /= ref_norm;

          std::sort(adj_faces.begin(), adj_faces.end(),
                    [&](unsigned int a, unsigned int b)
                    {
                      Point va = face_centers[a] - _midpoints[r];
                      Point vb = face_centers[b] - _midpoints[r];
                      // Project out edge component
                      va = va - (va * d) * d;
                      vb = vb - (vb * d) * d;
                      // Angle relative to reference direction
                      Real angle_a = std::atan2((va.cross(ref)) * d, va * ref);
                      Real angle_b = std::atan2((vb.cross(ref)) * d, vb * ref);
                      return angle_a < angle_b;
                    });
        }
      }

      // Dual face area vector:
      // polygon: m_r -> face_centers[adj[0]] -> x_c
      //               -> face_centers[adj[1]] -> m_r (for tet, 2 faces)
      // Split into triangles from x_c:
      // tri_k: {x_c, m_r, face_centers[adj[k]]}
      // and:   {x_c, face_centers[adj[k]], m_r} -- need consistent ordering

      // Compute area vector as sum of cross products
      Point area_vec;
      const Point& mr = _midpoints[r];

      unsigned int nf = adj_faces.size();

      if (nf == 2)
      {
        // Tet: dual face is quadrilateral {mr, fc0, x_c, fc1}
        Point fc0 = face_centers[adj_faces[0]];
        Point fc1 = face_centers[adj_faces[1]];
        // Split into 2 triangles from mr:
        // {mr, fc0, x_c} and {mr, x_c, fc1}
        area_vec = 0.5*(fc0-mr).cross(x_c-mr)
                 + 0.5*(x_c-mr).cross(fc1-mr);
      }
      else if (nf == 4)
      {
        // Hex: dual face is octagon {mr,fc0,x_c,fc1,mr,...}
        // Actually: {mr, fc0, x_c} + {mr, x_c, fc1} +
        //           {mr, fc2, x_c} + {mr, x_c, fc3}
        // Need correct ordering of face centers around edge
        // For now: sum triangles from mr to consecutive face centers via x_c
        for (unsigned int k=0; k<nf; ++k)
        {
          Point fc_k = face_centers[adj_faces[k]];
          area_vec += 0.5*(fc_k-mr).cross(x_c-mr);
        }
      }
      else
      {
        // General: sum triangles {mr, fc_k, x_c}
        for (unsigned int k=0; k<nf; ++k)
        {
          Point fc_k = face_centers[adj_faces[k]];
          area_vec += 0.5*(fc_k-mr).cross(x_c-mr);
        }
      }
      // After computing area_vec, check and fix orientation
      Point edge_vec;
      auto edge_ptr = elem.build_edge_ptr(r);
      for (unsigned int k = 0; k < 3; ++k)
        edge_vec(k) = edge_ptr->point(1)(k) - edge_ptr->point(0)(k);

      // R.C should be positive: area_vec should have positive
      // dot product with edge_vec (after metric application)
      // Check without metric first:
      Real dot = area_vec * edge_vec;
      if (dot < 0)
        area_vec *= -1.0;

      // Apply metric: R[r,k] = (mu * area_vec)[k]
      // For test field u=x^k: star(du)=star(e_k) is a 2-form
      // flux through dual face = (mu * area_vec) . e_k
      for (unsigned int k=0; k<3; ++k)
      {
        Real flux = 0.0;
        for (unsigned int j=0; j<3; ++j)
          flux += metric(k,j)*area_vec(j);
        R(r,k) = flux;
      }
    }
  }

  //----------------------------------------------------------------
  // MFD formula: H = R*C^dagger + alpha*P
  // C^dagger = (C^T*C)^{-1}*C^T
  // P = I - C*C^dagger
  //----------------------------------------------------------------

  // C^T * C (n_test x n_test, small matrix)
  DenseMatrix<Real> CtC(n_test, n_test);
  CtC.zero();
  for (unsigned int i=0; i<n_test; ++i)
    for (unsigned int j=0; j<n_test; ++j)
      for (unsigned int r=0; r<n_e; ++r)
        CtC(i,j) += C(r,i)*C(r,j);

  // Invert C^T*C
  DenseMatrix<Real> CtC_inv(n_test, n_test);
  CtC_inv = CtC;

  // Use the explicit inverse:
  if (n_test == 2)
  {
    Real det = CtC(0,0)*CtC(1,1)-CtC(0,1)*CtC(1,0);
    libmesh_assert_greater(std::abs(det), 1e-14);
    CtC_inv(0,0) =  CtC(1,1)/det;
    CtC_inv(0,1) = -CtC(0,1)/det;
    CtC_inv(1,0) = -CtC(1,0)/det;
    CtC_inv(1,1) =  CtC(0,0)/det;
  }
  else // n_test == 3
  {
    Real det = CtC(0,0)*(CtC(1,1)*CtC(2,2)-CtC(1,2)*CtC(2,1))
             - CtC(0,1)*(CtC(1,0)*CtC(2,2)-CtC(1,2)*CtC(2,0))
             + CtC(0,2)*(CtC(1,0)*CtC(2,1)-CtC(1,1)*CtC(2,0));
    libmesh_assert_greater(std::abs(det), 1e-14);
    CtC_inv(0,0) = (CtC(1,1)*CtC(2,2)-CtC(1,2)*CtC(2,1))/det;
    CtC_inv(0,1) = (CtC(0,2)*CtC(2,1)-CtC(0,1)*CtC(2,2))/det;
    CtC_inv(0,2) = (CtC(0,1)*CtC(1,2)-CtC(0,2)*CtC(1,1))/det;
    CtC_inv(1,0) = (CtC(1,2)*CtC(2,0)-CtC(1,0)*CtC(2,2))/det;
    CtC_inv(1,1) = (CtC(0,0)*CtC(2,2)-CtC(0,2)*CtC(2,0))/det;
    CtC_inv(1,2) = (CtC(0,2)*CtC(1,0)-CtC(0,0)*CtC(1,2))/det;
    CtC_inv(2,0) = (CtC(1,0)*CtC(2,1)-CtC(1,1)*CtC(2,0))/det;
    CtC_inv(2,1) = (CtC(0,1)*CtC(2,0)-CtC(0,0)*CtC(2,1))/det;
    CtC_inv(2,2) = (CtC(0,0)*CtC(1,1)-CtC(0,1)*CtC(1,0))/det;
  }

  // C_dagger = CtC_inv * C^T  (n_test x n_e)
  DenseMatrix<Real> C_dag(n_test, n_e);
  C_dag.zero();
  for (unsigned int i=0; i<n_test; ++i)
    for (unsigned int r=0; r<n_e; ++r)
      for (unsigned int k=0; k<n_test; ++k)
        C_dag(i,r) += CtC_inv(i,k)*C(r,k);

  // H_c = R * C_dag  (n_e x n_e)
  DenseMatrix<Real> H_c(n_e, n_e);
  H_c.zero();
  for (unsigned int r=0; r<n_e; ++r)
    for (unsigned int s=0; s<n_e; ++s)
      for (unsigned int k=0; k<n_test; ++k)
        H_c(r,s) += R(r,k)*C_dag(k,s);

  // P = I - C * C_dag  (n_e x n_e)
  DenseMatrix<Real> P(n_e, n_e);
  P.zero();
  for (unsigned int r=0; r<n_e; ++r) P(r,r) = 1.0;
  for (unsigned int r=0; r<n_e; ++r)
    for (unsigned int s=0; s<n_e; ++s)
      for (unsigned int k=0; k<n_test; ++k)
        P(r,s) -= C(r,k)*C_dag(k,s);

  // Stabilization: alpha = tr(H_c) / n_e
  Real alpha = 0.0;
  for (unsigned int r=0; r<n_e; ++r)
    alpha += H_c(r,r);
  alpha /= n_e;
  if (alpha < 1e-14) alpha = 1.0;


  // H = H_c + alpha * P
  for (unsigned int r=0; r<n_e; ++r)
    for (unsigned int s=0; s<n_e; ++s)
      H(r,s) = H_c(r,s) + alpha*P(r,s);

}