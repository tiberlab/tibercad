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

#include <fstream>
#include <cassert>

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

  // Reinit the Whitney interpolation object whithout a point
  _whip.reinit(elem, {});

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
    ne = _whip.get_1cells().size();
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
      unsigned int ni = _whip.get_1cells()[e].first;
      unsigned int nj = _whip.get_1cells()[e].second;
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
      // for a simplex, we can perform calculations in phycial coordinates.
      // Also, there is no need to use a mimetic Hodge.
      for (unsigned int e = 0; e < _primal.size(); ++e)
      {
        // we use the midpoint of the dual edge segment as integration point
        Point q_point = 0.5 * (_midpoints[e] + _center);
        RealGradient dual = _center - _midpoints[e];

        // Reinit the Whitney interpolation object
        _whip.reinit(*_elem, {q_point});

        auto &w1 = _whip.get_1forms();

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
      for (unsigned int e = 0; e < _primal.size(); ++e)
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

            // Reinit the Whitney interpolation object
            _whip.reinit(*_elem, {q_point});

            auto &w1 = _whip.get_1forms();

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
 * Compute the local Hodge matrix H for a quadrilateral element
 * using consistency with linear fields and graph compatibility.
 * 
 * Element nodes ordered anti-clockwise:
 *   3---2
 *   |   |
 *   0---1
 * 
 * Edges ordered anti-clockwise:
 *   e0: 0->1 (bottom)
 *   e1: 1->2 (right)
 *   e2: 2->3 (top)
 *   e3: 3->0 (left)
 * 
 * Dual edges connect edge midpoints to element center.
 * Non-adjacent pairs (opposite edges): (0,2) and (1,3)
 */
void
DEC::compute_quad_hodge_mfd(const libMesh::Elem& elem,
                        libMesh::DenseMatrix<libMesh::Real>& H,
                        const libMesh::RealTensor& metric) const
{
    assert(elem.n_nodes() == 4);
    H.resize(4, 4);
    H.zero();

    // Node coordinates
    const libMesh::Point& p0 = elem.point(0);
    const libMesh::Point& p1 = elem.point(1);
    const libMesh::Point& p2 = elem.point(2);
    const libMesh::Point& p3 = elem.point(3);

    // Element center
    Point center = 0.25 * (p0 + p1 + p2 + p3);

    // Edge midpoints and dual edge vectors
    std::vector<Point> midpoints(4), dual(4), primal(4);
    midpoints[0] = 0.5*(p0+p1); primal[0] = p1-p0;
    midpoints[1] = 0.5*(p1+p2); primal[1] = p2-p1;
    midpoints[2] = 0.5*(p2+p3); primal[2] = p3-p2;
    midpoints[3] = 0.5*(p3+p0); primal[3] = p0-p3;
    for (unsigned int r = 0; r < 4; ++r)
        dual[r] = center - midpoints[r];

    // Build C matrix (4x2): cochain values for u=x and u=y
    // c_r^(1) = primal[r].x, c_r^(2) = primal[r].y
    DenseMatrix<Real> C(4, 2);
    for (unsigned int r = 0; r < 4; ++r)
    {
        C(r, 0) = primal[r](0); // d(x) cochain
        C(r, 1) = primal[r](1); // d(y) cochain
    }

    // Build R matrix (4x2): exact dual fluxes
    // For mu*star(dx): flux through dual[r] = mu applied to star(dx)
    // star(dx) = dy, so mu*star(dx) has components (mu_yx, mu_yy)
    // flux = (mu_yx)*dual[r].x + (mu_yy)*dual[r].y  -- wait
    // More carefully: star(du) . dual[r] where du = (1,0) or (0,1)
    // With metric mu: (star du)_i = mu_ij (du)_j rotated 90 degrees
    // In 2D: star(a dx + b dy) = (mu_xx*a + mu_xy*b)dy 
    //                           -(mu_yx*a + mu_yy*b)dx
    // flux through dual[r] = (mu_xx*a+mu_xy*b)*dual[r].y
    //                       -(mu_yx*a+mu_yy*b)*dual[r].x
    DenseMatrix<Real> R(4, 2);
    for (unsigned int r = 0; r < 4; ++r)
    {
        // test field u=x: du=(1,0)
        R(r, 0) = (metric(0,0)*dual[r](1) - metric(1,0)*dual[r](0));
        // test field u=y: du=(0,1)
        R(r, 1) = (metric(0,1)*dual[r](1) - metric(1,1)*dual[r](0));
    }

    // Compute C^T C (2x2)
    DenseMatrix<Real> CtC(2, 2);
    CtC.zero();
    for (unsigned int i = 0; i < 2; ++i)
        for (unsigned int j = 0; j < 2; ++j)
            for (unsigned int r = 0; r < 4; ++r)
                CtC(i,j) += C(r,i) * C(r,j);

    // Invert C^T C
    DenseMatrix<Real> CtC_inv(2, 2);
    Real det = CtC(0,0)*CtC(1,1) - CtC(0,1)*CtC(1,0);
    libmesh_assert_greater(std::abs(det), 1e-14);
    CtC_inv(0,0) =  CtC(1,1)/det;
    CtC_inv(0,1) = -CtC(0,1)/det;
    CtC_inv(1,0) = -CtC(1,0)/det;
    CtC_inv(1,1) =  CtC(0,0)/det;

    // Compute C_dag = (C^T C)^{-1} C^T  (2x4)
    DenseMatrix<Real> C_dag(2, 4);
    C_dag.zero();
    for (unsigned int i = 0; i < 2; ++i)
        for (unsigned int r = 0; r < 4; ++r)
            for (unsigned int k = 0; k < 2; ++k)
                C_dag(i,r) += CtC_inv(i,k) * C(r,k);

    // Compute consistency part: H_c = R * C_dag  (4x4)
    DenseMatrix<Real> H_c(4, 4);
    H_c.zero();
    for (unsigned int r = 0; r < 4; ++r)
        for (unsigned int s = 0; s < 4; ++s)
            for (unsigned int k = 0; k < 2; ++k)
                H_c(r,s) += R(r,k) * C_dag(k,s);

    // Compute projection P = I - C * C_dag  (4x4)
    DenseMatrix<Real> P(4, 4);
    P.zero();
    for (unsigned int r = 0; r < 4; ++r)
        P(r,r) = 1.0;
    for (unsigned int r = 0; r < 4; ++r)
        for (unsigned int s = 0; s < 4; ++s)
            for (unsigned int k = 0; k < 2; ++k)
                P(r,s) -= C(r,k) * C_dag(k,s);


    // Stabilization parameter alpha
    // A common choice is the trace of H_c divided by the rank
    Real alpha = 0.0;
    for (unsigned int r = 0; r < 4; ++r)
        alpha += H_c(r,r);
    alpha /= 4.0;

    // ensure positive
    if (alpha < 1e-14)
        alpha = 1.0;


    // H = H_c + alpha * P
    for (unsigned int r = 0; r < 4; ++r)
        for (unsigned int s = 0; s < 4; ++s)
            H(r,s) = H_c(r,s) + alpha * P(r,s);
}


/*
 * Compute the local Hodge matrix H for a quadrilateral element
 * using piecewise Whitney interpolation on subtriangles.
 *
 * Quad nodes ordered anti-clockwise: 0, 1, 2, 3
 * (using 0-based indexing throughout)
 *
 * Edges (0-based):
 *   e0: 0->1 (bottom)
 *   e1: 1->2 (right)
 *   e2: 2->3 (top)
 *   e3: 3->0 (left)
 *
 * Subtriangles:
 *   T+ = {0, 1, 2}: supports basis for e0, e1
 *   T- = {0, 1, 3}: supports basis for e0, e3
 *   T2 = {1, 2, 3}: supports basis for e1, e2
 *   T3 = {0, 2, 3}: supports basis for e2, e3
 *
 * For each edge e_r, the support is the union of the two
 * subtriangles sharing that edge. The interpolant on the
 * support is the average of the Whitney expansions on
 * each subtriangle, with diagonal cochain eliminated via
 * Stokes' theorem on each subtriangle.
 *
 * The dual edge of e_r goes from the edge midpoint m_r
 * to the diagonal intersection x_D.
 */
void
DEC::compute_quad_hodge_interp(const libMesh::Elem& elem,
                               libMesh::DenseMatrix<libMesh::Real>& H,
                               const libMesh::RealTensor& metric) const
{
  assert(elem.n_nodes() == 4);
  H.resize(4, 4);
  H.zero();

  const Point p[4] = {elem.point(0), elem.point(1),
                      elem.point(2), elem.point(3)};

  Point center = diagonal_intersection(elem);

  // matrix representation of Hodge star in R^2
  RealTensor R;
  R(0, 1) = -1.0;
  R(1, 0) =  1.0;


  for (unsigned int e = 0; e < 4; ++e)
  {
    // Dual edge vector
    Point dual_r = center - _midpoints[e];

    // integration point
    Point q_point = 0.5 * (_midpoints[e] + center);

    // Subtriangles supporting edge e
    RealGradient w1[3];
    RealGradient w2[3];

    whitney_1forms(p[e], p[(e+1)%4], p[(e+2)%4], q_point, w1);
    whitney_1forms(p[e], p[(e+1)%4], p[(e+3)%4], q_point, w2);

    RealGradient lambda1 = 0.5 *(w1[0] + w2[0] - w1[2] - w2[1]); // eliminate diagonal cochain
    RealGradient lambda2 = 0.5 *(w1[1] - w1[2]);
    RealGradient lambda3 = 0.5 *(w2[2] - w2[1]);

    double contrib_r0 = (R * metric * lambda1) * dual_r;
    double contrib_r1 = (R * metric * lambda2) * dual_r;
    double contrib_r3 = (R * metric * lambda3) * dual_r;

    H(e, e) = contrib_r0;
    H(e, (e+1)%4) += contrib_r1;
    H(e, (e+3)%4) += contrib_r3;
  }
  
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



/*
 * Compute the three Whitney 1-forms on a triangle at a given point.
 * Triangle nodes: q0, q1, q2 (in given order, anti-clockwise assumed)
 * Returns Whitney 1-forms for edges:
 *   w[0] = lambda_01 (edge q0->q1)
 *   w[1] = lambda_12 (edge q1->q2)
 *   w[2] = lambda_20 (edge q2->q0)
 * All in physical coordinates.
 * 
 * TODO: adapt for 3D triangle in 3D space (currently assumes 2D triangle in 2D space)
 */
void
DEC::whitney_1forms(const libMesh::Point& q0,
                    const libMesh::Point& q1,
                    const libMesh::Point& q2,
                    const libMesh::Point& x,
                    libMesh::RealGradient w[3]) const
{
  // Signed area via cross product
  // area = 0.5 * (q1-q0) x (q2-q0)
  Real area2 = (q1(0)-q0(0))*(q2(1)-q0(1))
             - (q1(1)-q0(1))*(q2(0)-q0(0));

  libmesh_assert_greater(std::abs(area2), 1e-14);

  Real inv2A = 1.0 / area2;

  // Gradients of barycentric coordinates (constant on triangle)
  // grad lambda_i = (1/2A) * perp(opposite edge)
  RealGradient g[3];
  g[0](0) = (q1(1)-q2(1)) * inv2A;
  g[0](1) = (q2(0)-q1(0)) * inv2A;

  g[1](0) = (q2(1)-q0(1)) * inv2A;
  g[1](1) = (q0(0)-q2(0)) * inv2A;

  g[2](0) = (q0(1)-q1(1)) * inv2A;
  g[2](1) = (q1(0)-q0(0)) * inv2A;

  // Barycentric coordinates of x
  Real lam[3];
  lam[0] = ((q1(1)-q2(1))*(x(0)-q2(0))
           + (q2(0)-q1(0))*(x(1)-q2(1))) * inv2A;
  lam[1] = ((q2(1)-q0(1))*(x(0)-q2(0))
           + (q0(0)-q2(0))*(x(1)-q2(1))) * inv2A;
  lam[2] = 1.0 - lam[0] - lam[1];

  // Whitney 1-forms: lambda_ij = lam_i * grad_j - lam_j * grad_i
  // w[0] = lambda_01: edge q0->q1
  w[0](0) = lam[0]*g[1](0) - lam[1]*g[0](0);
  w[0](1) = lam[0]*g[1](1) - lam[1]*g[0](1);

  // w[1] = lambda_12: edge q1->q2
  w[1](0) = lam[1]*g[2](0) - lam[2]*g[1](0);
  w[1](1) = lam[1]*g[2](1) - lam[2]*g[1](1);

  // w[2] = lambda_20: edge q2->q0
  w[2](0) = lam[2]*g[0](0) - lam[0]*g[2](0);
  w[2](1) = lam[2]*g[0](1) - lam[0]*g[2](1);
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