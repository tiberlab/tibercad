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
         DEC::DualConstruction dual_constr)
: _elem(&elem),
  _dual_constr(dual_constr)
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

  unsigned int dim = elem.dim();
  unsigned int ne = (dim == 1) ? 1 : elem.n_edges();
  unsigned int nn = elem.n_nodes();
  
  _primal.resize(ne);
  _midpoints.resize(ne);
  _incidence.resize(ne, nn);
  _dual_volumes.resize(nn, 0.0);


  // circumcenter and thus Voronoi-construction works only
  // for triangles, otherwise we fall back to barycentric Hodge
  if ((nn != 3) || (_dual_constr == BARYCENTRIC))
    _center = elem.vertex_average();
  else
  {
    _center = circumcenter(elem);
    if ((_dual_constr == MIXED) && !elem.contains_point(_center))
      _center = elem.vertex_average();
  }

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
      unsigned int ni = elem.local_edge_node(e, 0);
      unsigned int nj = elem.local_edge_node(e, 1);
      _primal[e] = (elem.point(nj) - elem.point(ni));
      _midpoints[e] = 0.5 * (elem.point(ni) + elem.point(nj));

      _incidence(e, ni) = -1;
      _incidence(e, nj) =  1;
    }

    for (unsigned int e = 0; e < ne; ++e)
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

  // Reinit the Whitney interpolation object on the element center
  _whip.reinit(elem, {_center});
}



void
DEC::get_hodge(libMesh::DenseMatrix<double>& hodge,
        const libMesh::RealTensor& metric) const
{
  WhitneyInterpolation whip;

  unsigned int dim = _elem->dim();

  RealTensor R;
  R(0, 1) = -1.0;
  R(1, 0) =  1.0;
  R(2, 2) =  1.0;

  hodge.resize(_primal.size(), _primal.size());
  hodge.zero();

  if (dim == 1)
  {
    hodge(0, 0) =  metric(0, 0) / _primal[0].norm();
  }
  else if (dim == 2)
  {
    if (_elem->n_nodes() == 3)
    {
      // for a simplex, we can perform calculations in phycial coordinates
      for (unsigned int e = 0; e < _primal.size(); ++e)
      {
        // we use the midpoint of the dual edge segment as integration point
        // Point q_point = 0.5 * (_midpoints[e] + _center);
        Point q_point = _midpoints[e];
        RealGradient dual = _center - _midpoints[e];

        // Reinit the Whitney interpolation object
        whip.reinit(*_elem, {q_point});

        auto &w1 = whip.get_1forms();

        for (unsigned int i = 0; i < _primal.size(); ++i)
        {
          RealGradient w_a = w1[i][0];
          w_a = metric * w_a;

          w_a = w_a.cross(dual);

          hodge(e, i) = w_a(2);
        }
      }
    }
    else if (_elem->n_nodes() == 4)
    {
      compute_quad_hodge(*_elem, hodge, metric);
    }
    else
    {
      // code with explicit integration and pull back to reference
      // element, works for simplices only

      unique_ptr<libMesh::FEBase> fe = libMesh::FEBase::build(dim, libMesh::FEType(1, libMesh::LAGRANGE));
      const vector<vector<libMesh::Real>> &phi = fe->get_phi();
      const vector<vector<libMesh::RealGradient>> &dphi = fe->get_dphi();

      // Reference coordinates of dual edge endpoints
      Point ref_center = FEMap::inverse_map(dim, _elem, _center);

      for (unsigned int e = 0; e < _primal.size(); ++e)
      {
        Point ref_mid = FEMap::inverse_map(dim, _elem, _midpoints[e]);
        
        // Reference dual edge vector (straight in reference coords)
        Point ref_dual = ref_center - ref_mid;
        Point dual = _center - _midpoints[e];

        
        // Trapezoidal integration along reference dual edge
        unsigned int nq = 20;
        for (unsigned int k = 0; k <= nq; ++k)
        {

          double s = k * 1.0;
          double wq = (k == 0 || k == nq) ? 0.5 : 1.0;
          if (nq == 0)
            wq = 1.0;
          else
          {
            wq /= nq;
            s  /= nq;
          }

          // Quadrature point in reference coordinates
          Point ref_qp = ref_mid + s * ref_dual;
          Point phys_qp = FEMap::map(dim, _elem, ref_qp);

          // Get Jacobian J_T at this reference point
          // by reinitializing FEMap at ref_qp
          std::vector<Point> ref_qp_vec = {ref_qp};
          fe->reinit(_elem, &ref_qp_vec);

          auto fe_map = fe->get_fe_map();

          // J_T columns from FEMap
          // dxyzdxi = d(physical)/d(xi), dxyzdeta = d(physical)/d(eta)
          const auto &dxyzdxi = fe_map.get_dxyzdxi();
          const auto &dxyzdeta = fe_map.get_dxyzdeta();

          // Build J_T as 2x2 matrix at this quadrature point
          // J_T = [ dxyzdxi[0](x,y) | dxyzdeta[0](x,y) ]
          RealTensor J;
          J(0, 0) = dxyzdxi[0](0);
          J(0, 1) = dxyzdeta[0](0);
          J(1, 0) = dxyzdxi[0](1);
          J(1, 1) = dxyzdeta[0](1);
          J(2, 2) = 1.0;

          // Pulled-back metric g = J_T^T * J_T in reference coords
          RealTensor g = J.transpose() * J;

          // Reinit Whitney forms at reference quadrature point
          // whip must support reinit at reference coordinates
          whip.reinit(*_elem, ref_qp_vec, true);
          auto &w1 = whip.get_1forms();


          // Reference dual edge cross product direction
          // in 2D: a x b = a_x*b_y - a_y*b_x (z-component)
          for (unsigned int i = 0; i < _primal.size(); ++i)
          {
            /*
            // w1[i][0] is the Whitney 1-form in reference coords
            RealGradient w_a = w1[i][0];

            // Apply pulled-back metric and physical tensor mu
            // Combined: g^{-1} * mu or mu * g depending on convention
            w_a = g * w_a;      // pulled-back metric
            w_a = metric * w_a; // physical anisotropy tensor

            // Cross product with reference dual edge (z-component in 2D)
            double contrib = w_a(0) * ref_dual(1) - w_a(1) * ref_dual(0);
            hodge(e, i) += wq * contrib;
            */

            // Physical gradient from whip: J_T^{-T} * grad_xi N_i
            RealGradient grad_phys = w1[i][0]; // physical grad N_i
            grad_phys = R * metric * grad_phys; 

            // Recover reference gradient: J_T^T * grad_phys = grad_xi N_i
            RealGradient grad_ref = J.transpose() * grad_phys;
            
            // Now construct reference Whitney form contribution
            // w_a = mu * grad_ref (in reference coords)
            //RealGradient w_a = metric * g.inverse() * grad_ref;
            //RealGradient w_a = metric * grad_phys;

            Point dual_pushed = J * ref_dual; // push reference dual edge to physical space

            // Cross with reference dual edge (2D z-component)
            //double contrib = w_a(0) * ref_dual(1) - w_a(1) * ref_dual(0);
            //double contrib = w_a(0) * dual_pushed(1) - w_a(1) * dual_pushed(0);
            //double contrib = grad_ref * ref_dual;
            double contrib = grad_phys * dual;

            hodge(e, i) += wq * contrib;
            
          }
        }
      }
    }
  }
  else if (dim == 3)
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
          Point q_point = 1.0/3.0 * (_center + _midpoints[e] + c);

          // Reinit the Whitney interpolation object
          whip.reinit(*_elem, {q_point});

          auto& w1 = whip.get_1forms();

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
DEC::circumcenter(const libMesh::Elem& elem, int s) const
{
  Point x_i(0.0);

  unsigned int dim = elem.dim();

  if ((dim == 2) && (elem.n_nodes() == 3))
  {
    Point a, b, c;
    
    /*
    if ((s >= 0) && (elem.n_nodes() > 3))
    {
      auto side = elem.side_ptr(s);
      a = side.point(0);
      b = side.point(1);
      Point v1(b - a);

      // look for adjacent side that makes smallest angle
      unsigned int ns = elem.n_sides();

      auto s2 = elem.side_ptr((s+1)%ns);
      Point v2(s2.point(1) - s2.point(0));

      auto s3 = elem.side_ptr((s-1)%ns);
      Point v3(s3.point(0) - s3.point(1));

      double cosa = -(v1 * v2) / (v1.norm() * v2.norm());
      double cosb =  (v1 * v3) / (v1.norm() * v3.norm());

      if (cosa > cosb)
        c = b + v2;
      else
        c = a + v3;
      
    }*/
    //else
    {
      a = elem.point(0);
      b = elem.point(1);
      c = elem.point(2);
    }

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

    if (!elem.contains_point(x_i))
    {
      x_i = elem.vertex_average();
    }
  }
  else
    x_i = elem.vertex_average();

  return(x_i);
}


/**
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
DEC::compute_quad_hodge(const libMesh::Elem& elem,
                        libMesh::DenseMatrix<libMesh::Real>& H,
                        const libMesh::RealTensor& mu) const
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
        R(r, 0) = (mu(0,0)*dual[r](1) - mu(1,0)*dual[r](0));
        // test field u=y: du=(0,1)
        R(r, 1) = (mu(0,1)*dual[r](1) - mu(1,1)*dual[r](0));
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

