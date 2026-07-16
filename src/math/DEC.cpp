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


#include "libmesh/elem.h"

using namespace std;
using namespace libMesh;


void
DEC::reinit(const libMesh::Elem& elem, DEC::DualConstruction dual_constr)
{
  _elem = &elem;

  unsigned int dim = elem.dim();
  unsigned int ne = (dim == 1) ? 1 : elem.n_edges();
  unsigned int nn = elem.n_nodes();
  
  _primal.resize(ne);
  _midpoints.resize(ne);
  _incidence.resize(ne, nn);
  _dual_volumes.resize(nn, 0.0);

  // circumcenter and thus Voronoi-construction works only
  // for triangles, otherwise we fall back to barycentric Hodge
  if ((nn != 3) || (dual_constr == BARYCENTRIC))
    _center = elem.vertex_average();
  else
  {
    _center = circumcenter(elem);
    if (!elem.contains_point(_center) && (dual_constr == MIXED))
      _center = elem.vertex_average();
  }

  if (dim == 1)
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

  hodge.resize(_primal.size(), _primal.size());

  if (dim == 1)
  {
    hodge(0, 0) =  metric(0, 0) / _primal[0].norm();
  }
  else if (dim == 2)
  {
    for (unsigned int e = 0; e < _primal.size(); ++e)
    {
      // we use the midpoint of the dual edge segment as integration point
      Point q_point = 0.5 * (_midpoints[e] + _center);
      RealGradient dual = _center - _midpoints[e];

      // Reinit the Whitney interpolation object
      whip.reinit(*_elem, {q_point});

      auto& w1 = _whip.get_1forms();

      for (unsigned int i = 0; i < _primal.size(); ++i)
      {
        RealGradient w_a = w1[i][0];
        w_a = metric * w_a;

        w_a = w_a.cross(dual);

        hodge(e, i) = w_a(2);
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

          auto& w1 = _whip.get_1forms();

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




libMesh::Point
DEC::circumcenter(const libMesh::Elem& elem, int s) const
{
  Point x_i(0.0);

  unsigned int dim = elem.dim();

  // (It seems the centroid works better for quadrangles)
  //if (dim == 2)
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

