// $Id$

#include "VoronoiCell.h"
#include "InitFailedException.h"

#include "libmesh/elem.h"



VoronoiCell::VoronoiCell(const libMesh::Elem* elem, double scaling)
  : _elem(elem),
    _scaling(scaling)
{
  if (_elem == nullptr)
    throw InitFailedException("Cannot pass nullptr to VoronoiCell.");
}


VoronoiCell::VoronoiCell(const libMesh::Elem& elem, double scaling)
  : _elem(&elem),
    _scaling(scaling)
{
  if (_elem == nullptr)
    throw InitFailedException("Cannot pass nullptr to VoronoiCell.");
}


VoronoiCell::~VoronoiCell(void)
{

}


void
VoronoiCell::calculate(void)
{
  clear();

  switch (_elem->type())
  {
    case libMesh::NODEELEM:
    {
      _volumes.resize(1, 1.0);
      break;
    }

    case libMesh::EDGE2:
    {
      _edges.push_back(std::make_pair(0, 1));
      _edge_areas.push_back(1.0);
      _edge_lengths.push_back(_elem->volume() * _scaling);
      double vol = 0.5 * _elem->volume() * _scaling;
      _volumes.resize(2, vol);

      break;
    }

    case libMesh::TRI3:
    {
      _edges.push_back(std::make_pair(0, 1));
      _edges.push_back(std::make_pair(1, 2));
      _edges.push_back(std::make_pair(2, 0));

      libMesh::Point d(_elem->point(1) - _elem->point(0));
      double dx12 = d(0);
      double dy12 = d(1);
      _edge_lengths.push_back(d.norm() * _scaling);
      d = _elem->point(2) - _elem->point(1);
      double dx23 = d(0);
      double dy23 = d(1);
      _edge_lengths.push_back(d.norm() * _scaling);
      d = _elem->point(0) - _elem->point(2);
      double dx31 = d(0);
      double dy31 = d(1);
      _edge_lengths.push_back(d.norm() * _scaling);


      // node 0:
      //
      // | -dy12  -dy13 | |a|    1|dx23|
      // |              | | | =  -|    |
      // | dx12    dx13 | |b|    2|dy23|
      //
      // Voronoi face length:
      //  A12 = a*sqrt(dy12^2 + dx12^2)

      auto p = solve_2x2(-dy12, dy31, dx12, -dx31, 0.5*dx23, 0.5*dy23);
      double a = p.first;
      double c = p.second;

      // node 1:
      p = solve_2x2(-dy23, dy12, dx23, -dx12, 0.5*dx31, 0.5*dy31);
      a = 0.5 * (a + p.second);
      double b = p.first;

      // node 2:
      p = solve_2x2(-dy31, -dy23, dx31, dx23, 0.5*dx12, 0.5*dy12);
      c = 0.5 * (c + p.first);
      b = 0.5 * (b + p.second);

      _edge_areas.push_back(a * _edge_lengths[0]);
      _edge_areas.push_back(b * _edge_lengths[1]);
      _edge_areas.push_back(c * _edge_lengths[2]);

      _volumes.resize(3);
      _volumes[0] = 0.25 * (_edge_areas[0] * _edge_lengths[0] +
                            _edge_areas[2] * _edge_lengths[2]);
      _volumes[1] = 0.25 * (_edge_areas[0] * _edge_lengths[0] +
                            _edge_areas[1] * _edge_lengths[1]);
      _volumes[2] = 0.25 * (_edge_areas[1] * _edge_lengths[1] +
                            _edge_areas[2] * _edge_lengths[2]);

      break;
    }

    default:
      throw InitFailedException("Trying to calcuate Voronoi cell "
          "on unsupported element type.");
      break;
  }

}

void
VoronoiCell::clear(void)
{
  _edges.clear();
  _edge_areas.clear();
  _edge_lengths.clear();
  _volumes.clear();
  //_nodes_to_edges.clear();
}

std::pair<double, double>
VoronoiCell::solve_2x2(double a11, double a12,
                       double a21, double a22,
                       double b1, double b2) const
{
  double det = a11*a22 - a12*a21;
  double x1 =  a22 * b1 - a12 * b2;
  double x2 = -a21 * b1 + a11 * b2;

  return(std::make_pair(x1/det, x2/det));
}

