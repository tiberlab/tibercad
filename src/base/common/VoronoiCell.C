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
      break;

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



