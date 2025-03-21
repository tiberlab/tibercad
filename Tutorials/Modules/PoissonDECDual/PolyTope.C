#include "PolyTope.h"

void
PolyTope::clear(void)
{
  _nodes.clear();
}


void
PolyTope::reserve(unsigned int n)
{
  _nodes.reserve(n);
}

void
PolyTope::add_node(const libMesh::Point& point)
{
  _nodes.push_back(point);
}

const std::vector<libMesh::Point>&
PolyTope::get_nodes(void) const
{
  return(_nodes);
}

double
PolyTope::get_interpolation(const libMesh::Point &x,
                            libMesh::RealGradient &grad) const
{
  double w = 0;
  grad = 0;

  return(w);
}