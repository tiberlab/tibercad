#include "DECdual.h"
#include "PolyTope.h"

#include "libmesh/elem.h"

#include <queue>
#include <set>

using libMesh::Elem;
using libMesh::Point;

DECdual::DECdual(const libMesh::Elem* elem, unsigned int node)
{
  build_dual_element(elem, node);
}

void
DECdual::build_dual_element(const libMesh::Elem* elem, unsigned int node)
{
  _center_id = elem->node_id(node);
  _star.clear();

  std::queue<const libMesh::Elem*> elems;
  std::set<const libMesh::Elem*> el_done;

  elems.push(elem);

  while (!elems.empty())
  {
    const libMesh::Elem* nextelem = elems.front();
    unsigned int node_local_id = nextelem->local_node(_center_id);

    for (unsigned int s = 0; s < nextelem->n_sides(); ++s)
    {
      if (nextelem->is_node_on_side(node_local_id, s))
      {
        const libMesh::Elem* neigh = nextelem->neighbor_ptr(s);

        if (!el_done.count(neigh))
        {
          _star.push_back(OneChain2D(nextelem, s));

          if (neigh != nullptr)
          {
            elems.push(neigh);
          }
        }
        
      }
    }
    el_done.insert(nextelem);
    elems.pop();
  }
}


void
DECdual::get_polytope(PolyTope &poly) const
{
  poly.clear();
  poly.reserve(_star.size());

  // if there is no polytope, we return immediately
  if (_star.size() < 3) return;

  // we have to add nodes in counter-clock wise order,
  // but we can use the fact that by construction the
  // entries in _star are almost in this order, since
  // 2D elements have exactly two sides touching the
  // center node.

  const OneChain2D& c = _star[0];
  Point p0(c.get_elem()->vertex_average());
  poly.add_node(p0);

  const Elem *neigh = _star[1].get_elem()->neighbor_ptr(_star[1].get_side());
  Point p1(_star[1].get_elem()->vertex_average());

  // check orientation
  Point s1 = p0 - c.get_elem()->point(c.get_elem()->local_node(_center_id));
  Point d1 = p1 - p0;
  double orientation = s1.cross(d1)(2);
std::cerr <<"orientation: " << orientation << "\n";
  if (orientation > 0)
  {
    poly.add_node(p1);

    for (unsigned int i = 3; i < _star.size(); ++i)
    {
      poly.add_node(_star[i].get_elem()->vertex_average());
    }

    if (_star[2].get_elem() != nullptr)
    {
      const Elem* neigh = _star[1].get_elem()->neighbor_ptr(_star[1].get_side());
      poly.add_node(neigh->vertex_average());
    }
  }
  else
  {
    if (_star[2].get_elem() != nullptr)
    {
      const Elem* neigh = _star[1].get_elem()->neighbor_ptr(_star[1].get_side());
      poly.add_node(neigh->vertex_average());
    }

    for (unsigned int i = _star.size() - 1; i >= 3; --i)
    {
      poly.add_node(_star[i].get_elem()->vertex_average());
    }

    poly.add_node(p1);
  }

for (unsigned int i = 0; i < poly.get_nodes().size(); ++i)
{
  std::cerr << poly.get_nodes().at(i) - 
      c.get_elem()->point(c.get_elem()->local_node(_center_id)) << "\n";
}
std::cerr << "\n";
}