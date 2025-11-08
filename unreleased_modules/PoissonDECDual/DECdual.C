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

  // a reasonable guess for triangular grids
  _star.reserve(6);

  std::queue<const libMesh::Elem*> elems;

  elems.push(elem);

  std::set<const libMesh::Elem*> el_done;
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
          if (el_done.count(nextelem) && (neigh != nullptr))
          {
            for (unsigned int ns = 0; ns < neigh->n_sides(); ++ns)
            {
              if (neigh->neighbor_ptr(ns) == nextelem)
              {
                _star.push_back(OneChain2D(neigh, ns));
                el_done.insert(neigh);
              }
            }
          }
          else
          {
            _star.push_back(OneChain2D(nextelem, s));
            el_done.insert(nextelem);
          }

          if (neigh != nullptr)
          {
            elems.push(neigh);
          }
        }
      }
    }
    elems.pop();
  }

  _star.shrink_to_fit();
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

  for (unsigned int i = 0; i < _star.size(); ++i)
  {
    poly.add_node(_star[i].get_elem()->vertex_average());
  }
}