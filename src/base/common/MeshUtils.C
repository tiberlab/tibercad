// $Id$


#include "MeshUtils.h"
#include "HashMap.h"
#include "HashSet.h"
#include "RuntimeException.h"

#include "mesh.h"
#include "elem.h"
#include "mesh_tools.h"
#include "mesh_base.h"

#include <cassert>
#include <list>


using namespace std;


void
MeshUtils::get_subdomain_ids(MeshBase& mesh, std::set<ID>& subdomain_ids)
{
  subdomain_ids.clear();

  MeshBase::element_iterator it = mesh.local_elements_begin();
  const MeshBase::element_iterator end = mesh.local_elements_end();

  for ( ; it != end; ++it)
  {
    Elem* elem = *it;

    ID id = static_cast<ID>(elem->subdomain_id());
    subdomain_ids.insert(id);
  }
}


bool MeshUtils::may_belong_to_element(const Elem* element, const Point& point)
{
  Point vertex(element->point(0));
  double min_x = vertex(0);
  double min_y = vertex(1);
  double min_z = vertex(2);
  double max_x = min_x;
  double max_y = min_y;
  double max_z = min_z;

  const size_t n = element->n_nodes();
  for (size_t i = 1 ; i < n ; i++)
  {
    vertex = element->point(i);
    double x = vertex(0);
    double y = vertex(1);
    double z = vertex(2);

    if (min_x > x) min_x = x;
    if (min_y > y) min_y = y;
    if (min_z > z) min_z = z;

    if (max_x < x) max_x = x;
    if (max_y < y) max_y = y;
    if (max_z < z) max_z = z;

  }

  bool result = true;
  if ((point(0) > max_x) ||  (point(0) < min_x) ||
      (point(1) > max_y) ||  (point(1) < min_y) ||
      (point(2) > max_z) ||  (point(2) < min_z))
  {
    result = false;
  }

  return result;
}




Point
MeshUtils::get_outer_normal(const Elem* elem, int side)
{
  assert(elem != NULL);

  Point normal;

  AutoPtr<DofObject> sobj = elem->side(side);
  const Elem* side_el = dynamic_cast<Elem*>(sobj.get());
  const Point& centroid = elem->centroid();

  switch (elem->dim())
  {
    case 0:
      break;

    case 1:
      // side is a node
      normal = *side_el->get_node(0) - centroid;
      break;

    case 2:
    {
      // side should always be an Edge2
      // normal direction is: p0 + t*(p1 - p0) - centroid, where
      // t gives the intersection between the side and the perpendicular
      // through the element centroid
      Point p10((*side_el->get_node(1) - *side_el->get_node(0)).unit());
      Point p03(*side_el->get_node(0) - centroid);
      double t = p03 * p10;
      normal = p03 - t * p10;
      break;
    }

    case 3:
    {
      // TODO
      throw RuntimeException("MeshUtils::get_outer_normal() is not implemented yet for 3D");
      break;
    }

  }

  return normal.unit();
}



AutoPtr<MeshBase>
MeshUtils::create_boundary_mesh(const MeshBase& mesh)
{

  MeshBase *bdmesh = new SerialMesh(mesh.mesh_dimension() - 1);

  using namespace std;

  // extract region boundary
  unsigned int node_ctr = 0;
  HashMap<unsigned int, unsigned int>::Type node_id_map;
  set<pair<const Elem*, const Elem*>> sides_added;


  MeshBase::const_element_iterator it = mesh.level_elements_begin(0);
  const MeshBase::const_element_iterator end = mesh.level_elements_end(0);
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    ID subdomain = elem->subdomain_id();

    for (int i = 0; i < elem->n_sides(); ++i)
    {
      const Elem* nb = elem->neighbor(i);
      pair<const Elem*, const Elem*> el_pair(elem, nb);
      if (nb > elem) { el_pair.first = nb; el_pair.second = elem; }

      if ((nb == NULL) || ((nb->subdomain_id() != subdomain) && !sides_added.count(el_pair)))
      {
        sides_added.insert(el_pair);

        DofObject* sobj = elem->side(i).release();
        Elem* side_el = dynamic_cast<Elem*>(sobj);

        HashMap<unsigned int, unsigned int>::Type::iterator mit;

        for (unsigned int n = 0; n < side_el->n_nodes(); ++n)
        {
          mit = node_id_map.find(side_el->node(n));

          unsigned int id = side_el->get_node(n)->id();
          if (mit == node_id_map.end())
          {
            node_id_map[side_el->node(n)] = node_ctr;
            const Node* p = side_el->get_node(n);
            bdmesh->add_point(*p, p->id(), 0);
          }

          side_el->set_node(n) = bdmesh->node_ptr(id);
        }

        side_el->processor_id() = 0;
        bdmesh->add_elem(side_el);

      }
    }
  }



  //bdmesh->prepare_for_use();

  return AutoPtr<MeshBase>(bdmesh);
}



const Elem*
MeshUtils::search_element(const MeshBase* mesh, const Point& point)
{
  return GridMapper::get_mapper(mesh).get_element(point);
}



MeshUtils::GridMapper::GridMapper(const MeshBase* mesh)
: _mesh(mesh)
{
  setup();
}


void
MeshUtils::GridMapper::setup(void)
{
  MeshTools::BoundingBox bb(MeshTools::bounding_box(*_mesh));
  _tensor_grid.setup(bb.min(), bb.max(), 50, 50, 50);
}



MeshUtils::GridMapper&
MeshUtils::GridMapper::get_mapper(const MeshBase* mesh)
{
  std::map<const MeshBase*, GridMapper>::iterator it(_mappers.find(mesh));
  if (it == _mappers.end())
  {
    it = (_mappers.insert(make_pair(mesh, GridMapper(mesh)))).first;
  }

  return it->second;
}


const Elem*
MeshUtils::GridMapper::get_element(const Point& point) const
{
  const Elem* el = NULL;

  int tgrid_el = _tensor_grid.find_element(point);
  // we can assume that _elem_list is assembled when getting here
  vector<const Elem*>& list = _elem_list[tgrid_el];

  for (int i = 0; i < list.size(); ++i)
  {
    const Elem* elem = list[i];
    if (MeshUtils::may_belong_to_element(elem, point))
    {
      if (elem->contains_point(point))
      {
        el = elem;
        break;
      }
    }
  }

  return el;
}



