// $Id$


#include "MeshUtils.h"
#include "Utils.h"
#include "Messages.h"
#include "HashMap.h"
#include "HashSet.h"
#include "RuntimeException.h"
#include "TiberCad.h"

#include "mesh.h"
#include "elem.h"
#include "mesh_tools.h"
#include "mesh_base.h"


#include <cassert>
#include <list>
#include <algorithm>

#include "libMeshDefs.h"

using namespace std;

USELIBMESHTYPE(UniquePtr);
USELIBMESHTYPE(DofObject);



void
MeshUtils::get_subdomain_ids(libMesh::MeshBase& mesh, std::set<ID>& subdomain_ids)
{
  subdomain_ids.clear();

  libMesh::MeshBase::element_iterator it = mesh.local_elements_begin();
  const libMesh::MeshBase::element_iterator end = mesh.local_elements_end();

  for ( ; it != end; ++it)
  {
    libMesh::Elem* elem = *it;

    ID id = static_cast<ID>(elem->subdomain_id());
    subdomain_ids.insert(id);
  }
}


bool MeshUtils::may_belong_to_element(const libMesh::Elem* element, libMesh::Point& point)
{
  libMesh::Point vertex(element->point(0));
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




libMesh::Point
MeshUtils::get_outer_normal(const libMesh::Elem* elem, int side)
{
  assert(elem != NULL);

  libMesh::Point normal;



    libMesh::UniquePtr<libMesh::Elem> side_el(elem->side(side));
  //const libMesh::Elem* side_el = dynamic_cast<libMesh::Elem*>(sobj.get());
  const libMesh::Point& centroid = elem->centroid();

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
      libMesh::Point p10((*side_el->get_node(1) - *side_el->get_node(0)).unit());
      libMesh::Point p03(*side_el->get_node(0) - centroid);
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



libMesh::UniquePtr<libMesh::MeshBase>
MeshUtils::create_boundary_mesh(const libMesh::MeshBase& mesh)
{
  libMesh::MeshBase* bdmesh = new libMesh::SerialMesh(TiberCad::get_mpi_comm() , mesh.mesh_dimension());

  using namespace std;

  // extract region boundary
  unsigned int node_ctr = 0;
  HashMap<unsigned int, unsigned int>::Type node_id_map;
  set<pair<const libMesh::Elem*, const libMesh::Elem*>> sides_added;


  libMesh::MeshBase::const_element_iterator it = mesh.level_elements_begin(0);
  const libMesh::MeshBase::const_element_iterator end = mesh.level_elements_end(0);
  for ( ; it != end; ++it)
  {
    const libMesh::Elem* elem = *it;
    ID subdomain = elem->subdomain_id();

    for (int i = 0; i < elem->n_sides(); ++i)
    {
      const libMesh::Elem* nb = elem->neighbor(i);
      pair<const libMesh::Elem*, const libMesh::Elem*> el_pair(elem, nb);
      if (nb > elem) { el_pair.first = nb; el_pair.second = elem; }

      if ((nb == NULL) || ((nb->subdomain_id() != subdomain) && !sides_added.count(el_pair)))
      {
        sides_added.insert(el_pair);

        libMesh::Elem* side_el = elem->side(i).release();
        //libMesh::Elem* side_el = dynamic_cast<libMesh::Elem*>(sobj);

        HashMap<unsigned int, unsigned int>::Type::iterator mit;

        for (unsigned int n = 0; n < side_el->n_nodes(); ++n)
        {
          mit = node_id_map.find(side_el->node(n));

          unsigned int id = side_el->get_node(n)->id();
          if (mit == node_id_map.end())
          {
            node_id_map[side_el->node(n)] = node_ctr;
            const libMesh::Node* p = side_el->get_node(n);
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

  return libMesh::UniquePtr<libMesh::MeshBase>(bdmesh);
}



const libMesh::Elem*
MeshUtils::search_element(const libMesh::MeshBase* mesh, const libMesh::Point& point)
{
  return GridMapper::get_mapper(mesh).get_element(point);
}



multimap<const libMesh::MeshBase*, MeshUtils::GridMapper*>
MeshUtils::GridMapper::_mappers;



MeshUtils::GridMapper::GridMapper(const libMesh::MeshBase* mesh, const set<ID>& regions)
: _mesh(mesh),
  _regids(regions)
{
  setup();
}


MeshUtils::GridMapper::~GridMapper(void)
{
  map<const libMesh::MeshBase*, GridMapper*>::iterator it(_mappers.begin());
  map<const libMesh::MeshBase*, GridMapper*>::iterator end(_mappers.end());
  for ( ; it != end; ++it)
  {
    delete it->second;
  }
}



void
MeshUtils::GridMapper::setup(void)
{
  libMesh::MeshTools::BoundingBox bb(libMesh::MeshTools::bounding_box(*_mesh));

  int nx = 50, ny = 1, nz = 1;
  switch (_mesh->mesh_dimension())
  {
    case 3:
      nz = 50;

    case 2:
      ny = 50;
      break;
  }

  _tensor_grid.setup(bb.min(), bb.max(), nx, ny, nz);

  _elem_list.resize(_tensor_grid.num_elements());

  libMesh::MeshBase::const_element_iterator it = _mesh->local_elements_begin();
  const libMesh::MeshBase::const_element_iterator end = _mesh->local_elements_end();

  Messages::newline();
  Messages::info("Setup of grid mapper: ", false);
  Utils::Timer timer;

  for ( ; it != end; ++it)
  {
    const libMesh::Elem* elem = *it;

    // only if the subdomain ID is in the required subset we proceed
    if (!_regids.empty() && !_regids.count(elem->subdomain_id()))
      continue;

    // we decide whether an element touches a certain tensor grid element
    // by looking at the bounding box

    // get the bounding box
    libMesh::Point p0(elem->point(0));
    libMesh::Point p1(p0);
    for (unsigned int n = 1; n < elem->n_nodes(); n++)
    {
      const libMesh::Point& pn = elem->point(n);
      if (pn(0) < p0(0)) p0(0) = pn(0);
      if (pn(1) < p0(1)) p0(1) = pn(1);
      if (pn(2) < p0(2)) p0(2) = pn(2);
      if (pn(0) > p1(0)) p1(0) = pn(0);
      if (pn(1) > p1(1)) p1(1) = pn(1);
      if (pn(2) > p1(2)) p1(2) = pn(2);
    }


    int tg0[3];
    int tg1[3];
    _tensor_grid.find_element(p0, tg0);
    _tensor_grid.find_element(p1, tg1);

    // all tensor grid elements in the cube [tg0, tg1] may be touched
    // by elem

    for (unsigned int k = tg0[0]; k <= tg1[0]; ++k)
      for (unsigned int l = tg0[1]; l <= tg1[1]; ++l)
        for (unsigned int m = tg0[2]; m <= tg1[2]; ++m)
        {
          int tgrid_el = _tensor_grid.index_to_element(k, l, m);


          vector<const libMesh::Elem*>& ellist = _elem_list[tgrid_el];
          vector<const libMesh::Elem*>::iterator it(find(ellist.begin(), ellist.end(), elem));

          if (it == ellist.end())
            ellist.push_back(elem);
        }
  }

  long unsigned int mem = 0;
  for (unsigned int i = 0; i < _elem_list.size(); ++i)
    mem += _elem_list[i].size();

  mem *= sizeof(const libMesh::Elem*);

  ostringstream os;
  os << timer.elapsed_string() << ", memory usage: "
      << mem / (1024*1024) << " MB " << mem;
  Messages::info(os.str());
}



MeshUtils::GridMapper&
MeshUtils::GridMapper::get_mapper(const libMesh::MeshBase* mesh, const set<ID>& regions)
{
  multimap<const libMesh::MeshBase*, GridMapper*>::iterator it(_mappers.lower_bound(mesh));

  while ((it != _mappers.upper_bound(mesh)) && (it->second->_regids != regions))
    ++it;

  if (it == _mappers.upper_bound(mesh))
  {
    it = _mappers.insert(make_pair(mesh, new GridMapper(mesh, regions)));
  }

  return(*(it->second));
}


MeshUtils::GridMapper&
MeshUtils::GridMapper::get_mapper(const libMesh::MeshBase& mesh, const set<ID>& regions)
{
  return(get_mapper(&mesh, regions));
}

const libMesh::Elem*
MeshUtils::GridMapper::get_element(const libMesh::Point& point) const
{
  const libMesh::Elem* el = NULL;

  int tgrid_el = _tensor_grid.find_element(point);

  if (tgrid_el >= 0)
  {
    // we can assume that _elem_list is assembled when getting here
    const vector<const libMesh::Elem*>& list = _elem_list[tgrid_el];

    for (int i = 0; i < list.size(); ++i)
    {
      const libMesh::Elem* elem = list[i];
      // note: by construction it is inside the bounding box of elem
      //if (MeshUtils::may_belong_to_element(elem, point))
      {
        if (elem->contains_point(point))
        {
          el = elem;
          break;
        }
      }
    }
  }

  return(el);
}



