// $Id$


#include "MeshUtils.h"
#include "Utils.h"
#include "Messages.h"
#include "HashMap.h"
#include "HashSet.h"
#include "RuntimeException.h"
#include "TiberCad.h"

#include "libmesh/mesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_base.h"
#include "libmesh/face_tri3.h"


#include <cassert>
#include <list>
#include <algorithm>
#include <numeric>

#include "libMeshDefs.h"

using namespace std;
using namespace libMesh;

USELIBMESHTYPE(DofObject);


namespace
{
  // trinagle in terms of node indices
  struct triangle_t
  {
    size_t v1, v2, v3;
  };

  // edge in terms of node indices
  struct edge_t
  {
    size_t v1, v2;

    bool operator==(const edge_t& e) const
    {
      return(((v1 == e.v1) && (v2 == e.v2)) ||
             ((v1 == e.v2) && (v2 == e.v1)));
    }

    bool operator!=(const edge_t& e) const
    {
      return(!(*this == e));
    }
  };

  double cross(const Point& O, const Point& A, const Point& B)
  {
    double cr = ((A(0) - O(0)) * (B(1) - O(1)) - (A(1) - O(1)) * (B(0) - O(0)));
    return(cr);
  }

  class Circle
  {
    public:
      Circle(const Point& p1, const Point& p2, const Point& p3)
        : _center(p1), _radius(0)
      {
        // we assume that z coordinate of p1-3 is zero (or the same at least)

        // shift p1 to origin
        Point p2s = p2 - p1;
        Point p3s = p3 - p1;

        double d = 2 * (p2s(0)*p3s(1) - p2s(1)*p3s(0));

        double p2sq = p2s.norm_sq();
        double p3sq = p3s.norm_sq();
        double cx = (p3s(1) * p2sq - p2s(1) * p3sq) / d;
        double cy = (p2s(0) * p3sq - p3s(0) * p2sq) / d;

        _center(0) += cx;
        _center(1) += cy;

        _radius = sqrt(cx*cx + cy*cy);
      }

      const Point& center(void) const
      {
        return(_center);
      }

      double radius(void) const
      {
        return(_radius);
      }

      // returns true if p is inside or on boundary
      bool is_inside(const Point& p)
      {
        Point d = p - _center;

        return((d.norm() - _radius) < 1e-12);
      }


    private:

      Point _center;

      double _radius;
  };
}


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



  unique_ptr<const libMesh::Elem> side_el(elem->side_ptr(side));
  //const libMesh::Elem* side_el = dynamic_cast<libMesh::Elem*>(sobj.get());
  const libMesh::Point& centroid = elem->vertex_average();

  switch (elem->dim())
  {
    case 0:
      break;

    case 1:
      // side is a node
      normal = *side_el->node_ptr(0) - centroid;
      break;

    case 2:
    {
      // side should always be an Edge2
      // normal direction is: p0 + t*(p1 - p0) - centroid, where
      // t gives the intersection between the side and the perpendicular
      // through the element centroid
      libMesh::Point p10((*side_el->node_ptr(1) - *side_el->node_ptr(0)).unit());
      libMesh::Point p03(*side_el->node_ptr(0) - centroid);
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



unique_ptr<libMesh::MeshBase>
MeshUtils::create_boundary_mesh(const libMesh::MeshBase& mesh)
{
  libMesh::MeshBase* bdmesh = new libMesh::ReplicatedMesh(TiberCad::get_mpi_comm(), mesh.mesh_dimension());

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

    for (unsigned int i = 0; i < elem->n_sides(); ++i)
    {
      const libMesh::Elem* nb = elem->neighbor_ptr(i);
      pair<const libMesh::Elem*, const libMesh::Elem*> el_pair(elem, nb);
      if (nb > elem) { el_pair.first = nb; el_pair.second = elem; }

      if ((nb == NULL) || ((nb->subdomain_id() != subdomain) && !sides_added.count(el_pair)))
      {
        sides_added.insert(el_pair);

        // The const_cast is needed here, because elem is const 
        libMesh::Elem* side_el = const_cast<libMesh::Elem*>(elem->build_side_ptr(i, false).release());

        HashMap<unsigned int, unsigned int>::Type::iterator mit;

        for (unsigned int n = 0; n < side_el->n_nodes(); ++n)
        {
          mit = node_id_map.find(side_el->node_id(n));

          unsigned int id = side_el->node_ptr(n)->id();
          if (mit == node_id_map.end())
          {
            node_id_map[side_el->node_id(n)] = node_ctr;
            const libMesh::Node* p = side_el->node_ptr(n);
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

  return unique_ptr<libMesh::MeshBase>(bdmesh);
}



const libMesh::Elem*
MeshUtils::search_element(const libMesh::MeshBase* mesh, const libMesh::Point& point)
{
  return GridMapper::get_mapper(mesh).get_element(point);
}



void
MeshUtils::triangulate_point_set(libMesh::MeshBase& mesh)
{
  size_t n = mesh.n_nodes();

  if (n < 3)
    return;

  if (mesh.mesh_dimension() != 2)
    throw InitFailedException("Can create triangulation only for 2D point sets.");

  // first sort points (we do this on an index vector)
  vector<size_t> node_ids(n);
  iota(node_ids.begin(), node_ids.end(), 0);

  auto compare = [&](size_t i, size_t j)
  {
    const Point& pi = mesh.point(i);
    const Point& pj = mesh.point(j);

    return ((pi(0) < pj(0)) ||
            ((pi(0) == pj(0)) && (pi(1) < pj(1))));

  };

  sort(node_ids.begin(), node_ids.end(), compare);


  // now construct the supertriangle
  auto bb = libMesh::MeshTools::create_nodal_bounding_box(mesh);
  Point diag = bb.max() - bb.min();
  double dmax = (diag(0) > diag(1)) ? diag(0) : diag(1);
  Point mid = 0.5 * (bb.max() + bb.min());

  Point t1(mid);
  t1(0) -= 20*dmax;
  t1(1) -= dmax;

  Point t2(mid);
  t2(1) += 20*dmax;

  Point t3(mid);
  t3(0) += 20*dmax;
  t3(1) -= dmax;


  mesh.add_point(t1);
  mesh.add_point(t3);
  mesh.add_point(t2);

  // now the last 3 points make the supertriangle

  // the list of elements
  vector<triangle_t> elems;
  // a pessimistic guess, since we do not calculate the convex hull
  elems.reserve(2*mesh.n_nodes());

  vector<bool> complete;
  complete.reserve(elems.size());


  elems.push_back({mesh.n_nodes() - 3,
                   mesh.n_nodes() - 2,
                   mesh.n_nodes() - 1});

  complete.push_back(false);


  // we add one point at a time, looping on the ordered id list
  for (size_t i = 0; i < node_ids.size(); ++i)
  {
    dof_id_type node_id = node_ids[i];
    const Point& pi = mesh.point(node_id);

    // the edges that will form the hull around the new point
    vector<edge_t> edges;
    edges.reserve(100);

    size_t ctr = elems.size();

    for (size_t j = 0; j < ctr; ++j)
    {
      if (complete[j])
        continue;

      const Point& p1 = mesh.point(elems[j].v1);
      const Point& p2 = mesh.point(elems[j].v2);
      const Point& p3 = mesh.point(elems[j].v3);
      Circle c(p1, p2, p3);
      double dist = pi(0) - c.center()(0);
      if (dist > c.radius())
        complete[j] = true;

      if (c.is_inside(pi))
      {
        if (edges.size() + 3 > edges.capacity())
          edges.reserve(edges.capacity() + 100);

        const triangle_t& t = elems[j];
        edges.push_back({t.v1, t.v2});
        edges.push_back({t.v2, t.v3});
        edges.push_back({t.v3, t.v1});

        elems[j] = elems[ctr-1];
        elems.pop_back();
        complete[j] = complete[ctr-1];
        complete.pop_back();
        --ctr;
        --j;
      }
    }

    // look for edges multiple times appearing
    // (they are internal, single edges define the hull)
    vector<edge_t> hull;
    set<size_t> skip;
    for (size_t i = 0; i < edges.size(); ++i)
    {
      if (skip.count(i))
        continue;

      bool once = true;

      for (size_t j = i+1; j < edges.size(); ++j)
      {
        if (edges[i] == edges[j])
        {
          skip.insert(j);
          once = false;
          break;
          // because an edge is at most shared by two triangles
        }
      }

      if (once)
        hull.push_back(edges[i]);
    }

    // now add edges from point to hull points
    for (auto&& e : hull)
    {
      elems.push_back({node_id, e.v1, e.v2});
      complete.push_back(false);
    }

  }

  // now add all elements beside the ones with supertriangle
  // vertices

  for (auto&& el : elems)
  {
    if ((el.v1 < node_ids.size()) &&
        (el.v2 < node_ids.size()) &&
        (el.v3 < node_ids.size()))
    {
      Elem* elem = mesh.add_elem(new libMesh::Tri3);
      elem->set_node(0) = mesh.node_ptr(el.v1);
      elem->set_node(1) = mesh.node_ptr(el.v2);
      elem->set_node(2) = mesh.node_ptr(el.v3);
    }
  }

  //mesh.delete_node(mesh.node_ptr(mesh.n_nodes() - 1));
  //mesh.delete_node(mesh.node_ptr(mesh.n_nodes() - 2));
  //mesh.delete_node(mesh.node_ptr(mesh.n_nodes() - 3));

  mesh.allow_renumbering(false);
  mesh.prepare_for_use();


/*

  // upper and lower part of hull
  vector<size_t> U(mesh.n_nodes());
  vector<size_t> L(mesh.n_nodes());
  vector<size_t> T;

  size_t k = 0;

  for (size_t i = 0; i < mesh.n_nodes(); ++i)
  {
    size_t nid = node_ids[i];
    while (k >= 2 && cross(mesh.point(L[k-2]),
                           mesh.point(L[k-1]),
                           mesh.point(nid)) <= 1e-6)
    {
      T.push_back(L[k-2]);
      T.push_back(L[k-1]);
      T.push_back(nid);
      k--;
    }

    L[k++] = nid;
  }
  L.resize(k-1);

  k = 0;

  // note on the check: i >= 0 does not work because it underflows at the next --i
  for (size_t i = mesh.n_nodes()-1; i != static_cast<size_t>(-1); --i)
  {
    size_t nid = node_ids[i];
    while (k >= 2 && cross(mesh.point(U[k-2]),
                           mesh.point(U[k-1]),
                           mesh.point(nid)) <= 0)
    {
      T.push_back(U[k-2]);
      T.push_back(U[k-1]);
      T.push_back(nid);
      k--;
    }
    U[k++] = nid;
  }
  U.resize(k-1);


  for (size_t i = 0; i < T.size(); i += 3)
  {
    Elem* elem = mesh.add_elem(new libMesh::Tri3);
    elem->set_node(0) = mesh.node_ptr(T[i]);
    elem->set_node(1) = mesh.node_ptr(T[i+1]);
    elem->set_node(2) = mesh.node_ptr(T[i+2]);
  }
  //cerr << "# elem : " << mesh.n_elem() << endl;
  */
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
  libMesh::BoundingBox bb(libMesh::MeshTools::create_bounding_box(*_mesh));

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

  libMesh::MeshBase::const_element_iterator it = _mesh->elements_begin();
  const libMesh::MeshBase::const_element_iterator end = _mesh->elements_end();

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

    for (int k = tg0[0]; k <= tg1[0]; ++k)
      for (int l = tg0[1]; l <= tg1[1]; ++l)
        for (int m = tg0[2]; m <= tg1[2]; ++m)
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
      << mem / (1024*1024) << " MB ";
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
    libMesh::Point p(point);
    switch (_mesh->mesh_dimension())
    {
      case 1:
        p(1) = 0;
        [[fallthrough]];

      case 2:
        p(2) = 0;
        break;
    }

    for (int i = 0; i < list.size(); ++i)
    {
      const libMesh::Elem* elem = list[i];
      // note: by construction it is inside the bounding box of elem
      //if (MeshUtils::may_belong_to_element(elem, point))
      {
        if (elem->contains_point(p))
        {
          el = elem;
          break;
        }
      }
    }
  }



  return(el);
}



