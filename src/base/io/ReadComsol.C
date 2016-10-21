// $Id: ReadComsol.C 2176 2010-11-29 11:28:27Z maufder $

#include "ReadComsol.h"
#include "MeshRegionInfo.h"
#include "BoundaryRegions.h"
#include "Utils.h"
#include "InitFailedException.h"


// C++ includes
#include <fstream>
#include <sstream>
#include <set>
#include <cstring> // memcpy, strncmp

// Local includes
#include "elem.h"
#include "mesh_base.h"


using namespace std;


// anonymous namespace to hold local data
namespace
{

  /**
   * Defines a structure to hold boundary element information.
   *
   * We use a set because it keeps the nodes unique and ordered, and can be
   * easily compared to another set of nodes (the ones on the element side)
   */
  struct BoundaryElementInfo {
      set<unsigned int> nodes;
      unsigned int id;
  };

  struct Point2D {
    Point2D(double xcoord = 0, double ycoord = 0) :
      x(xcoord), y(ycoord) {}

    int index;
    double x;
    double y;
  };

  struct Line {
    int index;
    int point1;
    int point2;
    int domain;
  };

  struct Triangle {

    Triangle(int a, int b, int c) :
      point1(a), point2(b), point3(c) {};

    int index;
    int point1;
    int point2;
    int point3;
    int domain;
  };


} // end anonymous namespace


// ReadComsol  members
void ReadComsol::read(const string& name)
{
  ifstream in(name.c_str());

  if (!in.good())
    throw InitFailedException(string("Cannot read mesh file (") + name + ").");

  this->read_mesh(in);
}


void ReadComsol::read_mesh(istream& in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert(libMesh::global_processor_id() == 0);

  libmesh_assert(in.good());


  // clear any data in the mesh
  libMesh::MeshBase& mesh = MeshInput<libMesh::MeshBase>::mesh();
  mesh.clear();
  _reg_info.clear();
  _bd_regions.clear();


  // we will try to get the mesh dimension from the file
  //unsigned int dim = mesh.mesh_dimension();

  // some variables
  string  buf;
  int nnodes, nelem, nsides;

  // map to hold the physical names and dimensions (if found)
  map<unsigned int, string> phys_names;



  getline(in, buf);
  while ((buf.find("# sdim") == string::npos) && !in.eof())
    getline(in, buf);

  int dim = Utils::convert<int>(buf.substr(0, buf.find(" ")));
  mesh.set_mesh_dimension(dim);

  //
  // read mesh nodes
  //

  getline(in, buf);
  while ((buf.find("number of mesh points") == string::npos) && !in.eof())
    getline(in, buf);

  nnodes = Utils::convert<int>(buf.substr(0, buf.find(" ")));
  mesh.reserve_nodes(nnodes);

  while ((buf.compare(0, 24, "# Mesh point coordinates") != 0) && !in.eof())
    getline(in, buf);

  for (int n = 0; n < nnodes; ++n)
  {
    double x, y, z = 0;
    in >> x >> y;
    mesh.add_point(Point(x, y, z), n);
  }

  //
  // read lines
  //

  while ((buf.compare(0, 17, "3 edg # type name") != 0) && !in.eof())
    getline(in, buf);

  getline(in, buf);
  while ((buf.find("# number of elements") == string::npos) && !in.eof())
    getline(in, buf);

  nsides = Utils::convert<int>(buf.substr(0, buf.find(" ")));

  vector<BoundaryElementInfo> bd_elems(nsides);

  while ((buf.compare(0, 10, "# Elements") != 0) && !in.eof())
    getline(in, buf);

  for (int n = 0; n < nsides; ++n)
  {
    int a, b;

    in >> a >> b;
    BoundaryElementInfo& binfo = bd_elems[n];
    binfo.nodes.insert(a);
    binfo.nodes.insert(b);
  }

  getline(in, buf);
  while ((buf.find("# number of geometric entity indices") == string::npos) && !in.eof())
    getline(in, buf);

  int nindices = Utils::convert<int>(buf.substr(0, buf.find(" ")));
  if (nindices != nsides)
    throw InitFailedException("Comsol mesh has inconsistent number of boundary domain indices");

  getline(in, buf);
  while ((buf.compare(0, 26, "# Geometric entity indices") != 0) && !in.eof())
    getline(in, buf);

  for (int n = 0; n < nsides; ++n)
  {
    int id;
    in >> id;
    BoundaryElementInfo& binfo = bd_elems[n];
    binfo.id = id;
  }

  //
  // read elements
  //

  getline(in, buf);
  while ((buf.compare(0, 17, "3 tri # type name") != 0) && !in.eof())
    getline(in, buf);

  getline(in, buf);
  while ((buf.find("# number of elements") == string::npos) && !in.eof())
    getline(in, buf);

  nelem = Utils::convert<int>(buf.substr(0, buf.find(" ")));
  mesh.reserve_elem(nelem);

  getline(in, buf);
  while ((buf.compare(0, 10, "# Elements") != 0) && !in.eof())
    getline(in, buf);

  for (int n = 0; n < nelem; ++n)
  {
    int a, b, c;
    in >> a >> b >> c;
    Elem* elem = Elem::build(libMeshEnums::TRI3).release();
    elem->set_id(n);
    mesh.add_elem(elem);
    elem->set_node(0) = mesh.node_ptr(a);
    elem->set_node(1) = mesh.node_ptr(b);
    elem->set_node(2) = mesh.node_ptr(c);
  }

  // we put the domain IDs into a set, so we can compare
  // afterwards boundary IDs with bulk domain IDs
  set<int> domain_ids;

  getline(in, buf);
  while ((buf.find("# number of geometric entity indices") == string::npos) && !in.eof())
    getline(in, buf);

  nindices = Utils::convert<int>(buf.substr(0, buf.find(" ")));
  if (nindices != nelem)
    throw InitFailedException("Comsol mesh has inconsistent number of domain indices");

  getline(in, buf);
  while ((buf.compare(0, 26, "# Geometric entity indices") != 0) && !in.eof())
    getline(in, buf);

  for (int n = 0; n < nelem; ++n)
  {
    int id;
    in >> id;
    _reg_info.add_id(id);
    domain_ids.insert(id);
    mesh.elem(n)->subdomain_id() = id;
  }


  if (bd_elems.size() > 0)
  {
    // create a index of the boundary nodes to easily locate which
    // element might have that boundary
    HashMap<unsigned int, vector<unsigned int> >::Type node_index;
    for (unsigned int i = 0; i < bd_elems.size(); i++)
    {
      BoundaryElementInfo& binfo = bd_elems[i];
      set<unsigned int>::iterator iter = binfo.nodes.begin();
      for ( ; iter != binfo.nodes.end(); ++iter)
        node_index[*iter].push_back(i);
    }


    libMesh::MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const libMesh::MeshBase::const_element_iterator end = mesh.active_elements_end();

    // iterate over all elements and see which boundary element has
    // the same set of nodes as one of the boundary elements previously read
    for ( ; it != end; ++it)
    {
      const Elem* elem = *it;

      for (unsigned int s = 0; s < elem->n_sides(); s++)
        //if (elem->neighbor(s) == NULL)
      {
        libMesh::UniquePtr<Elem> side (elem->build_side(s));
        set<unsigned int> side_nodes;
        set<unsigned int>::iterator iter = side_nodes.begin();

        // make a set with all nodes from this side
        // this allows for easy comparison
        for (unsigned int ns = 0; ns < side->n_nodes(); ns++)
          side_nodes.insert(iter, side->node(ns));

        // See whether one of the side node occurs in the list
        // of tagged nodes. If we would loop over all side
        // nodes, we would just get multiple hits, so taking
        // node 0 is enough to do the job
        unsigned int sn = side->node(0);
        if (node_index.count(sn) > 0)
        {
          // Loop over all tagged ("physical") "sides" which
          // contain the node sn (typically just 1 to
          // three). For each of these the set of nodes is
          // compared to the current element's side nodes
          for (unsigned int n = 0; n < node_index[sn].size(); n++)
          {
            unsigned int bidx = node_index[sn][n];
            if (bd_elems[bidx].nodes == side_nodes)
            {
              unsigned int id = bd_elems[bidx].id;
              _bd_regions.add_side(elem, s, id);
            }
          }
        }

      }
    }
  }



  // set the physical region names
  map<unsigned int, string>::iterator it = phys_names.begin();
  for ( ; it != phys_names.end(); ++it)
  {
    ID id = it->first;

    _reg_info.set_name(id, it->second);
    _bd_regions.set_name(id, it->second);
  }


}




