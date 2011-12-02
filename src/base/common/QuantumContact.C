/*
 * QuantumContact.C
 *
 *  Created on: Nov 3, 2011
 *      Author: fpalomba
 */
#include "QuantumContact.h"
#include "Device.h"
#include "BoundaryRegions.h"
#include "DataOutput.h"
#include "Messages.h"
#include "MeshRegionInfo.h"
#include "Material.h"

#include "mesh.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "node.h"
#include "edge_edge2.h"
#include "face_quad4.h"
#include "elem.h"
#include "cell_prism6.h"
#include "cell_hex8.h"
#include "gmsh_io.h"
#include "dof_map.h"


#include <fstream>
#include <set>


QuantumContact::QuantumContact(void)
: _device(NULL),
  _mesh(NULL),
  _bd_regions(NULL),
  _id(0)
{

}


QuantumContact*
QuantumContact::create(void)
{
  QuantumContact* qc = new QuantumContact();
  return qc;
}

void
QuantumContact::init(ID id,
    const std::string& name,
    Device* device,
    BoundaryRegions* bd_regions,
    const std::vector<ID>& rg_ids_v,
    const std::vector<ID>& bd_ids_v,
    double length)
{
  _device = device;
  _mesh = &(_device->get_mesh());
  _bd_regions = bd_regions;
  _id = id;
  _name = name;

  _bd_ids = vec2set(bd_ids_v);

  _rg_ids = vec2set(rg_ids_v);

  _length = length;

  unsigned int num_sides;

  _normal = get_normal(num_sides);

  std::cout<<"normal: "<<_normal(0)<<" "<<_normal(1)<<" "<<_normal(2)<<std::endl; // plot normal

  extend_mesh();

  std::cout<<"active elements "<< _mesh->n_active_elem()<<std::endl;
}

void
QuantumContact::plot(const std::string& name)
{
}

QuantumContact::~QuantumContact(void)
{
}


Point
QuantumContact::get_normal(unsigned int& count)
{
  unsigned int dim = _mesh->mesh_dimension();

  AutoPtr<FEBase> fe(FEBase::build(dim, FEType()));

  QGauss qrule(dim-1, CONSTANT); // Order 0 rule because in this way we take centroid's normal

  fe->attach_quadrature_rule(&qrule);

  const std::vector<Point>& normal = fe->get_normals(); // normal definition like vector of points

  Point normal1;
  count=0;

  MeshBase::const_element_iterator it = _mesh->level_elements_begin(0);
  const MeshBase::const_element_iterator end = _mesh->level_elements_end(0);

  bool found_first=false;

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    ID elid = elem->subdomain_id();

    if(_rg_ids.count(elid))
    {
      for (ID ns = 0; ns<elem->n_sides(); ns++)  // loop over element side
      {
        ID side_id = _bd_regions->get_side_id(elem, ns);
        if ( side_id != INVALID_ID && _bd_ids.count(side_id) ) //normal computation
        {
          count++;
          fe->reinit(elem, ns);
          if (!found_first)
          {
            found_first=true;
            normal1=normal[0];
          }
          else
          {
            if(!Utils::almost_equal::compare(normal1*normal[0],1.0)) //parallelism check
            {
              Messages::error("Non parallel face normals");
            }
          }
        }
      }
    }
  }
  return normal1;

}

void
QuantumContact::activate_elements(void)
{
  MeshBase::element_iterator it = _mesh->level_elements_begin(0);
  const MeshBase::element_iterator end = _mesh->level_elements_end(0);

  for ( ; it != end; ++it)
  {
    Elem* elem = *it;
    ID elid = elem->subdomain_id();

    if(elid==_id)
    {
      elem->set_refinement_flag(Elem::DO_NOTHING);
    }
  }
}

void
QuantumContact::inactivate_elements(void)
{
  MeshBase::element_iterator it = _mesh->level_elements_begin(0);
  const MeshBase::element_iterator end = _mesh->level_elements_end(0);

  for ( ; it != end; ++it)
  {
    Elem* elem = *it;
    ID elid = elem->subdomain_id();

    if(elid==_id)
    {
     elem->set_refinement_flag(Elem::INACTIVE);
    }
  }
}


void
QuantumContact::extend_mesh(void)
{
  std::map<ID, ID> nodemap;
  std::map<const Elem*, Elem*> elemmap;
  std::map<ID, ID>::iterator nodeit;
  std::vector<ID> nodevec;

  BoundaryRegions::side_iterator it =  _bd_regions->sides_begin(_bd_ids);
  const BoundaryRegions::side_iterator end =  _bd_regions->sides_end(_bd_ids);

  unsigned int numelem = _mesh->n_elem(); // number of element in the mesh
  unsigned int numnode = _mesh->n_nodes();// number of node in the mesh
  //std::cout<<"num node: "<<numnode<<std::endl;
  //std::cout<<"num elem: "<<numelem<<std::endl;
  unsigned int cnt = 1;

  // 1D CASE -------------------------------------------------------
  //
  //  --------------------------------------------------------------
  if(_mesh->mesh_dimension() == 1)
  {
    for ( ; it != end; ++it)
    {
      // Adding first node
      const ElementSide& elemside = *it;

      ID side = elemside.side();

      const Elem* elem = elemside.elem();

      Elem* newelem = _mesh->add_elem(new Edge2);

      newelem->subdomain_id()=_id;

      elemmap[elem] = newelem;

      ID id = 0; ID inc = 1;

      for(ID nde = 0; nde < elem->n_nodes(); nde++)
      {

        if (elem->is_node_on_side(nde, side))
        {
          ID nodeid = elem->node(nde);
          nodeit = nodemap.find(nodeid);

          if(nodeit == nodemap.end())
          {
            _mesh->add_point(_mesh->point(nodeid)+(_normal*_length),numnode);
            nodeit = nodemap.insert(std::make_pair(nodeid, numnode)).first;
            numnode++;
          }

          ID num = nodeit->second;
          newelem->set_node(id) = _mesh->node_ptr(nodeid);
          newelem->set_node(id+inc) = _mesh->node_ptr(num);
        }
      }
    }
    // Adding second node
    nodemap.clear();
    {
      BoundaryRegions::side_iterator it =  _bd_regions->sides_begin(_bd_ids);
      const BoundaryRegions::side_iterator end =  _bd_regions->sides_end(_bd_ids);

      for ( ; it != end; ++it)
      {
        const ElementSide& elemside = *it;

        const Elem* elem = elemmap[elemside.elem()];

        Elem* newelem = _mesh->add_elem(new Edge2);

        newelem->subdomain_id()=_id;

        ID id = 0; ID inc = 1;

        ID nodeid = elem->node(0);

        _mesh->add_point(_mesh->point(nodeid)+(_normal*_length),numnode);
        nodeit = nodemap.insert(std::make_pair(nodeid, numnode)).first;
        numnode++;

        ID num = nodeit->second;
        newelem->set_node(id) = _mesh->node_ptr(nodeid);
        newelem->set_node(id+inc) = _mesh->node_ptr(num);
      }
    }
  }

  // 2D CASE -------------------------------------------------------
  //
  //  --------------------------------------------------------------
  if(_mesh->mesh_dimension() == 2)
  {
    for ( ; it != end; ++it)
    {
      const ElementSide& elemside = *it;

      ID side = elemside.side();

      const Elem* elem = elemside.elem();

      if (!_rg_ids.count(elem->subdomain_id())) continue;

      Elem* newelem = _mesh->add_elem(new Quad4);

      newelem->subdomain_id()=_id;

      elemmap[elem] = newelem;

      _elemmap[newelem] = &elemside;

      ID id = 0; ID inc = 1;

      //iterates on nodes of element
      for(ID nde = 0; nde < elem->n_nodes(); nde++)
      {
        //check whether node is on the boundary
        if (elem->is_node_on_side(nde, side))
        {
          ID nodeid = elem->node(nde);

          nodeit = nodemap.find(nodeid);

          if(nodeit == nodemap.end())
          {
            _mesh->add_point(_mesh->point(nodeid)+(_normal*_length),numnode);
            nodeit = nodemap.insert(std::make_pair(nodeid, numnode)).first; //insert new node in nodemap
            numnode++;
          }
          ID num = nodeit->second;
          newelem->set_node(id) = _mesh->node_ptr(nodeid);
          newelem->set_node(id+inc) = _mesh->node_ptr(num);
          id+=3; inc=-1; // for this increment we look the id node disposition in Libmesh class Quad4
        }
      }
    }
    // Adding second row of nodes
    nodemap.clear();
    {
      BoundaryRegions::side_iterator it =  _bd_regions->sides_begin(_bd_ids);
      const BoundaryRegions::side_iterator end =  _bd_regions->sides_end(_bd_ids);

      for ( ; it != end; ++it)
      {
        const ElementSide& elemside = *it;

        const Elem* elem = elemmap[elemside.elem()];

        Elem* newelem = _mesh->add_elem(new Quad4);

        _elemmap[newelem] = &elemside;

        newelem->subdomain_id()=_id;

        ID id = 0; ID inc = 1;

        for(ID nde = 0; nde < elem->n_nodes(); nde++)
        {

          if (elem->is_node_on_side(nde, 1))
          {
            ID nodeid = elem->node(nde);
            nodeit = nodemap.find(nodeid);

            if(nodeit == nodemap.end())
            {
              _mesh->add_point(_mesh->point(nodeid)+(_normal*_length),numnode);
              nodeit = nodemap.insert(std::make_pair(nodeid, numnode)).first;
              numnode++;
            }

            ID num = nodeit->second;
            newelem->set_node(id) = _mesh->node_ptr(nodeid);
            newelem->set_node(id+inc) = _mesh->node_ptr(num);
            id+=3; inc=-1;
          }
        }
      }
    }
  }
  // 3D CASE -------------------------------------------------------------
  //
  // ---------------------------------------------------------------------

  if(_mesh->mesh_dimension() == 3)
  {

    ID countnodes = 0;

    nodevec.reserve(4);
    // iterates on boundary elements
    for ( ; it != end; ++it)
    {
      const ElementSide& elemside = *it;

      ID side = elemside.side();

      const Elem* elem = elemside.elem();

      if (!_rg_ids.count(elem->subdomain_id())) continue;

      nodevec.clear();
      countnodes = 0;

      //iterates on nodes of element
      for(ID nde = 0; nde < elem->n_nodes(); nde++)
      {
        //check whether node is on the boundary
        //new nodes are added to the mesh and stored in nodevec[]
        if (elem->is_node_on_side(nde, side))
        {
          ID nodeid = elem->node(nde);

          nodeit = nodemap.find(nodeid);

          if(nodeit == nodemap.end())
          {
            _mesh->add_point(_mesh->point(nodeid)+(_normal*_length),numnode);
            nodeit = nodemap.insert(std::make_pair(nodeid, numnode)).first;
            numnode++;
          }
          ID num = nodeit->second;
          nodevec[countnodes]=nodeid;
          countnodes++;
        }
      }

      Elem* newelem;

      if (countnodes==3)
      {
        newelem = _mesh->add_elem(new Prism6);

        // reorder nodes to have a closed loop counter clockwise
        // determinat must be negative, otherwise swap points
        if (Deter(_mesh->point(nodevec[0]), _mesh->point(nodevec[1]), _mesh->point(nodevec[2]))<0)
          std::swap(nodevec[1], nodevec[2]);
      }

      if (countnodes==4)
      {
        newelem = _mesh->add_elem(new Hex8);

        // reorder nodes to have a closed loop counter clockwise
        // determinat must be negative, otherwise swap points
        if (Deter(_mesh->point(nodevec[0]), _mesh->point(nodevec[1]), _mesh->point(nodevec[2]))<0)
          std::swap(nodevec[1], nodevec[2]);

        if (Deter(_mesh->point(nodevec[0]), _mesh->point(nodevec[1]), _mesh->point(nodevec[3]))<0)
          std::swap(nodevec[1], nodevec[3]);

        if (Deter(_mesh->point(nodevec[0]), _mesh->point(nodevec[2]), _mesh->point(nodevec[3]))<0)
          std::swap(nodevec[2], nodevec[3]);
      }

      newelem->subdomain_id()=_id;

      elemmap[elem] = newelem;

      _elemmap[newelem] = &elemside;

      for (ID id = 0; id < countnodes; id++)
      {
        newelem->set_node(id) = _mesh->node_ptr(nodevec[id]);
        newelem->set_node(id+countnodes) = _mesh->node_ptr(nodemap[nodevec[id]]);
      }
    }
    // Add second row

    nodemap.clear();
    {
      BoundaryRegions::side_iterator it =  _bd_regions->sides_begin(_bd_ids);
      const BoundaryRegions::side_iterator end =  _bd_regions->sides_end(_bd_ids);

      for ( ; it != end; ++it)
      {

        const ElementSide& elemside = *it;

        std::map<const Elem*, Elem*>::iterator itmap(elemmap.find(elemside.elem()));

        if ( itmap!=elemmap.end() )
        {
          const Elem* elem = itmap->second;

          Elem* newelem;

          if (countnodes == 3)
          {
            newelem = _mesh->add_elem(new Prism6);
          }

          if (countnodes == 4)
          {
            newelem = _mesh->add_elem(new Hex8);
          }

          _elemmap[newelem] = &elemside;

          newelem->subdomain_id()=_id;

          ID id = 0;

          for(ID nde = 0; nde < elem->n_nodes(); nde++)
          {
            // check if node is on 'upper' side of previous element
            // is 4 for Prism6; 5 for Hex8
            if (elem->is_node_on_side(nde, countnodes+1))
            {
              ID nodeid = elem->node(nde);

              nodeit = nodemap.find(nodeid);

              if(nodeit == nodemap.end())
              {
                _mesh->add_point(_mesh->point(nodeid)+(_normal*_length),numnode);
                nodeit = nodemap.insert(std::make_pair(nodeid, numnode)).first;
                numnode++;
              }

              ID num = nodeit->second;
              newelem->set_node(id) = _mesh->node_ptr(nodeid);
              newelem->set_node(id+countnodes) = _mesh->node_ptr(num);
              id++;
            }
          }
        }
      }
    }
  }


  // Prepare for use (1D, 2D, 3D)
  _mesh->prepare_for_use(true);


  // Adding neighbors to new elements (1D, 2D, 3D)
  {
    BoundaryRegions::side_iterator it =  _bd_regions->sides_begin(_bd_ids);
    const BoundaryRegions::side_iterator end =  _bd_regions->sides_end(_bd_ids);

    for ( ; it != end; ++it)
    {
      const ElementSide& elemside = *it;

      ID side = elemside.side();

      std::map<const Elem*, Elem*>::iterator itmap(elemmap.find(elemside.elem()));

      if ( itmap!=elemmap.end() )
      {
        if( _mesh->mesh_dimension() == 1 || _mesh->mesh_dimension() == 3)
        {
          (itmap->second)->set_neighbor(0, const_cast<Elem*>(itmap->first)); //elem);
        }
        else if(_mesh->mesh_dimension() == 2 )
        {
          (itmap->second)->set_neighbor(3, const_cast<Elem*>(itmap->first)); //elem);
        }
      }
    }
  }

  /*
  //neighbor dei nuovi elementi aggiunti
  std::cerr<<"New elements : "<<std::endl;
  {
    for ( ID el=numelem; el< _mesh->n_elem(); ++el)
    {

      //std::cerr<<_mesh->elem(el)->id()<<" : " ;

      for (ID nn=0; nn<_mesh->elem(el)->n_neighbors() ;nn++)
      {
        if (_mesh->elem(el)->neighbor(nn)==NULL)
        {
          std::cerr<< "0" << " ";
        }
        else
        {
          std::cerr<< _mesh->elem(el)->neighbor(nn)->id() << " ";
        }

      }
      std::cerr<<std::endl;

     }
  }
   */
  //neighbor mapping
  {
    BoundaryRegions::side_iterator it =  _bd_regions->sides_begin(_bd_ids);
    const BoundaryRegions::side_iterator end =  _bd_regions->sides_end(_bd_ids);

    for ( ; it != end; ++it)
    {
      const ElementSide& elemside = *it;

      const Elem* elem = elemside.elem();

      if (_rg_ids.count(elem->subdomain_id()))
      {
        //std::cerr<<elem->id()<<" : " ;
        for (ID nn=0; nn<elem->n_neighbors() ;nn++)
        {
          if (elem->neighbor(nn)==NULL)
          {
            //std::cerr<< "0" << " ";
          }
          else
          {
            //std::cerr<< elem->neighbor(nn)->id() << " ";
          }

        }
        //std::cerr<<std::endl;
      }
    }
  }

}


double
QuantumContact::Deter (const Point& P1, const Point& P2, const Point& P3) //computation of a 3X3 matrix's determinant
{

  double d =   P1(0)*P2(1)*P3(2) + P2(0)*P3(1)*P1(2)
                 + P3(0)*P1(1)*P2(2) - P1(2)*P2(1)*P3(0)
                 - P2(2)*P3(1)*P1(0) - P3(2)*P1(1)*P2(0);
  return d;
}

std::set<ID>
QuantumContact:: vec2set(const std::vector<ID>& vec) // Transform a vector in a set
{
  std::set<ID> set;

  for(unsigned int n = 0; n < vec.size(); n++)
  {
    set.insert(vec[n]);
  }

  return set;
}

std::vector<ID>
QuantumContact:: set2vec(const std::set<ID>& set) // Transform a set in a vector
{
  std::vector<ID> vec;

  std::set<ID>::iterator it = set.begin();
  std::set<ID>::iterator end = set.end();

  for( ; it!=end; ++it)
  {
    vec.push_back(*it);
  }

  return vec;
}

std::pair<const Elem*, Point>
QuantumContact::project_on_boundary(const Elem* elem,const Point& point )
{
  std::vector<Point> p;
  Point out;

  p.reserve(4);

  const ElementSide* elemside = _elemmap[elem];

  const Elem* sidelem = elemside->elem();

  ID side = elemside->side();

  p.clear();

  for(ID nde = 0; nde < sidelem->n_nodes(); nde++)
  {
    if (sidelem->is_node_on_side(nde, side))
    {
      p.push_back(sidelem->point(nde));
    }
  }

  Point a = p[1]-p[0];
  Point b = point - p[0];

  if(_mesh->mesh_dimension() == 2)
  {
    double c = a*b;
    double d = a*a;
    Point w = a*c;
    Point y = w/d;
    out = p[0] + y;
  }

  if(_mesh->mesh_dimension() == 3)
  {
    Point e = p[2]-p[0];
    Point f = a.cross(e);
    double g = b*f;
    double h = f*f;
    Point i = f*g;
    Point l = i/h;
    out = b-l+p[0];
  }
 return std::make_pair(sidelem,out);
}







