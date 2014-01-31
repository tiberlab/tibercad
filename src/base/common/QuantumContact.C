/* QuantumContact.C
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

using namespace std;



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
QuantumContact::init(const ID id,
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

  _normal = get_normal(_area);

  //std::cout<<"normal: "<<_normal(0)<<" "<<_normal(1)<<" "<<_normal(2)<<std::endl; // plot normal

  extend_mesh();
}

void
QuantumContact::plot(const std::string& name)
{
}

QuantumContact::~QuantumContact(void)
{
}


Point
QuantumContact::get_normal(double& area)
{
  unsigned int dim = _mesh->mesh_dimension();

  AutoPtr<FEBase> fe(FEBase::build(dim, FEType()));

  QGauss qrule(dim-1, CONSTANT); // Order 0 rule because in this way we take centroid's normal

  fe->attach_quadrature_rule(&qrule);

  const std::vector<Point>& normal = fe->get_normals(); // normal definition like vector of points

  const std::vector<Real>& JxW = fe->get_JxW(); 

  Point normal1;
  area = 0.0;

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
          fe->reinit(elem, ns);
	  area += JxW[0];
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
QuantumContact::extend_mesh(void)
{
  std::map<ID, ID> nodemap;
  std::map<ID, ID>::iterator nodeit;
  std::vector<ID> nodevec;

  BoundaryRegions::side_iterator it =  _bd_regions->sides_begin(_bd_ids);
  const BoundaryRegions::side_iterator end =  _bd_regions->sides_end(_bd_ids);

  unsigned int numelem = _mesh->n_elem(); // number of element in the mesh
  unsigned int numnode = _mesh->n_nodes();// number of node in the mesh

  // 1D CASE -------------------------------------------------------
  //
  //  --------------------------------------------------------------
  if(_mesh->mesh_dimension() == 1)
  {
    for (unsigned int i = 0 ; i < _mesh->mesh_dimension(); ++i)
    {
      // Adding first node
      const ElementSide& elemside = *++it;

      ID side = elemside.side();

      const Elem* elem = elemside.elem();

      Elem* newelem = _mesh->add_elem(new Edge2);

      newelem->subdomain_id() = _id;

      _elemmap[elem] = newelem;

      _elemsidemap[newelem] = &elemside;

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

      for (unsigned int i = 0 ; i < _mesh->mesh_dimension(); ++i)
      {
        const ElementSide& elemside = *++it;

        const Elem* elem = _elemmap[elemside.elem()];

        Elem* newelem = _mesh->add_elem(new Edge2);

        newelem->subdomain_id() = _id;

        _elemsidemap[newelem] = &elemside;

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
          }
        }
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

      _elemmap[elem] = newelem;

      _elemsidemap[newelem] = &elemside;

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

        std::map<const Elem*, Elem*>::iterator itmap(_elemmap.find(elemside.elem()));

        if ( itmap!=_elemmap.end() )

        {

          const Elem* elem = itmap->second;

          Elem* newelem = _mesh->add_elem(new Quad4);

          _elemsidemap[newelem] = &elemside;

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
      const Point centr = elem->centroid();


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
        // determinat must be positive, otherwise swap points
        if (Deter(_mesh->point(nodevec[0])-centr, _mesh->point(nodevec[1])-centr, _mesh->point(nodevec[2])-centr)<0)
          std::swap(nodevec[1], nodevec[2]);
      }

      if (countnodes==4)
      {
        newelem = _mesh->add_elem(new Hex8);

        // reorder nodes to have a closed loop counter clockwise
        // determinat must be positive, otherwise swap points
        if (Deter(_mesh->point(nodevec[0])-centr, _mesh->point(nodevec[1])-centr, _mesh->point(nodevec[2])-centr)<0)
          std::swap(nodevec[1], nodevec[2]);

        if (Deter(_mesh->point(nodevec[0])-centr, _mesh->point(nodevec[1])-centr, _mesh->point(nodevec[3])-centr)<0)
          std::swap(nodevec[1], nodevec[3]);

        if (Deter(_mesh->point(nodevec[0])-centr, _mesh->point(nodevec[2])-centr, _mesh->point(nodevec[3])-centr)<0)
          std::swap(nodevec[2], nodevec[3]);
      }

      newelem->subdomain_id()=_id;

      _elemmap[elem] = newelem;

      _elemsidemap[newelem] = &elemside;

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

        std::map<const Elem*, Elem*>::iterator itmap(_elemmap.find(elemside.elem()));

        if ( itmap!=_elemmap.end() )
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

          _elemsidemap[newelem] = &elemside;

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

}

void
QuantumContact::set_neighbor_map(void)
{
  // Adding neighbors to new elements (1D, 2D, 3D)
  // Elements on quantum regions have 2 neighobours (one in the original mesh and a new one
  // on quantum contact. This is not a problem since negf is solved with outside elements inactive
  // _elemmap maps elements in the original mesh to quantum_contact

  BoundaryRegions::side_iterator it =  _bd_regions->sides_begin(_bd_ids);
  const BoundaryRegions::side_iterator end =  _bd_regions->sides_end(_bd_ids);

  for ( ; it != end; ++it)
  {
    const ElementSide& elemside = *it;

    ID side = elemside.side();

    std::map<const Elem*, Elem*>::iterator itmap(_elemmap.find(elemside.elem()));

    if ( itmap!=_elemmap.end() )
    {
      if( _mesh->mesh_dimension() == 1 || _mesh->mesh_dimension() == 3)
      {
        (itmap->second)->set_neighbor(0, const_cast<Elem*>(itmap->first));
      }
      else if(_mesh->mesh_dimension() == 2 )
      {
        (itmap->second)->set_neighbor(3, const_cast<Elem*>(itmap->first));
      }
    }
  }

  //write_neighbors();

  //neighbor mapping of the original mesh
  /*{
    std::cout<<"Neighbours of original mesh"<<std::endl;
    BoundaryRegions::side_iterator it =  _bd_regions->sides_begin(_bd_ids);
    const BoundaryRegions::side_iterator end =  _bd_regions->sides_end(_bd_ids);

    for ( ; it != end; ++it)
    {
      const ElementSide& elemside = *it;

      const Elem* elem = elemside.elem();

      if (_rg_ids.count(elem->subdomain_id()))
      {
        std::cout<<elem->id()<<" : " ;
        for (ID nn=0; nn<elem->n_sides() ;nn++)
        {
          if (elem->neighbor(nn)==NULL)
          {
            std::cout<< "0" << ", ";
          }
          else
          {
            std::cout<< elem->neighbor(nn)->id() << ", ";
          }

        }
        std::cout<<std::endl;
      }
    }
  }*/

}

void
QuantumContact::write_neighbors(void) const
{
  MeshBase::element_iterator it = _mesh->level_elements_begin(0);
  const MeshBase::element_iterator end = _mesh->level_elements_end(0);

  for ( ; it != end; ++it)
  {
    Elem* elem = *it;
    ID elid = elem->subdomain_id();

    if(elid==_id)
    {
      std::cout<<elem->id()<<" : " ;
      for (ID nn=0; nn<elem->n_sides() ;nn++)
      {
        if (elem->neighbor(nn)==NULL)
        {
          std::cout<< "0" << ", ";
        }
        else
        {
          std::cout<< elem->neighbor(nn)->id() << ", ";
        }

      }
      std::cout<<std::endl;


    }
  }

}

std::pair<const Elem*, Point >
QuantumContact::project_on_boundary(const Elem* elem, const Point& point ) 
{
  const Elem* sidelem = NULL;
  Point out;

  map<const Elem*, const ElementSide*>::iterator elit(_elemsidemap.find(elem));

  if (elit != _elemsidemap.end())
  {
      std::vector<Point> p;
          
      const ElementSide* elemside = elit->second;
      
      sidelem = elemside->elem();
      
      ID side = elemside->side();

      p.reserve(4);
      p.clear();
      
      for(ID nde = 0; nde < sidelem->n_nodes(); nde++)
      {
        if (sidelem->is_node_on_side(nde, side))
        {
          p.push_back(sidelem->point(nde));
        }
      }
      

      // 1D case=====================================
      if(_mesh->mesh_dimension() == 1)
      {
        out = p[0];
      }
      
      // 2D case=====================================
      //             point
      //            / |
      //         b /  |
      //          /   |
      //         /    |
      //      p[0]---out----p[1]
      //
      //  3D case=============================================
      //
      //
      //
      else if(_mesh->mesh_dimension() == 2)
      {
        
        Point a = p[1]-p[0];
        Point b = point - p[0];     
        double c = a*b;
        double d = a*a;
        Point w = a*c;
        Point y = w/d;
        out = p[0] + y;
      }

      // 3D case=====================================
      else if(_mesh->mesh_dimension() == 3)
      {   
        Point a = p[1]-p[0];
        Point b = point - p[0];
        Point e = p[2]-p[0];
        Point f = a.cross(e);
        double g = b*f;
        double h = f*f;
        Point i = f*g;
        Point l = i/h;
        out = b-l+p[0];
      }
  }
  

  return(make_pair(sidelem, out));
}


//Need to project a point on boundary element
std::pair<const Elem*, std::vector<Point> >
QuantumContact::project_on_boundary(const Elem* elem, const std::vector<Point>& point ) 
{
  std::vector<Point> out(point.size());
  const Elem* sideel = NULL;

  for (unsigned int qp=0; qp<point.size(); qp++)
  {
    std::pair<const Elem*, Point> pp( project_on_boundary(elem, point[qp]) );
 
    out[qp] = pp.second;
    sideel = pp.first;
  }

  return(make_pair(sideel, out));

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
