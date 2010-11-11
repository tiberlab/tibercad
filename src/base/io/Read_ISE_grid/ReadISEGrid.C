// $Id$

#include "ReadISEGrid.h"
#include "MeshRegionInfo.h"
#include "BoundaryRegions.h"
#include "ISE_Vertex.h"
#include "ISE_Edge.h"
#include "ISE_Face.h"
#include "ISE_Element_3D.h"
#include "ISE_Element_2D.h"
#include "ISE_Element_1D.h"
#include "ISE_Element_0D.h"

#include "mesh.h"
#include "elem.h"

#include <cassert>


namespace
{

  // A table to convert ISE element ID to libMesh ID
  vector<libMeshEnums::ElemType> translate;

  void setup_translator(void)
  {
    if (translate.size() == 0)
    {
      translate.resize(11);

      translate[0] = INVALID_ELEM;
      translate[1] = EDGE2;
      translate[2] = TRI3;
      translate[3] = QUAD4;
      translate[4] = INVALID_ELEM;
      translate[5] = TET4;
      translate[6] = PYRAMID5;
      translate[7] = PRISM6;
      translate[8] = HEX8;
      translate[9] = INVALID_ELEM;
      translate[10] = INVALID_ELEM;
    }
  }

}




void ReadISEGrid::integrity_check(float ver, string tp)
{
  if (ver >= 1.0)
  {
    if (tp != "grid")
      throw InitFailedException("Unsupported format for mesh file.");
  }
  else
    throw InitFailedException("Unsupported format for mesh file.");
}


void ReadISEGrid::read(const std::string& name)
{
  std::ifstream in(name.c_str());

  if (!in.good())
    throw InitFailedException(std::string("Cannot read mesh file (") + name + ").");

  this->scan_grid_file(in);

  in.close();
}


void ReadISEGrid::scan_grid_file(std::istream& ISE_INPUT)
{

  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();
  _reg_info.clear();
  _bd_regions.clear();

  string dummy;

  //
  // Read Info block
  //

  do
  {
    ISE_INPUT >> dummy;
  }
  while (dummy != "Info");


  string qualifier, equal, type;
  float version;
  unsigned int dimension,
               nb_vertices,
               nb_edges,
               nb_faces,
               nb_elements,
               nb_regions;

  ISE_INPUT >> dummy >> qualifier >> equal >> version;
  ISE_INPUT >> qualifier >> equal >> type;

  integrity_check(version, type);

  ISE_INPUT >> qualifier >> equal >> dimension;
  ISE_INPUT >> qualifier >> equal >> nb_vertices;
  ISE_INPUT >> qualifier >> equal >> nb_edges;
  ISE_INPUT >> qualifier >> equal >> nb_faces;
  ISE_INPUT >> qualifier >> equal >> nb_elements;
  ISE_INPUT >> qualifier >> equal >> nb_regions;



  // the connection between region names and IDs
  std::map<std::string, unsigned int> regions;

  // the materials
  vector<string> materials(nb_regions);

  ISE_INPUT >> qualifier >> equal >> dummy;

  for (unsigned int i = 0; i < nb_regions; i++)
  {
    string rname;
    ISE_INPUT >> rname;
    // region names are quoted
    rname.erase(0, 1);
    rname.erase(rname.size() - 1, 1);

    // regions are counted from 1
    regions[rname] = i + 1;
  }

  ISE_INPUT >> dummy >> qualifier >> equal >> dummy;

  for (unsigned int i = 0; i < nb_regions; i++)
  {
    ISE_INPUT >> materials[i];
  }

  //
  // End of Info block
  //

  setup_translator();
  mesh.set_mesh_dimension(dimension);


  std::vector<ISE_Vertex*> vertices;
  std::vector<ISE_Edge*> edges;
  std::vector<ISE_Face*> faces;
  std::vector<ISE_Element*> elements_list;

  // *****************************
  // VERTICES  section
  // *****************************
  {

    do
    {
      ISE_INPUT >> dummy;
    }
    while (dummy != "Vertices");

    do
    {
      ISE_INPUT >> dummy;
    }
    while (dummy != "{");


    vector<double> node_coord(3, 0.0);

    vertices.reserve(nb_vertices);

    for (unsigned int i = 0; i < nb_vertices ; i++)
    {

      for (unsigned int j = 0; j < dimension; j++)
      {
        double value;
        ISE_INPUT >> value;
        node_coord[j] =  value;
      }

      vertices.push_back(new ISE_Vertex(node_coord, i));
      mesh.add_point(Point(node_coord[0], node_coord[1], node_coord[2]), i);
    }

  }



  // *****************************
  // EDGES   section
  // *****************************
  if (dimension > 1)
  {

    do
    {
      ISE_INPUT >> dummy;
    }
    while (dummy != "Edges");

    do
    {
      ISE_INPUT >> dummy;
    }
    while (dummy != "{");

    edges.reserve(nb_edges);

    for (unsigned int i = 0; i < nb_edges; i++)
    {
      unsigned int id_1;
      unsigned int id_2;
      ISE_INPUT >> id_1;
      ISE_INPUT >> id_2;

      edges.push_back(new ISE_Edge(vertices[id_1], vertices[id_2]));
    }

  }




  // *****************************
  // FACES   section
  // *****************************
  if (dimension == 3)
  {

    do
    {
      ISE_INPUT >> dummy;
    }
    while (dummy != "Faces");

    do
    {
      ISE_INPUT >> dummy;
    }
    while (dummy != "{");

    int edge_id;
    vector<ISE_Edge*> face_edgs;
    vector<bool> neg_edges;

    faces.reserve(nb_faces);

    for (unsigned int f = 0; f < nb_faces; f++)
    {
      unsigned int n_edg;
      ISE_INPUT >> n_edg;

      face_edgs.resize(n_edg);
      neg_edges.resize(n_edg);

      for (unsigned int e = 0; e < n_edg; e++)
      {
        ISE_INPUT >> edge_id;

        if (edge_id < 0)
        {
          edge_id = -edge_id-1;
          neg_edges[e] = true;
        }
        else
        {
          neg_edges[e] = false;
        }

        face_edgs[e] = edges[edge_id];
      }

      faces.push_back(new ISE_Face(face_edgs, neg_edges));
    }

  }




  // *****************************
  // Elements section
  // *****************************
  {

    // skip Locations section
    do
    {
      ISE_INPUT >> dummy;
    }
    while (dummy != "{");

    do
    {
      ISE_INPUT >> dummy;
    }
    while (dummy != "{");


    elements_list.reserve(nb_elements);

    vector<ISE_Edge*> edge_list;  // local vector
    vector<ISE_Face*> face_list;  // local vector

    vector<bool> negative_edges; // local vector

    for (unsigned int i = 0; i < nb_elements; i++)
    {
      // element type
      unsigned int  elem_type;
      ISE_INPUT >> elem_type ;


      switch (elem_type)
      {
        case 0: // point
        {

          int vertex;
          ISE_INPUT >> vertex;

          ISE_Element* el = new ISE_Element_0D(vertices[vertex]);
          el->set_type(elem_type);
          el->set_dimension(0);

          elements_list.push_back(el);

          break;
        }

        case 1: // segment
        {

          int vertex_0, vertex_1;
          ISE_INPUT >> vertex_0;
          ISE_INPUT >> vertex_1;

          ISE_Element* el = new ISE_Element_1D(vertices[vertex_0], vertices[vertex_1]);
          el->set_type(elem_type);
          el->set_dimension(1);

          elements_list.push_back(el);

          break;
        }

        case 2: // Triangle
        case 3: // Rectangle
        {

          // # of edges of the element = elem_type + 1
          unsigned int num_edges = elem_type + 1;

          negative_edges.resize(num_edges);
          edge_list.resize(num_edges);

          for (unsigned int j = 0; j< num_edges; j++)
          {
            int id;
            ISE_INPUT >> id;

            if (id < 0)
            {
              id = (-id-1);  //  ISE code for inverted edge
              negative_edges[j] = true;
            }
            else
            {
              negative_edges[j] = false;
            }

            edge_list[j] = edges[id];
          } //   edges cycle

          ISE_Element* el = new ISE_Element_2D(edge_list, negative_edges);

          el->set_type(elem_type);
          el->set_dimension(2);

          elements_list.push_back(el);

          break;
        }

        case 4: // 2D polygon
        {

          unsigned int elem_edges;
          ISE_INPUT >> elem_edges;

          negative_edges.resize(elem_edges);
          edge_list.resize(elem_edges);


          for (unsigned int j = 0; j< elem_edges; j++)
          {
            int id;
            ISE_INPUT >> id; // edge id

            if (id < 0)
             {
               id = (-id-1);  //  ISE code for inverted edge
               negative_edges[j] = true;
             }
             else
             {
               negative_edges[j] = false;
             }

             edge_list[j] = edges[id];

          } //   edges cycle

          ISE_Element* el = new ISE_Element_2D(edge_list, negative_edges);

          el->set_type(elem_type);
          el->set_dimension(2);

          elements_list.push_back(el);

          break;
        }

        case 5: //  Tetrahedron
        {
          // this can only be a bulk element

          // tetrahedron has 4 faces
          negative_edges.resize(4);
          face_list.resize(4);

          for (unsigned int j = 0; j < 4; j++)
          {
            int id;
            ISE_INPUT >> id; // face id

            if (id < 0)  //  negative  face  id
            {
              id = (-id-1);   //  ISE code  for inverted face
              //  INVERT THE  ORDER (AND ORIENTATION) OF THE  EDGES OF THIS (NEGATIVE) FACE:
              //  (1,2,3,4) ->   (-4,-3,-2,-1)
              negative_edges[j] = true;
            }
            else
            {
              negative_edges[j] = false;
            }

            face_list[j] = faces[id];

          } // faces cycle


          ISE_Element* el = new  ISE_Element_3D(face_list, negative_edges, elem_type);

          el->set_dimension(3);

          elements_list.push_back(el);

          break;
        }

        case 6:
        case 7: //  Pyramid/Prism
        {

          // they have 5 faces
          negative_edges.resize(5);
          face_list.resize(5);

          for (unsigned int j = 0; j < 5; j++)
          {
            int id;
            ISE_INPUT >> id;  // face id

            if (id < 0)  //  negative  face  id
            {
              id = (-id-1);   //  ISE code  for inverted face
              //  INVERT THE  ORDER (AND ORIENTATION) OF THE EDGES OF THIS (NEGATIVE) FACE :
              //  (1,2,3,4,5) ->   (-5,-4,-3,-2,-1)
              negative_edges[j] = true;
            }
            else
            {
              negative_edges[j] = false;
            }

            face_list[j] = faces[id];

          } //   faces cycle

          ISE_Element* el = new ISE_Element_3D(face_list, negative_edges, elem_type);

          el->set_dimension(3);
          elements_list.push_back(el);

          break;
        }

        case 8: //  Brick
        {

          // has 6 faces
          negative_edges.resize(6);
          face_list.resize(6);

          for (unsigned int j = 0; j < 6; j++)
          {
            int id;
            ISE_INPUT >> id;   // face id

            if (id < 0)  //  negative  face  id
            {
              id = (-id-1);   //  ISE code  for inverted face
              // INVERT THE  ORDER (AND ORIENTATION) OF THE EDGES OF THIS (NEGATIVE) FACE :
              // (1,2,3,4,5,6) ->   (-6,-5,-4,-3,-2,-1)
              negative_edges[j] = true;
            }
            else
            {
              negative_edges[j] = false;
            }

            face_list[j] = faces[id];

          } //   faces cycle

          ISE_Element* el = new ISE_Element_3D(face_list, negative_edges, elem_type);

          el->set_dimension(3);
          elements_list.push_back(el);

          break;
        }

        case 9: // Tetrabrick
        {

          // has 7 faces
          negative_edges.resize(7);
          face_list.resize(7);

          for (unsigned int j = 0; j < 7; j++)
          {
            int id;
            ISE_INPUT >> id; // face id

            if (id < 0)  //  negative  face  id
            {
              id = (-id-1);   //  ISE code  for inverted face
              // INVERT THE ORDER (AND ORIENTATION) OF THE EDGES OF THIS (NEGATIVE) FACE :
              // (1,2,3,4,5,6,7) ->   (-7,-6,-5,-4,-3,-2,-1)
              negative_edges[j] = true;
            }
            else
            {
              negative_edges[j] = false;
            }

            face_list[j] = faces[id];

          } //   faces cycle


          ISE_Element* el = new ISE_Element_3D(face_list, negative_edges, elem_type);

          el->set_dimension(3);
          elements_list.push_back(el);


          break;
        }

        case 10: //  Polyhedron
        {

          unsigned int elem_faces;
          ISE_INPUT >> elem_faces;
          negative_edges.resize(elem_faces);
          face_list.resize(elem_faces);

          for (unsigned int j = 0; j < elem_faces; j++)
          {
            int id;
            ISE_INPUT >> id; // face id

            if (id < 0)  // negative  face  id
            {
              id = (-id-1);   //  ISE code  for inverted face
              // INVERT THE ORDER (AND ORIENTATION) OF THE EDGES OF THIS (NEGATIVE) FACE :
              // (1,2,3,...,n) ->   (-n,...,-3,-2,-1)
              negative_edges[j] = true;
            }
            else
            {
              negative_edges[j] = false;
            }

            face_list[j] = faces[id];

          } // faces cycle

          ISE_Element* el = new ISE_Element_3D(face_list, negative_edges, elem_type);

          el->set_dimension(3);
          elements_list.push_back(el);

          break;
        }

        default:
          throw InitFailedException("Read unexistent ISE element type. Check your mesh file.");
          break;
      }	//  end  switch
    } //  next  element

    ISE_INPUT >> dummy; //  read closing  bracket } of  element section
  }


  // *****************************************************************************
  //                           PHYSICAL REGIONS
  // *****************************************************************************

  vector<unsigned int> boundary_elem;
  vector<unsigned int> edge_elem;

  for (unsigned int i = 0; i < nb_regions; i++)
  {

    // read and discard "Region" keyword
    ISE_INPUT >> dummy;

    // region name in the  form: ("region_1")
    string region_name;
    ISE_INPUT >> region_name;

    // TODO could be done in amore robust way
    string::size_type loc = region_name.find_first_of( '(', 0 );
    if(loc == string::npos)
      throw InitFailedException("in reading ISE grid file "
          " (could not read physical region name).");

    region_name.erase(loc, loc+2);

    loc = region_name.find( '"', 0 );
    if(loc == string::npos)
      throw InitFailedException("in reading ISE grid file "
          " (could not read physical region name).");

    region_name.erase(loc);

    // we assume that we find only names present in the list read in the Info block
    unsigned int phys_reg_id = regions[region_name];

    do
    {
      //  discard "material =..." and "Elements" keyword
      ISE_INPUT >> dummy;
    }
    while (dummy != "(");

    // find number of elements in region
    int numb_region_elements;
    ISE_INPUT >> numb_region_elements;

    do
    {
      ISE_INPUT >> dummy;
    }
    while (dummy != "{");

    for  (unsigned int j = 0; j < numb_region_elements; j++)
    {
      int id;
      ISE_INPUT >> id;

      // set the physical region ID
      elements_list[id]->set_physical_region(phys_reg_id);

      switch (dimension - elements_list[id]->get_dimension())
      {
        case 1: // boundary element
          boundary_elem.push_back(id);
          break;

        case 2: // edge element
          edge_elem.push_back(id);
          break;

        case 3: // node element (only in 3D)
        {
          assert(elements_list[id]->n_nodes() == 1);
          unsigned int nod = (elements_list[id]->get_nodes_id())[0];
          _bd_regions.add_node(mesh.node_ptr(nod), phys_reg_id);
          break;
        }

        default: // bulk element
          // add the subdomain id to the MeshRegionInfo
          _reg_info.add_id(phys_reg_id);
          break;
      }

    } // end of  present region elements

    ISE_INPUT >> dummy; //  read first closing  bracket } of  physical region section
    ISE_INPUT >> dummy; //  read second closing  bracket } of  physical region section

  }  //  next  phys. region


  // add the elements to the mesh
  mesh.reserve_elem(nb_elements);

  unsigned int elem_id_counter = 0;
  for (unsigned int i = 0; i < nb_elements; i++)
  {
    ISE_Element* el = elements_list[i];

    if (el->get_dimension() == dimension)
    {
      libMeshEnums::ElemType eltype = translate[el->get_type()];
      Elem* elem = Elem::build(eltype).release();

      elem->set_id(elem_id_counter);
      mesh.add_elem(elem);

      unsigned int nnodes = el->n_nodes();
      const vector<unsigned int>& elnodes = el->get_nodes_id();
      for (unsigned int n = 0; n < nnodes; n++)
      {
        elem->set_node(n) = mesh.node_ptr(elnodes[n]);
      }

      elem->subdomain_id() =
        static_cast<subdomain_id_type>(el->get_physical_region());

      // this is different from iel: lower dimensional elems aren't added
      elem_id_counter++;

    }

  } // nb_elements


  /**
   * If any lower dimensional elements have been found in the file,
   * try to add them to the mesh.boundary_info as sides and nodes with
   * the respecitve id's (called "physical" in Gmsh).
   */
  if (boundary_elem.size() > 0)
  {
    // create a index of the boundary nodes to easily locate which
    // element might have that boundary
    HashMap<unsigned int, std::vector<unsigned int> >::Type node_index;
    for (unsigned int i = 0; i < boundary_elem.size(); i++)
    {
      const ISE_Element* el = elements_list[boundary_elem[i]];
      for (unsigned int n = 0; n < el->n_nodes(); n++)
        node_index[el->get_nodes_id()[n]].push_back(i);
    }


    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    // iterate over all elements and see which boundary element has
    // the same set of nodes as one of the boundary elements previously read
    for ( ; it != end; ++it)
    {
      const Elem* elem = *it;

      for (unsigned int s = 0; s < elem->n_sides(); s++)
      {
        AutoPtr<Elem> side (elem->build_side(s));
        std::set<unsigned int> side_nodes;
        std::set<unsigned int>::iterator iter = side_nodes.begin();

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
            const ISE_Element* ise_el = elements_list[boundary_elem[bidx]];

            std::set<unsigned int> tmp_set;
            iter = tmp_set.begin();
            for (unsigned int ns = 0; ns < ise_el->n_nodes(); ns++)
              tmp_set.insert(iter, ise_el->get_nodes_id()[ns]);

            if (tmp_set == side_nodes)
            {
              unsigned int id = ise_el->get_physical_region();
              _bd_regions.add_side(elem, s, id);
            }
          }
        }
      } // side loop
    } // element loop
  } // if boundary_elem.size() > 0


  // We do the same for the edges
  if (edge_elem.size() > 0)
  {
    // create a index of the edge nodes to easily locate which
    // element might have that boundary
    HashMap<unsigned int, std::vector<unsigned int> >::Type node_index;
    for (unsigned int i = 0; i < edge_elem.size(); i++)
    {
      const ISE_Element* el = elements_list[edge_elem[i]];
      for (unsigned int n = 0; n < el->n_nodes(); n++)
        node_index[el->get_nodes_id()[n]].push_back(i);
    }

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    // iterate over all elements and see which edge element has
    // the same set of nodes as one of the edge elements previously read
    for ( ; it != end; ++it)
    {
      const Elem* elem = *it;

      for (unsigned int s = 0; s < elem->n_edges(); s++)
      {
        AutoPtr<Elem> edge (elem->build_edge(s));
        std::set<unsigned int> edge_nodes;
        std::set<unsigned int>::iterator iter = edge_nodes.begin();

        // make a set with all nodes from this edge
        // this allows for easy comparison
        for (unsigned int ns = 0; ns < edge->n_nodes(); ns++)
          edge_nodes.insert(iter, edge->node(ns));

        // See whether one of the side node occurs in the list
        // of tagged nodes. If we would loop over all edge
        // nodes, we would just get multiple hits, so taking
        // node 0 is enough to do the job
        unsigned int sn = edge->node(0);
        if (node_index.count(sn) > 0)
        {
          // Loop over all tagged ("physical") "edges" which
          // contain the node sn (typically just 1 to
          // three). For each of these the set of nodes is
          // compared to the current element's edge nodes
          for (unsigned int n = 0; n < node_index[sn].size(); n++)
          {
            unsigned int bidx = node_index[sn][n];
            const ISE_Element* ise_el = elements_list[edge_elem[bidx]];

            std::set<unsigned int> tmp_set;
            iter = tmp_set.begin();
            for (unsigned int ns = 0; ns < ise_el->n_nodes(); ns++)
              tmp_set.insert(iter, ise_el->get_nodes_id()[ns]);


            if (tmp_set == edge_nodes)
            {
              unsigned int id = ise_el->get_physical_region();
              _bd_regions.add_edge(elem, s, id);
            }
          }
        }
      } // edge loop
    } // element loop
  } // if edge_elem.size() > 0


  // set the physical region names
  std::map<std::string, unsigned int>::const_iterator it = regions.begin();
  for ( ; it != regions.end(); ++it)
  {
    ID id = it->second;
    _reg_info.set_name(id, it->first);
    _bd_regions.set_name(id, it->first);
  }

  // clean up
  for (size_t i = 0; i < vertices.size(); i++)
    delete vertices[i];

  for (size_t i = 0; i < edges.size(); i++)
    delete edges[i];

  for (size_t i = 0; i < faces.size(); i++)
    delete faces[i];

  for (size_t i = 0; i < elements_list.size(); i++)
    delete elements_list[i];

}

