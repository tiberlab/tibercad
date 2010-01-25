// $Id$

//
// Quite a bit of this code has been taken from libmesh
//

#include "ReadGMSH.h"
#include "MeshRegionInfo.h"
#include "BoundaryRegions.h"
#include "InitFailedException.h"


// C++ includes
#include <fstream>
#include <sstream>
#include <set>
#include <cstring> // std::memcpy, std::strncmp

// Local includes
#include "elem.h"
#include "mesh_base.h"


// anonymous namespace to hold local data
namespace
{

  /**
   * Defines a structure to hold boundary element information.
   *
   * We use a set because it keeps the nodes unique and ordered, and can be
   * easily compared to another set of nodes (the ones on the element side)
   */
  struct boundaryElementInfo {
      std::set<unsigned int> nodes;
      unsigned int id;
  };

  /**
   * Defines mapping from libMesh element types to Gmsh element types.
   */
  struct elementDefinition {
      std::string label;
      std::vector<unsigned int> nodes;
      ElemType type;
      unsigned int exptype;
      unsigned int dim;
      unsigned int nnodes;
  };


  // maps from a libMesh element type to the proper
  // Gmsh elementDefinition.  Placing the data structure
  // here in this anonymous namespace gives us the
  // benefits of a global variable without the nasty
  // side-effects
  std::map<ElemType, elementDefinition> eletypes_exp;
  std::map<unsigned int, elementDefinition> eletypes_imp;



  // ------------------------------------------------------------
  // helper function to initialize the eletypes map
  void init_eletypes ()
  {
    if (eletypes_exp.empty() && eletypes_imp.empty())
    {
      // This should happen only once.  The first time this method
      // is called the eletypes data struture will be empty, and
      // we will fill it.  Any subsequent calls will find an initialized
      // eletypes map and will do nothing.

      //==============================
      // setup the element definitions
      elementDefinition eledef;

      // use "swap trick" from Scott Meyer's "Effective STL" to initialize
      // eledef.nodes vector

      // POINT (only Gmsh)
      {
        eledef.exptype = 15;
        eledef.dim     = 0;
        eledef.nnodes  = 1;
        eledef.nodes.clear();

        // import only
        eletypes_imp[15] = eledef;
      }

      // EDGE2
      {
        eledef.type    = EDGE2;
        eledef.dim     = 1;
        eledef.nnodes  = 2;
        eledef.exptype = 1;
        eledef.nodes.clear();

        eletypes_exp[EDGE2] = eledef;
        eletypes_imp[1]     = eledef;
      }

      // EDGE3
      {
        eledef.type    = EDGE3;
        eledef.dim     = 1;
        eledef.nnodes  = 3;
        eledef.exptype = 8;
        eledef.nodes.clear();

        eletypes_exp[EDGE3] = eledef;
        eletypes_imp[8]     = eledef;
      }

      // TRI3
      {
        eledef.type    = TRI3;
        eledef.dim     = 2;
        eledef.nnodes  = 3;
        eledef.exptype = 2;
        eledef.nodes.clear();

        eletypes_exp[TRI3] = eledef;
        eletypes_imp[2] = eledef;
      }

      // TRI6
      {
        eledef.type    = TRI6;
        eledef.dim     = 2;
        eledef.nnodes  = 6;
        eledef.exptype = 9;
        eledef.nodes.clear();

        eletypes_exp[TRI6] = eledef;
        eletypes_imp[9]    = eledef;
      }

      // QUAD4
      {
        eledef.type    = QUAD4;
        eledef.dim     = 2;
        eledef.nnodes  = 4;
        eledef.exptype = 3;
        eledef.nodes.clear();

        eletypes_exp[QUAD4] = eledef;
        eletypes_imp[3]     = eledef;
      }

      // QUAD8
      {
        eledef.type    = QUAD8;
        eledef.dim     = 2;
        eledef.nnodes  = 8;
        eledef.exptype = 100;
        const unsigned int nodes[] = {1,2,3,4,5,6,7,8};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[QUAD8] = eledef;
        eletypes_imp[10]    = eledef;
      }

      // QUAD9
      {
        eledef.type    = QUAD9;
        eledef.dim     = 2;
        eledef.nnodes  = 9;
        eledef.exptype = 10;
        eledef.nodes.clear();

        eletypes_exp[QUAD9] = eledef;
        eletypes_imp[10]    = eledef;
      }

      // HEX8
      {
        eledef.type    = HEX8;
        eledef.dim     = 3;
        eledef.nnodes  = 8;
        eledef.exptype = 5;
        eledef.nodes.clear();

        eletypes_exp[HEX8] = eledef;
        eletypes_imp[5]    = eledef;
      }

      // HEX20
      {
        eledef.type    = HEX20;
        eledef.dim     = 3;
        eledef.nnodes  = 20;
        eledef.exptype = 101;
        const unsigned int nodes[] = {1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,16};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[HEX20] = eledef;
        eletypes_imp[12]    = eledef;
      }

      // HEX27
      {
        eledef.type    = HEX27;
        eledef.dim     = 3;
        eledef.nnodes  = 27;
        eledef.exptype = 12;
        const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,
            15,16,19,17,18,20,21,24,22,23,25,26};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[HEX27] = eledef;
        eletypes_imp[12]    = eledef;
      }

      // TET4
      {
        eledef.type    = TET4;
        eledef.dim     = 3;
        eledef.nnodes  = 4;
        eledef.exptype = 4;
        eledef.nodes.clear();

        eletypes_exp[TET4] = eledef;
        eletypes_imp[4]    = eledef;
      }

      // TET10
      {
        eledef.type    = TET10;
        eledef.dim     = 3;
        eledef.nnodes  = 10;
        eledef.exptype = 11;
        const unsigned int nodes[] = {0,1,2,3,4,5,6,7,9,8};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);
        eletypes_exp[TET10] = eledef;
        eletypes_imp[11]    = eledef;
      }

      // PRISM6
      {
        eledef.type    = PRISM6;
        eledef.dim     = 3;
        eledef.nnodes  = 6;
        eledef.exptype = 6;
        eledef.nodes.clear();

        eletypes_exp[PRISM6] = eledef;
        eletypes_imp[6]      = eledef;
      }

      // PRISM15
      {
        eledef.type    = PRISM15;
        eledef.dim     = 3;
        eledef.nnodes  = 15;
        eledef.exptype = 103;
        eledef.nodes.clear();

        eletypes_exp[PRISM15] = eledef;
        eletypes_imp[13] = eledef;
      }

      // PRISM18
      {
        eledef.type    = PRISM18;
        eledef.dim     = 3;
        eledef.nnodes  = 18;
        eledef.exptype = 13;
        const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,
            12,14,13,15,17,16};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[PRISM18] = eledef;
        eletypes_imp[13]      = eledef;
      }

      // PYRAMID5
      {
        eledef.type    = PYRAMID5;
        eledef.dim     = 3;
        eledef.nnodes  = 5;
        eledef.exptype = 7;
        eledef.nodes.clear();

        eletypes_exp[PYRAMID5] = eledef;
        eletypes_imp[7]        = eledef;
      }

      //==============================
    }
  }

} // end anonymous namespace


// ReadGMSH  members
void ReadGMSH::read(const std::string& name)
{
  std::ifstream in(name.c_str());

  if (!in.good())
    throw InitFailedException(std::string("Cannot read mesh file (") + name + ").");

  this->read_mesh(in);
}


void ReadGMSH::read_mesh(std::istream& in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert(libMesh::processor_id() == 0);

  libmesh_assert(in.good());

  // initialize the map with element types
  init_eletypes();

  // clear any data in the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();
  _reg_info.clear();
  _bd_regions.clear();


  // we will try to get the mesh dimension from the file
  unsigned int dim = mesh.mesh_dimension();

  // some variables
  const int  bufLen = 256;
  char       buf[bufLen+1];
  int        format=0, size=0;
  Real       version = 1.0;

  // map to hold the node numbers for translation
  // note the the nodes can be non-consecutive
  std::map<unsigned int, unsigned int> nodetrans;

  // map to hold the physical names and dimensions (if found)
  std::map<unsigned int, std::string> phys_names;

  {
    while (!in.eof())
    {
      in >> buf;

      if (!std::strncmp(buf,"$MeshFormat",11))
      {
        in >> version >> format >> size;

        if (version > 2.1)
          throw InitFailedException("Unsupported msh file version.");

        if (format)
          throw InitFailedException("Unknown data format for mesh.");
      }

      else if (!std::strncmp(buf, "$PhysicalNames", 14))
      {
        // we get the physical names and try to guess the mesh dimension
        int num_phys_names;
        in >> num_phys_names;

        for (int i = 0; i < num_phys_names; i++)
        {
          int id, d = -1;
          std::string name;
          if (version >= 2.1)
          {
            in >> d >> id >> name;
            dim = (d > dim) ? d : dim;
          }
          else
            in >> id >> name;

          // name is double quoted!
          name.erase(0, 1);
          name.erase(name.size() - 1, 1);

          phys_names[id] = name;
        }

        if (version >= 2.1)
          mesh.set_mesh_dimension(dim);

      }

      // read the node block
      else if (!std::strncmp(buf,"$NOD",4) ||
          !std::strncmp(buf,"$NOE",4) ||
          !std::strncmp(buf,"$Nodes",6))
      {

        // check the dimension for reasonable value
        if ((dim < 1) || (dim > 3))
        {
          std::ostringstream os;
          os << "mesh dimension of " << dim << " is invalid." << std::endl
            << "Hint: check the option \'dimension\' in the "
            << "device options block.";
          throw InitFailedException(os.str());
        }

        unsigned int numNodes = 0;
        in >> numNodes;
        mesh.reserve_nodes (numNodes);

        // read in the nodal coordinates and form points.
        Real x, y, z;
        unsigned int id;

        // add the nodal coordinates to the mesh
        for (unsigned int i = 0; i < numNodes; ++i)
        {
          in >> id >> x >> y >> z;
          mesh.add_point (Point(x, y, z), i);
          nodetrans[id] = i;
        }
        // read the $ENDNOD delimiter
        in >> buf;
      }

      /**
       * Read the element block
       *
       * If the element dimension is smaller than the mesh dimension, it
       * is added to the BoundaryRegions object.
       *
       * Because the elements might not yet exist, the sides are put on hold
       * until the elements are created, and inserted once reading elements is
       * finished
       */
      else if (!std::strncmp(buf,"$ELM",4) ||
          !std::strncmp(buf,"$Elements",9))
      {
        unsigned int numElem = 0;
        std::vector< boundaryElementInfo > boundary_elem;
        std::vector< boundaryElementInfo > edge_elem;

        // read how many elements are there, and reserve space in the mesh
        in >> numElem;
        mesh.reserve_elem(numElem);

        // read the elements
        unsigned int elem_id_counter = 0;
        for (unsigned int iel = 0; iel < numElem; ++iel)
        {
          unsigned int id, type, physical, elementary,
          /* partition = 1,*/ nnodes, ntags;
          // note - partition was assigned but never used - BSK
          if(version <= 1.0)
          {
            in >> id >> type >> physical >> elementary >> nnodes;
          }
          else
          {
            in >> id >> type >> ntags;
            elementary = physical = /* partition = */ 1;
            for(unsigned int j = 0; j < ntags; j++)
            {
              int tag;
              in >> tag;
              if(j == 0)
                physical = tag;
              else if(j == 1)
                elementary = tag;
              // else if(j == 2)
              //  partition = tag;
              // ignore any other tags for now
            }
          }

          // consult the import element table which element to build
          const elementDefinition& eletype = eletypes_imp[type];
          nnodes = eletype.nnodes;

          // only elements that match the mesh dimension are added
          // if the element dimension is less than dim, the nodes and
          // sides are added to the BoundaryRegions
          if (eletype.dim == dim)
          {
            // add the elements to the mesh
            Elem* elem = Elem::build(eletype.type).release();
            elem->set_id(elem_id_counter);
            mesh.add_elem(elem);

            // this is different from iel: lower dimensional elems aren't added
            elem_id_counter++;

            // check number of nodes. We cannot do that for version 2.0
            if (version <= 1.0)
            {
              if (elem->n_nodes() != nnodes)
              {
                std::cerr << "Number of nodes for element " << id
                << " of type " << eletypes_imp[type].type
                << " (Gmsh type " << type
                << ") does not match Libmesh definition. "
                << "I expected " << elem->n_nodes()
                << " nodes, but got " << nnodes << "\n";
                libmesh_error();
              }
            }

            // add node pointers to the elements
            int nod = 0;
            // if there is a node translation table, use it
            if (eletype.nodes.size() > 0)
              for (unsigned int i=0; i<nnodes; i++)
              {
                in >> nod;
                elem->set_node(eletype.nodes[i]) = mesh.node_ptr(nodetrans[nod]);
              }
            else
            {
              for (unsigned int i = 0; i < nnodes; i++)
              {
                in >> nod;
                elem->set_node(i) = mesh.node_ptr(nodetrans[nod]);
              }
            }

            // Finally, set the subdomain ID to physical
            elem->subdomain_id() = static_cast<subdomain_id_type>(physical);

            // add the subdomain id to the MeshRegionInfo
            _reg_info.add_id(elem->subdomain_id());

          } // if element.dim == dim

          else if (eletype.dim == dim-1)
          {
            // this is a boundary
            /**
             * add the boundary element nodes to the set of nodes
             */

            boundaryElementInfo binfo;
            std::set<unsigned int>::iterator iter = binfo.nodes.begin();
            int nod = 0;
            for (unsigned int i = 0; i < nnodes; i++)
            {
              in >> nod;
              binfo.nodes.insert(iter, nodetrans[nod]);
            }
            binfo.id = physical;
            boundary_elem.push_back(binfo);
          }

          else if (eletype.dim == dim-2)
          {
            // this is an edge
            boundaryElementInfo binfo;
            std::set<unsigned int>::iterator iter = binfo.nodes.begin();
            int nod = 0;
            for (unsigned int i = 0; i < nnodes; i++)
            {
              in >> nod;
              binfo.nodes.insert(iter, nodetrans[nod]);
            }
            binfo.id = physical;
            edge_elem.push_back(binfo);
          }

          else if (eletype.dim == dim-3)
          {
            // this is a node (we get here only in 3D)

            int nod = 0;
            for (unsigned int i = 0; i < nnodes; i++)
            {
              in >> nod;
              _bd_regions.add_node(mesh.node_ptr(nodetrans[nod]), physical);
            }
          }

          else
          {
            // this means eletype.dim > dim and is an error
            std::ostringstream os;
            os << "Trying to add a " << eletype.dim << "D element "
              << "into a " << dim << "D mesh! " << std::endl
              << "Hint: check the option \'dimension\' in the "
              << "device options block.";
            throw InitFailedException(os.str());
          }
        }//element loop

        // read the $ENDELM delimiter
        in >> buf;

        // check for consistency regarding mesh dimension
        // elem_id_counter has to be > 0, otherwise mesh dim was assumed too big
        if (elem_id_counter == 0)
        {
          std::ostringstream os;
          os << "The expected mesh dimension of " << dim
            << " seems to be bigger than the actual one. " << std::endl
            << "Hint: check the option \'dimension\' in the "
            << "device options block.";
            throw InitFailedException(os.str());
        }

        /**
         * If any lower dimensional elements have been found in the file,
         * try to add them to the mesh.boundary_info as sides and nodes with
         * the respecitve id's (called "physical" in Gmsh).
         */
        if (boundary_elem.size() > 0)
        {
          // create a index of the boundary nodes to easily locate which
          // element might have that boundary
          TiberCad::HashMap<unsigned int, std::vector<unsigned int> >::Type node_index;
          for (unsigned int i = 0; i < boundary_elem.size(); i++)
          {
            boundaryElementInfo binfo = boundary_elem[i];
            std::set<unsigned int>::iterator iter = binfo.nodes.begin();
            for ( ; iter != binfo.nodes.end(); ++iter)
              node_index[*iter].push_back(i);
          }


          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          // iterate over all elements and see which boundary element has
          // the same set of nodes as one of the boundary elements previously read
          for ( ; it != end; ++it)
          {
            const Elem* elem = *it;

            for (unsigned int s = 0; s < elem->n_sides(); s++)
            //if (elem->neighbor(s) == NULL)
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
                  if (boundary_elem[bidx].nodes == side_nodes)
                  {
                    unsigned int id = boundary_elem[bidx].id;
                    _bd_regions.add_side(elem, s, id);
                  }
                }
              }
            } // if elem->neighbor(s) == NULL
          } // element loop
        } // if boundary_elem.size() > 0


        // We do the same for the edges
        if (edge_elem.size() > 0)
        {
          // create a index of the boundary nodes to easily locate which
          TiberCad::HashMap<unsigned int, std::vector<unsigned int> >::Type node_index;
          for (unsigned int i = 0; i < edge_elem.size(); i++)
          {
            boundaryElementInfo binfo = edge_elem[i];
            std::set<unsigned int>::iterator iter = binfo.nodes.begin();
            for ( ; iter != binfo.nodes.end(); ++iter)
              node_index[*iter].push_back(i);
          }

          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          // iterate over all elements and see which boundary element has
          // the same set of nodes as one of the boundary elements previously read
          for ( ; it != end; ++it)
          {
            const Elem* elem = *it;

            for (unsigned int s = 0; s < elem->n_edges(); s++)
            {
              AutoPtr<Elem> side (elem->build_edge(s));
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
                  if (edge_elem[bidx].nodes == side_nodes)
                  {
                    unsigned int id = edge_elem[bidx].id;
                    _bd_regions.add_edge(elem, s, id);
                  }
                }
              }
            } // if elem->neighbor(s) == NULL
          } // element loop
        } // if boundary_elem.size() > 0
      } // if $ELM

    } // while !in.eof()

  }

  // set the physical region names
  std::map<unsigned int, std::string>::iterator it = phys_names.begin();
  for ( ; it != phys_names.end(); ++it)
  {
    ID id = it->first;

    _reg_info.set_name(id, it->second);
    _bd_regions.set_name(id, it->second);
  }


}




