// $Id$

#include "boost/tokenizer.hpp"

#include "VTKIO.h"

#include "elem.h"
#include "mesh_base.h"

// C++ includes
#include <fstream>
#include <sstream>
#include <map>
#include <stdexcept>

VTKIO::VTKIO(const MeshBase& mesh)
  : MeshOutput<MeshBase>(mesh)
{
}


void VTKIO::write_nodal_data(const std::string& fname,
    const std::vector<Number>& soln,
    const std::vector<std::string>& names)
{

  if(libMesh::processor_id() != 0)
    return;

  if (names.size() == 0)
    return;

  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();


  std::set<unsigned int> node_ids;

  // count nodes ...
  {
    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

    for ( ; it != end; ++it)
    {
      const Elem* elem = *it;

      unsigned int nn = elem->n_nodes();
      for (unsigned int i = 0; i < nn; i++)
        node_ids.insert(elem->node(i));
    }
  }

  unsigned int n_nodes = node_ids.size();



  // Create an output stream for script file
  std::ofstream out(fname.c_str());

  if (!out.good())
    throw std::runtime_error("Could not open " + fname);

  // The number of variables in the equation system
  const unsigned int n_vars = names.size();
  
  // Write header to stream
  out << "# vtk DataFile Version 2.0\n"
      << "Nodal data\n"
      << "ASCII\n\n";
  out << "DATASET UNSTRUCTURED_GRID\n";
  out << "POINTS " << n_nodes << " double\n";

  std::map<unsigned int, unsigned int> vtk_node_ids;
  
  std::set<unsigned int>::iterator nodeit(node_ids.begin());
  const std::set<unsigned int>::iterator nodeend(node_ids.end());
  for (unsigned int vtk_id = 0; nodeit != nodeend; ++nodeit, vtk_id++)
  {
    const Node& node = mesh.node(*nodeit);
    out << node(0) << " " << node(1) << " " << node(2) << "\n";

    vtk_node_ids[*nodeit] = vtk_id;
  }

  out << "\n";


  unsigned int n_active_elem = mesh.n_active_elem();

  unsigned int size = n_active_elem;


  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    size += elem->n_nodes();
  }

  out << "CELLS " << n_active_elem << " " << size << "\n";


  it  = mesh.active_elements_begin();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    unsigned int nn = elem->n_nodes();
    
    out << nn;

    for (unsigned int i = 0; i < nn; i++)
      out << " " << vtk_node_ids[elem->node(i)];

    out << "\n";
  }



  out << "\n";
  out << "CELL_TYPES " << n_active_elem << "\n";

  it  = mesh.active_elements_begin();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    out << get_VTK_cell_type(elem) << "\n";
  }
  

      
  // get ordered nodal data using a map
  typedef std::pair<unsigned int, Number> key_value_pair;
  typedef std::map<unsigned int, Number> map_type;
  typedef map_type::iterator map_iterator;

  map_type node_map;


  out << "\n";
  out << "POINT_DATA " << n_nodes << "\n";

  for (unsigned int var = 0; var < n_vars; var++)
  {
    out << "SCALARS " <<  names[var] << " double\n";
    out << "LOOKUP_TABLE default\n";


    it  = mesh.active_elements_begin();

    for ( ; it != end; ++it)
    {
      const Elem* elem = *it;

      for(unsigned int i = 0; i < elem->n_nodes(); i++)
      {
        Number value;

        // Get the global id of the node
        unsigned int global_id = elem->node(i);

        value = soln[global_id*n_vars + var];

        node_map[vtk_node_ids[elem->node(i)]] = value;
      }
    }

    out << std::setprecision(10);

    map_iterator map_it = node_map.begin();
    const map_iterator end_map_it = node_map.end();

    for( ; map_it != end_map_it; ++map_it)
    {
      key_value_pair kvp = *map_it;
      Number value = kvp.second;

      out << value << "\n";
    }
  }
  out << "\n";

  out.close();
}




void VTKIO::write_elemental_data(const std::string& fname,
    const std::vector<Number>& soln,
    const std::vector<std::string>& names)
{

  if(libMesh::processor_id() != 0)
    return;

  if (names.size() == 0)
    return;

  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  unsigned int n_nodes = 0;

  // count nodes ...
  {
    std::set<unsigned int> node_ids;
    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

    for ( ; it != end; ++it)
    {
      const Elem* elem = *it;

      unsigned int nn = elem->n_nodes();
      for (unsigned int i = 0; i < nn; i++)
        node_ids.insert(elem->node(i));
    }
    n_nodes = node_ids.size();
  }



  // Create an output stream for script file
  std::ofstream out(fname.c_str());

  if (!out.good())
    throw std::runtime_error("Could not open " + fname);

  // The number of variables in the equation system
  const unsigned int n_vars = names.size();
  
  // Write header to stream
  out << "# vtk DataFile Version 2.0\n"
      << "Nodal data\n"
      << "ASCII\n\n";
  out << "DATASET UNSTRUCTURED_GRID\n";
  out << "POINTS " << n_nodes << " double\n";

  
  for (unsigned int i = 0; i < n_nodes; i++)
  {
    const Node& node = mesh.node(i);
    out << node(0) << " " << node(1) << " " << node(2) << "\n";
  }

  out << "\n";


  unsigned int n_active_elem = mesh.n_active_elem();

  unsigned int size = n_active_elem;


  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    size += elem->n_nodes();
  }

  out << "CELLS " << n_active_elem << " " << size << "\n";


  it  = mesh.active_elements_begin();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    unsigned int nn = elem->n_nodes();
    
    out << nn;

    for (unsigned int i = 0; i < nn; i++)
      out << " " << elem->node(i);

    out << "\n";
  }



  out << "\n";
  out << "CELL_TYPES " << n_active_elem << "\n";

  it  = mesh.active_elements_begin();

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    out << get_VTK_cell_type(elem) << "\n";
  }
  

      
  // get ordered nodal data using a map
  typedef std::pair<unsigned int, Number> key_value_pair;
  typedef std::pair<unsigned int, std::vector<Number> > vector_key_value_pair;
  typedef std::map<unsigned int, Number> map_type;
  typedef std::map<unsigned int, std::vector<Number> > vector_map_type;
  typedef map_type::iterator map_iterator;
  typedef vector_map_type::iterator vector_map_iterator;

  map_type node_map;
  vector_map_type vector_node_map;


  out << "\n";
  out << "CELL_DATA " << n_active_elem << "\n";

  unsigned int var = 0;
  for ( ; var < n_vars; var++)
  {
    // check if it is a vector
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep("_");

    tokenizer tokens(names[var], sep);
    tokenizer::iterator tokit(tokens.begin());
    tokenizer::iterator tokit_next(tokit);
    ++tokit_next;

    if (tokit_next != tokens.end())
    {
      // it is a vector
      std::string name(*tokit);
      std::string coord(*tokit_next);

      std::vector<int> indices(3, -1);

      if (coord == "x")
        indices[0] = var;
      else if (coord == "y")
        indices[1] = var;
      else if (coord == "z")
        indices[2] = var;
      else
        continue;

      // check the following variable
      if ((var + 1) < n_vars)
      {
        tokens.assign(names[var+1]);
        tokit = tokens.begin();
        tokit_next = tokit;
        ++tokit_next;

        if (tokit_next != tokens.end())
        {
          // it is a vector
          std::string name2(*tokit);
          std::string coord2(*tokit_next);


          if (name2 == name)
          {
            if (coord2 == "x")
              indices[0] = var + 1;
            if (coord2 == "y")
              indices[1] = var + 1;
            else if (coord2 == "z")
              indices[2] = var + 1;
            else
              continue;

            var = var + 1;

            // check the following variable
            if ((var + 2) < n_vars)
            {
              tokens.assign(names[var+1]);
              tokit = tokens.begin();
              tokit_next = tokit;
              ++tokit_next;

              if (tokit_next != tokens.end())
              {
                std::string name3(*tokit);
                std::string coord3(*tokit_next);

                if (name3 == name)
                {
                  if (coord3 == "x")
                    indices[0] = var + 1;
                  else if (coord3 == "y")
                    indices[1] = var + 1;
                  else if (coord3 == "z")
                    indices[2] = var + 1;
                  else
                    continue;

                  var = var + 1;
                }
              }
            }
          }
        }
      }

      out << "VECTORS " <<  name << " double\n";

      it  = mesh.active_elements_begin();

      unsigned int elem_number = 0;
      for ( ; it != end; ++it)
      {
        const Elem* elem = *it;

        unsigned int global_id = elem_number * n_vars;

        std::vector<double> value(3, 0.0);
        if (indices[0] != -1)
          value[0] = soln[global_id + indices[0]];
        if (indices[1] != -1)
          value[1] = soln[global_id + indices[1]];
        if (indices[2] != -1)
          value[2] = soln[global_id + indices[2]];

        vector_node_map[global_id] = value;

        elem_number++;
      }

      out << std::setprecision(10);

      vector_map_iterator map_it = vector_node_map.begin();
      const vector_map_iterator end_map_it = vector_node_map.end();

      for( ; map_it != end_map_it; ++map_it)
      {
        vector_key_value_pair kvp = *map_it;

        out << kvp.second[0] << " "
          << kvp.second[1] << " "
          << kvp.second[2] << "\n";
      }

    }
    else
    {
    
      out << "SCALARS " <<  names[var] << " double\n";
      out << "LOOKUP_TABLE default\n";


      it  = mesh.active_elements_begin();

      unsigned int elem_number = 0;
      for ( ; it != end; ++it)
      {
        const Elem* elem = *it;

        unsigned int global_id = elem_number * n_vars;

        Number value = soln[global_id + var];

        node_map[global_id] = value;

        elem_number++;
      }

      out << std::setprecision(10);

      map_iterator map_it = node_map.begin();
      const map_iterator end_map_it = node_map.end();

      for( ; map_it != end_map_it; ++map_it)
      {
        key_value_pair kvp = *map_it;
        Number value = kvp.second;

        out << value << "\n";
      }
    }
  }
  out << "\n";

}



VTKIO::VTKCellType
VTKIO::get_VTK_cell_type(const Elem* elem)
{
  assert(elem != NULL);

  using namespace libMeshEnums;
  
  VTKCellType vtk_type;

  switch (elem->type())
  {
    case EDGE2:
      vtk_type = VTK_LINE;
      break;
    case EDGE3:
    case EDGE4:
      vtk_type = VTK_POLY_LINE;
      break;
    case TRI3:
      vtk_type = VTK_TRIANGLE;
      break;
    case TRI6:
    case QUAD8:
    case QUAD9:
      vtk_type = VTK_POLYGON;
      break;
    case QUAD4:
      vtk_type = VTK_QUAD;
      break;
    case TET4:
      vtk_type = VTK_TETRA;
      break;
    case HEX8:
      vtk_type = VTK_HEXAHEDRON;
      break;
    case PRISM6:
      vtk_type = VTK_WEDGE;
      break;
    case PYRAMID5:
      vtk_type = VTK_PYRAMID;
      break;
  }

  return vtk_type;
}
