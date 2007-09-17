// $Id$

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

  unsigned int n_nodes = mesh.n_nodes();


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

    out << "5\n";
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

        node_map[global_id] = value;
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

  unsigned int n_nodes = mesh.n_nodes();


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

    out << "5\n";
  }
  

      
  // get ordered nodal data using a map
  typedef std::pair<unsigned int, Number> key_value_pair;
  typedef std::map<unsigned int, Number> map_type;
  typedef map_type::iterator map_iterator;

  map_type node_map;


  out << "\n";
  out << "CELL_DATA " << n_active_elem << "\n";

  for (unsigned int var = 0; var < n_vars; var++)
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
  out << "\n";

 
}
