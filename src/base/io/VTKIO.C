// $Id$

#include "boost/tokenizer.hpp"

#include "VTKIO.h"
#include "RuntimeException.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"

#include "b64/encode.h"



// C++ includes
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <typeinfo>

using namespace std;




TiberVTKIO::TiberVTKIO(const MeshBase& mesh)
  : DataOutput()
{
  set_mesh(mesh);
}


void TiberVTKIO::write_nodal_data(const string& fname,
    const vector<double>& soln,
    const vector<string>& names)
{

  const MeshBase& mesh = get_mesh();

  if(mesh.comm().rank() != 0)
    return;

  if (names.size() == 0)
    return;


  set<unsigned int> node_ids;

  // count nodes of active part of mesh ...
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
  ofstream out(fname.c_str());

  if (!out.good())
    throw runtime_error("Could not open " + fname);

  // The number of variables in the equation system
  const unsigned int n_vars = names.size();

  // Write header to stream
  out << "# vtk DataFile Version 2.0\n"
      << "Nodal data\n"
      << "ASCII\n\n";
  out << "DATASET UNSTRUCTURED_GRID\n";
  out << "POINTS " << n_nodes << " double\n";

  map<unsigned int, unsigned int> vtk_node_ids;

  set<unsigned int>::iterator nodeit(node_ids.begin());
  const set<unsigned int>::iterator nodeend(node_ids.end());
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
  typedef pair<unsigned int, Number> key_value_pair;
  typedef map<unsigned int, Number> map_type;
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

    out << setprecision(10);

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




void TiberVTKIO::write_elemental_data(const string& fname,
    const vector<double>& soln,
    const vector<string>& names)
{

  const MeshBase& mesh = get_mesh();

  if(mesh.comm().rank() != 0)
    return;

  if (names.size() == 0)
    return;


  set<unsigned int> node_ids;

  // count nodes  of active part of mesh ...
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
  ofstream out(fname.c_str());

  if (!out.good())
    throw runtime_error("Could not open " + fname);

  // The number of variables in the equation system
  const unsigned int n_vars = names.size();

  // Write header to stream
  out << "# vtk DataFile Version 2.0\n"
      << "Nodal data\n"
      << "ASCII\n\n";
  out << "DATASET UNSTRUCTURED_GRID\n";
  out << "POINTS " << n_nodes << " double\n";

  map<unsigned int, unsigned int> vtk_node_ids;

  set<unsigned int>::iterator nodeit(node_ids.begin());
  const set<unsigned int>::iterator nodeend(node_ids.end());
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
  typedef pair<unsigned int, Number> key_value_pair;
  typedef pair<unsigned int, vector<Number> > vector_key_value_pair;
  typedef map<unsigned int, Number> map_type;
  typedef map<unsigned int, vector<Number> > vector_map_type;
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

    bool is_vec = false;
    if (tokit_next != tokens.end())
    {
      // is it a vector?
      string coord(*tokit_next);

      if ((coord == "x") || (coord == "y") || (coord == "z"))
        is_vec = true;
    }

    if (is_vec)
    {
      // it is a vector
      string name(*tokit);
      string coord(*tokit_next);

      vector<int> indices(3, -1);


      if (coord == "x")
        indices[0] = var;
      else if (coord == "y")
        indices[1] = var;
      else if (coord == "z")
        indices[2] = var;

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
          string name2(*tokit);
          string coord2(*tokit_next);


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
                string name3(*tokit);
                string coord3(*tokit_next);

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

        vector<double> value(3, 0.0);
        if (indices[0] != -1)
          value[0] = soln[global_id + indices[0]];
        if (indices[1] != -1)
          value[1] = soln[global_id + indices[1]];
        if (indices[2] != -1)
          value[2] = soln[global_id + indices[2]];

        vector_node_map[elem_number] = value;

        elem_number++;
      }

      out << setprecision(10);

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

      out << setprecision(10);

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



TiberVTKIO::VTKCellType
TiberVTKIO::get_VTK_cell_type(const Elem* elem)
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


void
TiberVTKIO::do_write(bool force)
{
  const MeshBase& mesh = get_mesh();
  string file = get_filename() + ".vtu";


  if (mesh.comm().rank() == 0)
  {
    ofstream of(file.c_str());

    of << "<VTKFile type=\"UnstructuredGrid\" "
        << "version=\"0.1\" byte_order=\"";
    {
      short int word = 0x0001;
      char *byte = (char *) &word;
      if (byte[0])
        of << (byte[0] ? "LittleEndian" : "BigEndian");
    }
    of << "\">\n"
        << "<UnstructuredGrid>\n";

    of.flush();
    of.close();
  }


  for (unsigned int rank = 0; rank < mesh.comm().size(); ++rank)
  {

    if (mesh.comm().rank() == rank)
    {
      ofstream of(file.c_str(), std::ofstream::app);


      map<ID, vector<unsigned int> > nodes;
      map<ID, vector<VTKElem> > vtk_elem;
      create_pieces(nodes, vtk_elem);

      map<ID, vector<unsigned int> >::iterator it(nodes.begin());
      for ( ; it != nodes.end(); ++it)
      {
        ID id = it->first;

        // only pieces with data are written
        if (has_data(id) || force)
        {

          const vector<unsigned int>& nodevec = nodes[id];
          const vector<VTKElem>& elemvec = vtk_elem[id];
          unsigned int num_pt = nodevec.size();
          unsigned int num_el = elemvec.size();
          of << "<Piece NumberOfPoints=\"" << num_pt << "\" "
              << "NumberOfCells=\"" << num_el << "\">\n";

          vector<float> buffer;

          of << "<Points>\n";
          buffer.resize(3 * num_pt);
          for (unsigned int i = 0; i < num_pt; ++i)
          {
            const Point& p = mesh.point(nodevec[i]);
            buffer[3 * i] = static_cast<float>(p(0));
            buffer[3 * i + 1] = static_cast<float>(p(1));
            buffer[3 * i + 2] = static_cast<float>(p(2));
          }
          write_data_array("", 3, buffer, of);
          buffer.clear();
          of << "</Points>\n";

          vector<uint32_t> intbuf;
          of << "<Cells>\n";
          intbuf.reserve(4 * num_el);
          for (size_t i = 0; i < num_el; ++i)
          {
            const VTKElem& el = elemvec[i];
            size_t elconn = el.connectivity.size();
            if (intbuf.capacity() < (intbuf.size() + elconn))
              intbuf.reserve(2 * intbuf.capacity());

            for (size_t n = 0; n < el.connectivity.size(); n++)
              intbuf.push_back(static_cast<uint32_t>(el.connectivity[n]));
          }
          write_data_array("connectivity", 1, intbuf, of);

          intbuf.resize(num_el);
          size_t offset = 0;
          for (size_t i = 0; i < num_el; ++i)
          {
            const VTKElem& el = elemvec[i];
            offset += el.connectivity.size();
            intbuf[i] = static_cast<uint32_t>(offset);
          }
          write_data_array("offsets", 1, intbuf, of);
          intbuf.clear();


          vector<uint8_t> shortbuf(num_el);
          for (size_t i = 0; i < num_el; ++i)
          {
            const VTKElem& el = elemvec[i];
            shortbuf[i] = static_cast<uint8_t>(el.type);
          }
          write_data_array("types", 1, shortbuf, of);
          intbuf.clear();
          of << "</Cells>\n";

          if (has_data(id))
          {
            const DataMap& data = get_zone_data(id);

            of << "<PointData";
            DataMap::const_iterator it(data.begin());
            const DataMap::const_iterator end(data.end());

            bool scalar_ok = false;
            bool vector_ok = false;
            bool tensor_ok = false;
            for ( ; it != end; ++it)
            {
              const SolutionDescriptor& descr = it->first;
              if (descr.location() == SolutionDescriptor::NODES)
              {
                if (!scalar_ok && (descr.n_components() == 1))
                {
                  of << " Scalars=\"" << descr.name() << "\"";
                  scalar_ok = true;
                }
                else if (!vector_ok && (descr.n_components() == 3))
                {
                  of << " Vectors=\"" << descr.name() << "\"";
                  vector_ok = true;
                }
                else if (!tensor_ok && (descr.n_components() == 9))
                {
                  of << " Tensors=\"" << descr.name() << "\"";
                  tensor_ok = true;
                }
              }
            }
            of << ">\n";

            // loop over all data
            for (it = data.begin(); it != end; ++it)
            {
              const SolutionDescriptor& descr = it->first;
              if (descr.location() == SolutionDescriptor::NODES)
                write_data_array(descr.name(), descr.n_components(), it->second, of);
            }
            of << "</PointData>\n";

            of << "<CellData";
            scalar_ok = false;
            vector_ok = false;
            tensor_ok = false;
            for (it = data.begin(); it != end; ++it)
            {
              const SolutionDescriptor& descr = it->first;
              if (descr.location() == SolutionDescriptor::CELL)
              {
                if (!scalar_ok && (descr.n_components() == 1))
                {
                  of << " Scalars=\"" << descr.name() << "\"";
                  scalar_ok = true;
                }
                else if (!vector_ok && (descr.n_components() == 3))
                {
                  of << " Vectors=\"" << descr.name() << "\"";
                  vector_ok = true;
                }
                else if (!tensor_ok && (descr.n_components() == 9))
                {
                  of << " Tensors=\"" << descr.name() << "\"";
                  tensor_ok = true;
                }
              }
            }
            of << ">\n";
            for (it = data.begin(); it != end; ++it)
            {
              const SolutionDescriptor& descr = it->first;
              if (descr.location() == SolutionDescriptor::CELL)
                write_data_array(descr.name(), descr.n_components(), it->second, of);
            }
            of << "</CellData>\n";
          }

          of << "</Piece>\n";
        }
      }
      of.flush();
      of.close();
    }
    mesh.comm().barrier();
  }

  if (mesh.comm().rank() == 0)
  {
    ofstream of(file.c_str(), std::ofstream::app);
    of << "</UnstructuredGrid>\n"
        << "</VTKFile>\n";
  }

}


template <typename T>
void
TiberVTKIO::write_data_array(const string& name, int comp,
    const vector<T>& data, ostream& os)
{
  string type;
  int typelen;
  if ((typeid(T) == typeid(float)) || (typeid(T) == typeid(double)))
  {
    type = "Float32";
    typelen = 4;
  }
  else if (typeid(T) == typeid(uint32_t))
  {
    type = "UInt32";
    typelen = 4;
  }
  else if (typeid(T) == typeid(int32_t))
  {
    type = "Int32";
    typelen = 4;
  }
  else if (typeid(T) == typeid(uint8_t))
  {
    type = "UInt8";
    typelen = 1;
  }
  else if (typeid(T) == typeid(int8_t))
  {
    type = "Int8";
    typelen = 1;
  }
  else
  {
    ostringstream es;
    es << "VTK output is not implemented for data type \'"
        << typeid(T).name();
    throw RuntimeException(es.str());
  }

  os << "<DataArray type=\"" << type << "\" ";
  if (name.size() > 0)
    os << "Name=\"" << name << "\" ";
  os << "NumberOfComponents=\"" << comp << "\" "
      << "format=\"" << (is_ascii() ? "ascii" : "binary") << "\">\n";

  size_t n = data.size();

  if (is_ascii())
  {
    for (size_t i = 0; i < n; i++)
    {
      if ((typeid(T) == typeid(uint8_t)) || (typeid(T) == typeid(int8_t)))
        os << static_cast<int>(data[i]) << " ";
      else
        os << data[i] << " ";
    }
  }
  else
  {
    uint32_t nbytes = typelen * n;
    stringstream is(ios_base::in | ios_base::out | ios_base::binary);
    is.write(reinterpret_cast<char*>(&nbytes), 4);
    for (size_t i = 0; i < n; i++)
    {
      if (typeid(T) == typeid(double))
      {
        float tmp = static_cast<float>(data[i]);
        is.write(reinterpret_cast<const char*>(&tmp), typelen);
      }
      else
        is.write(reinterpret_cast<const char*>(&data[i]), typelen);
    }

    base64::encoder enc;
    enc.encode(is, os);
  }

  os << "\n</DataArray>\n";
}


void
TiberVTKIO::create_pieces(map<ID, vector<unsigned int> >& points,
    map<ID, vector<VTKElem> >& elems)
{
  std::set<libMesh::subdomain_id_type> subdomains;
  std::set<ID> subdomain_ids;
  this->get_zone_ids(subdomain_ids);
  for (auto&& id : subdomain_ids)
  {
    subdomains.insert(id);
  }

  {
    // # elements in each subdomain
    map<ID, unsigned int> n_elem;

    const MeshBase& mesh = get_mesh();

    MeshBase::const_element_iterator it  = mesh.active_local_subdomains_elements_begin(subdomains);
    MeshBase::const_element_iterator end = mesh.active_local_subdomains_elements_end(subdomains);

    for ( ; it != end; ++it)
    {
      const Elem* elem = *it;
      unsigned int id = elem->subdomain_id();
      map<ID, unsigned int>::iterator it(n_elem.find(id));
      if (it != n_elem.end())
        (it->second)++;
      else
        n_elem[id] = 0;
    }

    // translation table for node IDs
    map<ID, map<libMesh::dof_id_type, unsigned int> > ttable;

    for (it = mesh.active_local_subdomains_elements_begin(subdomains); it != end; ++it)
    {
      const Elem* elem = *it;
      unsigned int id = elem->subdomain_id();

      unsigned int nn = elem->n_nodes();
      for (unsigned int i = 0; i < nn; i++)
      {
        if (ttable[id].count(elem->node(i)) == 0)
        {
          unsigned int n_id = ttable[id].size();
          ttable[id][elem->node(i)] = n_id;
        }
      }
    }

    for (map<ID, unsigned int>::const_iterator i(n_elem.begin());
      i != n_elem.end(); ++i)
    {
      elems[i->first].reserve(i->second);
      points[i->first].resize(ttable[i->first].size());
    }


    for (it = mesh.active_local_subdomains_elements_begin(subdomains); it != end; ++it)
    {
      const Elem* elem = *it;
      unsigned int id = elem->subdomain_id();

      unsigned int nn = elem->n_nodes();

      VTKElem vtkel;// = elems[id][el_cnt];
      vtkel.type = get_VTK_cell_type(elem);
      vtkel.connectivity.resize(nn);

      for (unsigned int i = 0; i < nn; i++)
      {
        unsigned int n = ttable[id][elem->node(i)];
        points[id][n] = elem->node(i);
        vtkel.connectivity[i] = n;
      }

      elems[id].push_back(vtkel);
    }

  }
}
