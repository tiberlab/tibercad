// $Id$

#include "MeshInput.h"
#include "ReadISEGrid.h"
#include "Read_MSH.h"
#include "InitFailedException.h"

#include "mesh.h"
#include "mesh_data_elements.h"

#include <iostream> //for I/O interaction
#include <fstream>  //for file streaming

using namespace std;



void MeshInput::read_mesh(const string& file_name, int sim_dim,
                          Mesh*& mesh, MeshData_elements*& mesh_data,
                          map<unsigned int, vector<unsigned int> >& BoundCond ,
                          map<ID, string >& region_names_map,
                          map<ID, string >& BC_region_names_map)

{

  unsigned int   dimension;

  // See if the file exists.  
  ifstream in (file_name.c_str());
     
  if (!in.good())
  {
    ostringstream os;
    os << "Cannot locate specified mesh file: " << file_name;
    throw InitFailedException(os.str());
  }

  if  ( file_name.rfind(".grd") < file_name.size() )
  {

    //   ReadISEGrid ISE_mesh(file_name.c_str(), mesh, mesh_data);
    ReadISEGrid ISE_mesh(file_name.c_str());


    ISE_mesh.get_BC_data(BoundCond);
    ISE_mesh.get_region_names_map(region_names_map);
    ISE_mesh.get_BC_region_names_map(BC_region_names_map);
    //
    dimension  = ISE_mesh.get_dim();

    mesh = new Mesh(dimension);

    mesh_data = new MeshData_elements(*mesh);
    string  mesh_file_data = "elem_data.xta";
    string  mesh_file_inp = "in.xda";

    (*mesh_data).enable_compatibility_mode();

    //  mesh->read (mesh_file_inp,&mesh_data  ); 
    mesh->read (mesh_file_inp,mesh_data  ); 

  
    mesh_data->read(mesh_file_data);

    // ----------------------------------------------


  }
  else if ( file_name.rfind(".msh") < file_name.size() )
  {
    //  Read_MSH GMSH_mesh(file_name, sim_dim, mesh, mesh_data);
    Read_MSH GMSH_mesh(file_name,sim_dim);


    GMSH_mesh.get_BC_data(BoundCond);
    dimension = GMSH_mesh.get_simulation_dim();


    // NOTE: new   for  sim_dim read  from  mesh  file !!!!
    // ----------------------------------------------
   
    mesh = new Mesh(dimension);
    mesh_data = new MeshData_elements(*mesh);
    string  mesh_file_data = "elem_data.xta";
    string  mesh_file_inp = "in.xda";

    (*mesh_data).enable_compatibility_mode();

    //  mesh->read (mesh_file_inp,&mesh_data  ); 
    mesh->read (mesh_file_inp,mesh_data  ); 

  
    mesh_data->read(mesh_file_data);

    // ----------------------------------------------



    // assign the names
    map<ID, string> names;
    GMSH_mesh.get_physical_names_map(names);
    if (names.size() > 0)
    {
      map<ID, string>::const_iterator it(names.begin());
      const map<ID, string>::const_iterator end(names.end());
      for ( ; it != end; ++it)
      {
        if (BoundCond.find(it->first) != BoundCond.end())
          BC_region_names_map[it->first] = it->second;
        else
          region_names_map[it->first] = it->second;
      }
    }
  }
}
