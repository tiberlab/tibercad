// $Id$

#include "MeshInput.h"
#include "ReadISEGrid.h"
#include "Read_MSH.h"

#include "mesh.h"
#include "mesh_data_elements.h"

#include <iostream> //for I/O interaction
#include <fstream>  //for file streaming

using namespace std;



// void MeshInput::read_mesh(const string& file_name,unsigned int sim_dim,
//                           Mesh& mesh, MeshData_elements& mesh_data,
//                           map<unsigned int, vector<unsigned int> >& BoundCond  )

void MeshInput::read_mesh(const string& file_name,unsigned int sim_dim,
                          Mesh& mesh, MeshData_elements& mesh_data,
                          map<unsigned int, vector<unsigned int> >& BoundCond ,
                          map<ID, string >& region_names_map,
                          map<ID, string >& BC_region_names_map)


{

  // See if the file exists.  
  
  ifstream in (file_name.c_str());
     
  if (!in.good())
  {
    cerr << "ERROR: cannot locate specified mesh file:\n\t"
              << file_name
              << endl;
    error();
  }

  if  ( file_name.rfind(".grd") < file_name.size() )
  {

    ReadISEGrid  ISE_mesh( file_name.c_str() , mesh, mesh_data );
    ISE_mesh.get_BC_data (BoundCond );
    ISE_mesh.get_region_names_map (region_names_map );  //  not empty only  for  ISE  grid
    ISE_mesh.get_BC_region_names_map (BC_region_names_map ); //  not empty only  for  ISE  grid

  }

  else if ( file_name.rfind(".msh") < file_name.size() )
  {
    Read_MSH  GMSH_mesh( file_name,sim_dim, mesh, mesh_data );
    GMSH_mesh.get_BC_data (BoundCond );
  }



}
