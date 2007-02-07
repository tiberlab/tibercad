#include "MeshInput.h"
#include "ReadISEGrid.h"
#include "Read_MSH.h"

#include <iomanip>
#include <vector>
#include <iostream> //for I/O interaction
#include <fstream>  //for file streaming
#include <string>   //for strings
#include <set>
#include <map>
using namespace std;

MeshInput::MeshInput(void)
{

}

MeshInput::~MeshInput(void)
{

}


void MeshInput::read_mesh(const std::string& file_name,unsigned int sim_dim,
                          Mesh& mesh, MeshData_elements& mesh_data,
                          map<unsigned int, vector<unsigned int> >& BoundCond  )
{

  // See if the file exists.  
  
  std::ifstream in (file_name.c_str());
     
  if (!in.good())
  {
    std::cerr << "ERROR: cannot locate specified mesh file:\n\t"
              << file_name
              << std::endl;
    error();
  }

  if  ( file_name.rfind(".grd") < file_name.size() )
  {

    ReadISEGrid  ISE_mesh( file_name.c_str() , mesh, mesh_data );
    ISE_mesh.get_BC_data (BoundCond );
    //    BoundCond =  BoundCond_map;

  }

  else if ( file_name.rfind(".msh") < file_name.size() )
  {
    Read_MSH  GMSH_mesh( file_name,sim_dim, mesh, mesh_data );
    GMSH_mesh.get_BC_data (BoundCond );
  }



}
