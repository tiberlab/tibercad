#ifndef _READISEGRID_H_
#define _READISEGRID_H_


#include <iomanip>

#include <iostream> //for I/O interaction
#include <fstream>  //for file streaming
#include <InitFailedException.h>

/* #include <string>   //for strings */
 /* #include <set> */
 /* #include <map> */

#include "ISE_Vertex.h"
#include "ISE_Edge.h"
#include "ISE_Face.h"
#include "ISE_Element_3D.h"
#include "ISE_Element_2D.h"
#include "ISE_Element_1D.h"
#include "ISE_Element_0D.h"
#include "TypeDefs.h"

//  LibMesh  include files
#include "mesh_data_elements.h"
#include "mesh_data.h"
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"

//!ISE Reading and Libmesh Writing Class.
/*!
 *  Reads a '*.grd' (ISE T-Cad) file and converts it into
 * libmesh type ('*.xda'; '*.xta'). Extracts also boundary conditions.
 */
class ReadISEGrid 
{
		
 public:
  //!  Constructor
  /*!  Needs the '*.grd' file path, the \c Mesh and \c MeshData objects.
   *  Writes *.xda and *.xta files into default directory.
   * 
   */
 
  //  ReadISEGrid(const char* file_name , Mesh& mesh, MeshData_elements&  mesh_data );
  ReadISEGrid(const char* file_name);


  //! Virtual Destructor.
  /*! 
    Deallocate pointers.
  */
  virtual ~ReadISEGrid();
	

  //!  Gets Boundary region nodes map [(dimension-1) elements].
  /*!
   *Puts in map BoundCond_map  a  map which associates   each Boundary region ID  with its nodes.
   */ 
  void  get_BC_data (std::map<unsigned int, std::vector<unsigned int> >& BoundCond_map );
 

  /*!
    Returns simulation dimension.
  */
  unsigned int get_dim();


  /*!
    Gets the map which associates internal mesh IDs to user-defined physical region names
  */
  void  get_region_names_map (std::map<ID, std::string >& region_names_map );

  /*!
    Gets the map which associates internal mesh BC IDs to user-defined BC  region names
  */
  void  get_BC_region_names_map (std::map<ID, std::string >& region_names_map );





 private:
  

  //!  Writes  mesh file for  Libmesh in  DEAL  format  (.xda)
  /*!
   *  .xda file contains a DEAL header, a list of all (simulation-dimension) mesh elements
   *  and a list of all the nodes.
   */
  void write_xda();

  

  //! Writes meshdata data file (xta).
  /*!
   *  .xta file contains a list of all simulation-dimension elements
   *  associated to the relative physical region ID.
   */
  void write_xta();

  //! Utility    to  eliminate  repetitions  in  node  list.
  /*!
    Deletes repetitions in node list.
  */
  void  unique_nodes(std::vector<unsigned int>& v1);

 


  /*!
    Xda file name string.
  */  
  std::string fname_xda;

  /*!
    Xta file name string.
  */
  std::string fname_xta;

 




  /*!
    Write  mesh  and  meshdata from  .xda and  .xta  files
  */
  void   read_mesh_and_data(Mesh& mesh, MeshData_elements&  mesh_data );

  /*!
    Sets Boundary Condition Data.
  */
  void set_BC_data();

  /*!
   *  Integrity Check. Controls if file Version and Type are correct.
   *  Otherwise aborts program.
   */
  void integrity_check(int ver, std::string tp);

  /*!
    Reads Input file.
  */
  void scan_grid_file();
	
  /*!
    List of all ISE file types whose  dimension is equal to the simulation dimension.
  */
  std::vector<unsigned int> ISE_element_type_list;

  /*!
    Elements belonging to xda file (subset of the vector  of all elements "elements_list").
  */
  std::vector<ISE_Element*> xda_list_elements;  
	
  /*!
    Base variables.
  */
  unsigned int dimension, nb_vertices, nb_edges, nb_faces, nb_elements, nb_regions;


  /*!
    Vector of all vertices.
  */	
  std::vector<ISE_Vertex*> vertices; 

  /*!
    Input file name.
  */
  const char* ISE_file_name;
	
	


 

  /*!
    Temporary edge vector.
  */
  std::vector<ISE_Edge*> edges;

 

  /*!
    Temporary face vector.
  */
  std::vector<ISE_Face*> faces;

 

  /*!
    List of all elements.
  */
  std::vector<ISE_Element*> elements_list;  //  delete  these  pointers in  destructor !!
	
  /*!
    lists of ISE physical region IDs, respectively of 0D, 1D, 2D and 3D.
  */
  std::vector<unsigned int> regions_0D;
  std::vector<unsigned int> regions_1D;
  std::vector<unsigned int> regions_2D;
  std::vector<unsigned int> regions_3D;

 

  /*!	
    maps <ISE phis reg ID, region_elements_nD>
  */
  std::map <unsigned int , vector<ISE_Element*> >  map_0D_region_elements;
  std::map <unsigned int , vector<ISE_Element*> >  map_1D_region_elements;
  std::map <unsigned int , vector<ISE_Element*> >  map_2D_region_elements;
  std::map <unsigned int , vector<ISE_Element*> >  map_3D_region_elements;
  /*!	
    map <unsigned int tiber_BC_region, vector<unsigned int> region_nodes )
  */
  std::map <unsigned int , vector<unsigned int>  > map_BC_region_nodes;
	


  /*!	
    Map which associates internal mesh IDs to user-defined physical region names 
  */
  std::map<ID, std::string >  region_names_ID_map;


  /*!	
    Map which associates internal mesh BC IDs to user-defined BC region names 
  */
  std::map<ID, std::string > BC_region_names_ID_map;



	
};

inline unsigned int ReadISEGrid::get_dim()
{
  return (dimension);
}

#endif /* _READISEGRID_H_ */
