#ifndef _READ_ISE_GRID_H_
#define _READ_ISE_GRID_H_

#include <iomanip>

#include <iostream> //for I/O interaction
#include <fstream>  //for file streaming
#include <string>   //for strings
#include <set>
#include <map>

#include "ISE_Vertex.h"
#include "ISE_Edge.h"
#include "ISE_Face.h"
#include "ISE_Element_3D.h"
#include "ISE_Element_2D.h"
#include "ISE_Element_1D.h"
#include "ISE_Element_0D.h"
class Read_ISE_Grid
{
		
 public:

  Read_ISE_Grid(const char* file_name);

  virtual ~Read_ISE_Grid();
	
  //! Writes  mesh file for  Libmesh in  DEAL  format  (*.xda)
  //Writes .xda file with a  DEAL header, a  list  of  all (simulation-dimension) elements  
  //in the  mesh and  a  list  of all  nodes.

  void write_xda();

  //! Writes  data file for  meshdata (*.xta) 
  // Writes .xta file with a  list  of  all  (simulation-dimension) elements 
  // associated to  related physical  region ID

  void write_xta();

//! Utility    to  eliminate  repetitions  in  node  list

  void  unique_nodes(vector<unsigned int>& v1);

  //!  Writes map which associates  nodes belonging to Boundary regions
  // [ (dimension-1) elements] and relative regions


  void  get_BC_data (map<unsigned int, vector<unsigned int> >& BoundCond_map );


  
  string fname_xda;

  string fname_xta;


  unsigned int get_dim();


  private:


  void set_BC_data();

  void integrity_check(int ver, string tp);

  void scan_grid_file();
	
  vector<unsigned int> ISE_element_type_list;

  vector<ISE_Element*> xda_list_elements;
	
  unsigned int dimension, nb_vertices, nb_edges, nb_faces, nb_elements, nb_regions;
	
  vector<ISE_Vertex*> vertices;

  const char* ISE_file_name;
	

	


  ISE_Vertex* vertex_point;
		
  ISE_Edge* edge_point;
  vector<ISE_Edge*> edges;

  ISE_Face* face_point;
  vector<ISE_Face*> faces;

  ISE_Element* elements_list_point;
  vector<ISE_Element*> elements_list;
	

  // lists of ISE physical region IDs, respectively of 0D, 1D, 2D and 3D
  vector<unsigned int> regions_0D;
  vector<unsigned int> regions_1D;
  vector<unsigned int> regions_2D;
  vector<unsigned int> regions_3D;

  // lists of nD Elements with nD region
  vector<ISE_Element*>  region_elements_0D;
  vector<ISE_Element*>  region_elements_1D;
  vector<ISE_Element*>  region_elements_2D;
  vector<ISE_Element*>  region_elements_3D;

	
  //  maps <ISE phis reg ID, region_elements_nD>
  map <unsigned int , vector<ISE_Element*> >  map_0D_region_elements;
  map <unsigned int , vector<ISE_Element*> >  map_1D_region_elements;
  map <unsigned int , vector<ISE_Element*> >  map_2D_region_elements;
  map <unsigned int , vector<ISE_Element*> >  map_3D_region_elements;
	
  //map <unsigned int tiber_BC_region, vector<unsigned int> region_nodes )
  map <unsigned int , vector<unsigned int>  > map_BC_region_nodes;
	
	
};

inline unsigned int Read_ISE_Grid::get_dim()
{
	return (dimension);
};

#endif /*_READ_ISE_GRID_H_*/
