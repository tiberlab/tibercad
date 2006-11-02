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


/** ISE Reading and Libmesh Writing Class.
 *  Reads a '*.grd' (ISE T-Cad) file and converts it into
 * libmesh type ('*.xda'; '*.xta'). Extracts also boundary conditions.
 */
class Read_ISE_Grid
{
		
 public:


/** Constructor. Needs the '*.grd' file path.
 *  Writes *.xda and *.xta files into default directory.
 *  Gets Boundary conditions.
 */
  Read_ISE_Grid(const char* file_name);


/** Virtual Destructor. Dummy. */
  virtual ~Read_ISE_Grid();
	

/** Writes libmesh mesh file in DEAL format (xda).
 *  *.xda file contains DEAL header, list of all simulation-dimension mesh elements
 *  and list of all nodes.
 */
  void write_xda();


/** Writes meshdata data file (xta).
 *  *.xta file contains list of all simulation-dimension elements
 *  associated to physical region ID.
 */
  void write_xta();


/** Deletes repetitions in node list. */
  void  unique_nodes(vector<unsigned int>& v1);


/** Writes Boundary region nodes map [(dimension-1) elements] and relative region. */
  void  get_BC_data (map<unsigned int, vector<unsigned int> >& BoundCond_map );


/** Xda file name string. */  
  string fname_xda;


/** Xta file name string. */
  string fname_xta;


/** Returns simulation dimension. */
  unsigned int get_dim();


  private:


/** Sets Boundary Condition Data. */
  void set_BC_data();


/** Integrity Check. Controls if file Version and Type are correct.
 *  Otherwise aborts program.
 */
  void integrity_check(int ver, string tp);


/** Reads Input file. */
  void scan_grid_file();
	

/** List of all ISE file types. */
  vector<unsigned int> ISE_element_type_list;


/** Elements belonging to xda file. */
  vector<ISE_Element*> xda_list_elements;

	
/** Base variables. */
  unsigned int dimension, nb_vertices, nb_edges, nb_faces, nb_elements, nb_regions;
	

/** Vector of all vertices. */
  vector<ISE_Vertex*> vertices;


/** Input file name. */
  const char* ISE_file_name;
	

/** Temporary vertex. */
  ISE_Vertex* vertex_point;
		

/** Temporary edge. */
  ISE_Edge* edge_point;


/** Temporary edge vector. */
  vector<ISE_Edge*> edges;


/** Temporary face. */
  ISE_Face* face_point;


/** Temporary face vector. */
  vector<ISE_Face*> faces;


/** Temporary Element. */
  ISE_Element* elements_list_point;


/** List of all elements. */
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
	
  //map <unsigned int tiber_BC_region, vector<unsigned int> region_nodes>
  map <unsigned int , vector<unsigned int>  > map_BC_region_nodes;
	
	
};

inline unsigned int Read_ISE_Grid::get_dim()
{
	return (dimension);
};

#endif /*_READ_ISE_GRID_H_*/
