#include "ReadISEGrid.h"

using namespace std;
//ReadISEGrid::ReadISEGrid(const char* file_name)
 ReadISEGrid::ReadISEGrid(const char* file_name , Mesh& mesh, MeshData_elements&  mesh_data )
{
	
  ISE_file_name = file_name;
  scan_grid_file();
  write_xda();
  write_xta();
  set_BC_data();
	
  // read .xta and  .xda files  in  Libmesh  data structures
  read_mesh_and_data(mesh,mesh_data );
}

ReadISEGrid::~ReadISEGrid()
{
}




void ReadISEGrid::scan_grid_file()
{
	
	 
	 
  ifstream ISE_INPUT(ISE_file_name);
	
  string dummy;
     
      
  do
  {
    ISE_INPUT >> dummy;
  }
  while (dummy != "Info");
    
  ////     READ  INFO BLOCK VARIABLES
  //    
  string qualifier, equal, type;
  unsigned int version;
    
  ISE_INPUT >> dummy >> qualifier >> equal >> version;
  ISE_INPUT >> qualifier >> equal >> type;
  //    
  //    TO DO !!! : integrity_check (version, type);
  //    
  //    
     
  //  unsigned int dimension, nb_vertices, nb_edges, nb_faces, nb_elements, nb_regions;
  //  ->  private member data 
  ISE_INPUT >> qualifier >> equal >> dimension;
  ISE_INPUT >> qualifier >> equal >> nb_vertices;
  ISE_INPUT >> qualifier >> equal >> nb_edges;
  ISE_INPUT >> qualifier >> equal >> nb_faces;
  ISE_INPUT >> qualifier >> equal >> nb_elements;
  ISE_INPUT >> qualifier >> equal >> nb_regions;
    
    
    
  // array  of  strings  for  "region" and  "material" labels
  string  regions[nb_regions], materials[nb_regions];
    
  ISE_INPUT >> qualifier >> equal >> dummy;
    
  for (int i = 0; i < nb_regions; i++)
  {
    ISE_INPUT >> regions[i];
    //cout << regions[i]<< endl;
  }
    
  ISE_INPUT >> dummy >> qualifier >> equal >> dummy;
    
  for (int i = 0; i < nb_regions; i++)
  {
    ISE_INPUT >> materials[i];
   // cout << materials[i]<< endl;
  }
    
  //    
  ////END OF INFO BLOCK VARIABLES

  //
  //
  //
  ////DATA section will  be   neglected !!!

  do
  {
    ISE_INPUT >> dummy;
  }
  while (dummy != "Vertices");
    
  do
  {
    ISE_INPUT >> dummy;
  }
  while (dummy != "{");
    
  ///  
    
    // *****************************
    // VERTICES   section 
    // *****************************
    //    
    //    
    
    vector<double> node_coord ;
    double coord_id;
	
    //	ISE_Vertex*  vertices[nb_vertices];   // not  possible!

    //vector<ISE_Vertex*> vertices; // member  data
    //  ISE_Vertex* vertex_point;
    vertices.clear();
     
	
    for (unsigned int i=0; i < nb_vertices ; i++)
    {
		
      for (unsigned int j = 0; j < dimension; j++)
      {
        ISE_INPUT >> coord_id ;
        node_coord.push_back( coord_id);
            
      }
		
      //vertices[i] = new  ISE_Vertex( node_coord , i);	
		
      vertex_point = new  ISE_Vertex( node_coord , i);	
		
      vertices.push_back(vertex_point);
		
		
      node_coord.clear();
		
    } 
    
    
    
    //
    // *****************************
    ////EDGES   section
    // *****************************
    //******************************
        //
    
        //  ISE_Edge* edge_point; // member data
        //	vector<ISE_Edge*> edges; //member data
        edges.clear();
    
    if (dimension > 1)
    {
      do
      {
        ISE_INPUT >> dummy;  //  da graffa  che  chiude Vertices
      }
      while (dummy != "{");
                     
      unsigned int id_1=0;
      unsigned int id_2=0;
                   
                  
    
      for (int i = 0; i < nb_edges; i++)
      {
        ISE_INPUT >> id_1; 
        ISE_INPUT >> id_2; //
        
        //	 edges[i] = new  ISE_Edge(vertices[id_1] , vertices[id_2] );	
       					 
        edge_point = new  ISE_Edge(vertices[id_1] , vertices[id_2] );	
        edges.push_back(edge_point);
       					 
                         
      }
    }
    
    
    
    //    
    ////Faces : ********************  TO  BE  IMPLEMENTED  **************************************
    //
    //
    //
    
    
  
    // *****************************************************************************
    //  ELEMENTS :  2D case
    // *****************************************************************************  
    
    
    
    //   ISE_Element* elements_list[nb_elements]; // NO
  
    // ISE_Element* elements_list_point;//->  member  data
    //	vector<ISE_Element*> elements_list_point;  //->  member  data
    elements_list.clear();
  
    int  vertex_0, vertex_1;
    //	unsigned int num_of_elem=0;  //  number of  2D  elements
    num_of_elem=0;
    
    vector<ISE_Edge*> edge_list;  // local  vector
    edge_list.clear();
    
    vector<bool> negative_edges; // local  vector
    negative_edges.clear();
	
    //vector<unsigned int> ISE_element_type_list; //  ->  member  data
	
    vector<unsigned int> :: iterator p;  //  iterator for  ISE_element_type_list vector 
    ISE_element_type_list.clear();
  	
    unsigned int  elem_type; ////
    int id;
	
    if (dimension == 2)
    {

      do
      {
        ISE_INPUT >> dummy;
      }
      while (dummy != "{");  //  begin  of Locations section
                     
      // ................................. neglected
      do
      {
        ISE_INPUT >> dummy;
      }
      while (dummy != "{");  //  begin  of Elements section
                     
                     
      for (unsigned int i = 0; i < nb_elements; i++)
      {
        ISE_INPUT >> elem_type ; //read element type;
                     	  
                     	                  	  
        switch(elem_type) 
        {
        case 0:
          //  0D element = Point 
     						 	
          break;
        case 1:
          //  1D   element = Segment (<vertex0,vertex1>)
     						
          ISE_INPUT >> vertex_0;
          //  cout  <<  vertex_0 << "   " ;
          ISE_INPUT >> vertex_1;
          // cout  <<  vertex_1 << "   " ;
          //	elements_list[i] = new  ISE_Element_1D(vertices[vertex_0] , vertices[vertex_1]);
	     						
          elements_list_point = new  ISE_Element_1D(vertices[vertex_0] , vertices[vertex_1]);
          elements_list_point->set_type(elem_type);
          elements_list_point->set_dimension(1);
	     						
          elements_list.push_back(elements_list_point);
	     						
          // 	(elements_list[i])->set_type(elem_type);
          //	(elements_list[i])->set_dimension(1);
   						   		
     						 
          break;
     						 	
        case 2: // Triangle
        case 3: //  Rectangle
          // ***** 2D elements *******
     						 
          // update list of  ISE element types:
          p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(),elem_type  );

          if (p == ISE_element_type_list.end())   //  not  found
          {
            ISE_element_type_list.push_back(elem_type);   						 		 
          }
       							 
          num_of_elem++ ; // increase number  of  2D  elements (for .xda file)
       							 
          // cout  <<  endl;
     						 
          // for  all  edges  of  the  element -> elem_type+1
          for (unsigned int j = 0; j<(elem_type+1); j++)
          {
            ISE_INPUT >> id;   // edge id
            // cout  <<  id << "   " ;
	     						
            if (id < 0)  //  negative  edge  id 
            {
              id = (-id-1);   //  ISE code  for inverted edge
              //  !!!!! INVERT THE  ORDER OF  THE  NODES OF  THIS  (NEGATIVE) EDGE :
              //  (1,2) ->   (2,1)
              negative_edges.push_back(true);  //  flag  for ISE_Element_2D::set_element_nodes()
              //cout << "TRUE ********** "<<endl;
              //cout << "negative_edges[j] before" << negative_edges[j] << endl;
            }
            else {negative_edges.push_back(false);}
            //  cout << endl;
            //cout << "negative_edges[j]" << negative_edges[j] << endl;
     							
            edge_list.push_back(edges[id]); 
       							
          } //   edges cycle	
	     					 
          //  cout << endl;
 						      
          //  creates new  Element_2D
                     	  
          //  elements_list[i] = new  ISE_Element_2D( edge_list);
          //   elements_list[i] = new  ISE_Element_2D( edge_list,negative_edges );
          elements_list_point = new  ISE_Element_2D( edge_list,negative_edges );
   						     
          elements_list_point->set_type(elem_type);
          elements_list_point->set_dimension(2);
   						     
   						
          //   (elements_list[i])->set_type(elem_type);
          //  (elements_list[i])->set_dimension(2);
   						 
          elements_list.push_back(elements_list_point);
   						   
          edge_list.clear();
          negative_edges.clear();
   						   
          break;
        default:
          cout <<  "Error : NOT implemented  yet"; 
        }		//  end  switch
                     	   
      } //  next  element
   						
    }
    
    
      
    if ( (dimension == 3) || (dimension == 1) )
    { 
      cout <<  "Error : NOT implemented  yet"; 
    }
	    
	    
    // *******************************
    //  END  ELEMENT  SECTION 
    // *******************************
        
	    
	    
	    
	    
    // ************************************
    //   PHYSICAL  REGIONS 
    // ***********************************
   
   
    int numb_region_elements;
    unsigned int phys_reg_id;

    //vector<unsigned int> regions_1D;  ->  member data
    //vector<unsigned int> regions_2D;
    regions_1D.clear();
    regions_2D.clear();

    vector<ISE_Element*>  region_elements_1D;
    vector<ISE_Element*>  region_elements_2D;


    for  (unsigned int i = 0; i < nb_regions; i++)
    {
      region_elements_1D.clear();
      region_elements_2D.clear();
	
      phys_reg_id = i+1;  //  phys  reg =  {1,...,nb_regions}
	
      do
      {
    	ISE_INPUT >> dummy;
      }
      while (dummy != "{");  
      //  begin  region data
    
      ////ISE_INPUT >> dummy >> qualifier >> equal >> dummy;
      do
      {
    	ISE_INPUT >> dummy;
      }
      while (dummy != "(");  
      //  find  numb  elements of  regions
      ISE_INPUT >> numb_region_elements;
      do
      {
    	ISE_INPUT >> dummy;
      }
      while (dummy != "{"); 
      // find  list of  elements id
      // cout << endl;
    
      for  (unsigned int j = 0; j < numb_region_elements; j++)
      {
    	ISE_INPUT >> id;  //  id  is  element id ,  that  is  its  position in the  array elements_list
    	
        //	cout <<  id <<  "   ";
    	  	
    	// ******************************************
    	// associate current physical  region (phys_reg_id)  to  element elements_list[id]
    	elements_list[id]->set_physical_region(phys_reg_id);
    	
    	// makes a vector of all elements nD in the  ISE phys. reg.
    	if  ( (elements_list[id]->get_dimension()) == 1)
        {
          region_elements_1D.push_back(elements_list[id]);
    	 	
        }
        else if ( (elements_list[id]->get_dimension()) == 2)
        {
          region_elements_2D.push_back(elements_list[id]);
    	 	
        }
    	
     	    	
      }
    
      //  makes vectors with 1D and  2D  ISE phys reg ID
      if  ( (elements_list[id]->get_dimension()) == 1)
      {
        regions_1D.push_back(phys_reg_id) ;
    		 
        //  cout << "regions_1D" <<  "   " << phys_reg_id << endl;
    		 
        // makes map <phys reg, elements>
        // map <unsigned int phys_reg_id, vector<ISE_Element*> region_elements_1D >
    		     		 
        map_1D_region_elements.insert(make_pair(phys_reg_id, region_elements_1D) );
    		  
      }
      else if ( (elements_list[id]->get_dimension()) == 2)
      {
        regions_2D.push_back(phys_reg_id) ;
     	 	
        // map <unsigned int phys_reg_id, vector<ISE_Element*>  region_elements_2D>
        // makes map <phys reg, elements>
        map_2D_region_elements.insert(make_pair(phys_reg_id, region_elements_2D) );
     	 	    	 	
      }
    
      // cout <<  endl;	
    
    
      phys_reg_id = 0;
    	
    }  //  next  phys. region
    
    
    //  END PHYSICAL REGION       
    //
    //***************************************
        ////END OF ISE  GRID  FILE
        //**********************************
            //   
            //  
 
            //  close  file  in  input
            ISE_INPUT.close(); 
 	 	
}


void ReadISEGrid::write_xda() 
{
	
  // ************************************************************
  //  xda (libmesh)  ELEM TYPES  [see enum_elem_type.h] 
  // ************************************************************

  //                  EDGE2=0,    // 0
  //                  EDGE3,      // 1
  //                  EDGE4,      // 2

  //                  TRI3,       // 3
  //                  TRI6,       // 4

  //                  QUAD4,      // 5
  //                  QUAD8,      // 6
  //                  QUAD9,      // 7

  //                  TET4,       // 8
  //                  TET10,      // 9

  //                  HEX8,       // 10
  //                  HEX20,      // 11
  //                  HEX27,      // 12

  //                  PRISM6,     // 13
  //                  PRISM15,    // 14
  //                  PRISM18,    // 15

  //                  PYRAMID5,   // 16
  //
  //
    
    
  // ************************************************************
  //    ISE   types 
  // ****************************************************************


  //  0   Point  (1 node)
  //  1  Segment  ( 2  nodes)
  //  2  Triangle (3 edges)
  //  3  Rectangle (4 edges)
  //  4  Polygon  (n  edges)
  //  5  Tetrahedron (4 faces)
  //  6  Pyramid     (5 faces)
  //  7  Prism       (5 faces)
  //  8  Brick       (6  faces)
  //  9  Tetrabrick  (7 faces)
  // 10 Polyhedron   (n  faces)

  //ISE_element_type_list

  // *****************************************************************
  // elem_type_conversion from ISE_type to  xda (libmesh) type :
  // ****************************************************************
  //
 
  vector<unsigned int> xda_elem_type;
   
  xda_elem_type.clear();
  
  unsigned int ISE_element_type_list_size = ISE_element_type_list.size();
  
  
  for (unsigned int i =0; i< ISE_element_type_list_size ;++i)
  { //for
    switch(ISE_element_type_list[i]) 
    {
      //   switch() {

    case 0:
      // 0D :  not  implemented in Libmesh

    case 1:
      // 1D:  SEGMENT  -> line EDGE2 = 0 (2 nodes)
      xda_elem_type.push_back(0);

      break;
	
    case 2:
      //   TRIANGLE 3nodes       // 3
      //   
      xda_elem_type.push_back(3);
      break;

    case 3:
      // Rectangle 4nodes,      // 5
      //   xda_elem_type = 5; 
      xda_elem_type.push_back(5);
      break;

    case 4:
      //    Polygon n nodes ,       
      //     TO  BE  IMPLEMENTED !!
      /// xda_elem_type.push_back(???);
      break;

    case 5:
      // Tetrahedron  4faces    // 
      //  elem_type =  TET4; 
      xda_elem_type.push_back(8);
      break;

    case 6:
      //  Pyramid 5faces      
      //    PYRAMID5,   // 16   
      xda_elem_type.push_back(16);
      break;

    case 7:
      //  Prism 5faces 
      //  elem_type = PRISM6,     // 13;
      xda_elem_type.push_back(13);
      break;

    case 8:
      // Brick 6 faces 
      //  elem_type = HEX8     // 10  ; 
      xda_elem_type.push_back(10);
     
      break;

    case 9:
      //  Tetrabrick ??????????
      //   elem_type = 4???; 
      // ???????  UNKNOWN
      break;


    case 10:
      //  Polyhedron n  faces  
      //  elem_type = ; 
      //  TO  BE  IMPLEMENTED !! 
      break;

    default:
      cout <<  "Error : wrong  value  of ISE_elem_type  "; 

    }
  }

 
 
  // *****************************************************
  //   SYNTAX OF  FILE  .XDA   
  // *****************************************************
  //    see ...................

  //DEAL 003:003
  //	 # Num. Elements  =  num of  elements
  //	 # Num. Nodes  =  num of  nodes
  //	 # Sum of Element Weights= sum of nodes of  all the  
  //             elements = Sigma(i=N_el)[ nodes(i)]
  //0	 # Num. Boundary Conds.  (= 0 )
  //65536	 # String Size (ignore)
  //	 # Num. Element Blocks. = "number  of  mesh  blocks"  !
  //	 # Element types in each block. = vector of  element types 
  //              (see  elem_type.h)
  // 	 # Num. of elements in each block = total number of elements 
  //            of a given type
  //Id String
  //Title String
  //  for  (all elements)
  // { list  of  nodes}
  //  for  (all  nodes) 
  // { list of  node coordinates}
  //



  //  Open   .xda  file 
  //
 
 
  string fname_xda;
  fname_xda = "in.xda";

  // Open the output file stream
  std::ofstream out_xda (fname_xda.c_str());
 
  assert (out_xda.good());
    
  //0	 # Num. Boundary Conds.  (= 0 )
  unsigned int num_bc =0;
 	
 	
  //	 # Num. Element Blocks. = "number  of  mesh  blocks"  !
  //     * A mesh block by definition contains
  //     * only a single type of element.
  //     *
  //     * The number of mesh blocks.
  unsigned int num_mesh_block =  xda_elem_type.size();
    


  // # Num. of elements in each block = total number of elements 
  //            of a given type
  //* The size of each element block is
      //     * the total number of a given type of
      //     * element in the mesh.
      //     *
      //     * @return The vector of block sizes

      //  num_elem_per_type =  to  be  calculated  !!
 
      // ************************************************************
      //  calculation of  list  of  elements( per  type block) 
      //  and  of  num_elem_per_type =  num  of  elements of  each of the  types in elem_type

      // num_of_elem = number of  elements  2D  //  member  data
      //  nb_elements  =  total  number  of  elements  // member  data
 
 
      //  vector<ISE_Element*> xda_list_elements; // member  data !!!!!!!!!!!
 
      xda_list_elements.clear();
  vector <unsigned int> num_elem_per_type;
  num_elem_per_type.clear();
   
  vector<unsigned int> elem_nodes_list;
  elem_nodes_list.clear();
   
  unsigned int  count = 0;
  unsigned int el_weight = 0;
  unsigned int sum_of_element_nodes = 0;
    
    
     
  for (unsigned int i =0; i< ISE_element_type_list_size;++i)
  { //for i

    for (unsigned int j =0; j< nb_elements ;++j)
    {//for j

      unsigned int elem_dim  = elements_list[j]->get_dimension();
      if (elem_dim == dimension) 
      {
        //elements_list[j]->get_type();
        if ( (elements_list[j]->get_type() )== ISE_element_type_list[i])
        {
          //    cout << "All_elements[j].type " << All_elements[j].type << endl;

          xda_list_elements.push_back(elements_list[j] );
          count++;  //  counter of  element of  current  type
         		 
          // ****************************************************
          //  !!!!!!!!! calculate number of  nodes in element !!!!!!!!!!!!!!
          // *****************************************************
          elem_nodes_list = (dynamic_cast<ISE_Element_2D*> (  elements_list[j] ) )-> get_nodes();
		
          unique_nodes(elem_nodes_list);  // nodes = list  of  nodes   belonging to  current  element, 
          // (without  repetitions !)
          //nodes_size= nodes_list.size();		
          sum_of_element_nodes = elem_nodes_list.size(); // 
         		 
          // ****************************************************
          // END     calculate number of  nodes in element !!!!!!!!!!!!!!
          // *****************************************************
          el_weight = el_weight + sum_of_element_nodes;

        }
      }
		
      elem_nodes_list.clear();

    }

    num_elem_per_type.push_back(count);
    count = 0;
    sum_of_element_nodes = 0;

  }
    
    
    
  unsigned int   xda_elem_type_size = xda_elem_type.size();
  unsigned int num_of_nodes = nb_vertices;



  //  ***********************************************************
  // writes  HEADER  of  .xda  file
  // ************************************************************
  //
  out_xda << "DEAL 003:003\n";
  out_xda << num_of_elem  <<   "	 # Num. Elements  \n";
  out_xda << num_of_nodes  <<   "	 # Num. Nodes \n";
  out_xda <<  el_weight <<   " 	 # Sum of Element Weights \n";
  out_xda <<  num_bc  <<   "	 # Num. Boundary Conds. \n";
  out_xda << "65536	 # String Size (ignore) \n";
  out_xda <<  num_mesh_block  <<   "	 # Num. Element Blocks.	\n";
  for (unsigned int i =0; i< xda_elem_type_size;i++)
  {

    out_xda << xda_elem_type[i]<< "   ";

  }
  out_xda << " 	 # Element types in each block. \n";
  for (unsigned int i =0; i<xda_elem_type_size ;++i)
  {

    out_xda <<  num_elem_per_type[i]<< "   ";

  }
  out_xda << "  	 # Num. of elements in each block. \n"; 

  out_xda << "Id String  \n";
  out_xda <<"Title   String  \n";
 
 
  // ***********************************
  // write   elements  section :  node list for  each  element
  // ***********************************

  vector<unsigned int> nodes_list;
  unsigned int nodes_size;

  for (unsigned int i =0; i< num_of_elem ; i++)
  {

	
    nodes_list = (dynamic_cast<ISE_Element_2D*> (xda_list_elements[i]) )-> get_nodes();
		
    unique_nodes(nodes_list);  // nodes = list  of  nodes   belonging to  current  element,  (without  repetitions !)
    nodes_size= nodes_list.size();
    for (unsigned int j =0; j< nodes_size; j++)

    {
      out_xda << nodes_list[j]  << "  " ;
    }
    out_xda << "\n";
		
		
  }
  
  
 

 
  // ***********************************
  // write  nodes  section 
  // ***********************************
 
  vector<double>  node_coordinates;

  for (unsigned int i=0; i < nb_vertices ; i++)
  {
    node_coordinates = vertices[i]->get_coord();
		
    for (unsigned int j=0; j < dimension ; j++)
    {
      out_xda << std::setprecision(10) << node_coordinates[j]   << "  ";
			 
    }
    if (dimension < 3 )
    {
      for (unsigned int k=0; k < (3-dimension) ; k++)
      {
        out_xda << std::setprecision(10) << 0.00 << "  ";
      }
				 	
    }
			
		
    out_xda << "\n";
		
  }
		
		
 
}

//  END  write  xda


void  ReadISEGrid::write_xta()
{
  // *************************************
  //  write xta  file
  // *************************************
	
	
  string fname;
  fname = "elem_data.xta";

  // Open the output file stream
  std::ofstream out (fname.c_str());
 
  assert (out.good());
  
  double  reg_id= 0.0;  
  unsigned int el_id =0;
  
  unsigned int elem_dim = 0;
 
  // WRITES  HEADER  OF  .XTA FILE
  
  out << " # Data description \n";
  out << "REAL	 # type of values     \n";
  out << "0	 # No. of nodes for which data is stored     \n";
  out << num_of_elem  << "	 # No. of elements for which data is stored    \n";
  
  //  num_of_elem  =  number  of  2D  elements !!
	
  for (unsigned int i =0; i<num_of_elem; i++) //  num_of_elem  =  number  of  2D  elements 
  {
	
    el_id = i;
	
    reg_id	  =  (double) (  (xda_list_elements[i])->get_physical_region() );
	 	 
    out << el_id   << "     # Foreign element id ";
	
    out << "\n";
    out << "1	 # vector length    ";
    out << "\n";
    out << std::setprecision(10) << reg_id   << "  	   # Values   ";
    out << "\n";
	    	 
  }
		
  // *************************************
  // END  write xta
  // *************************************
		
}






void  ReadISEGrid::unique_nodes(vector<unsigned int>& v1)
{
	
  vector<unsigned int> v_temp;
  unsigned int v_size = v1.size();
  vector<unsigned int>::iterator find_iter;
		
  v_temp.clear();
		
  for (unsigned int i=0; i < v_size ; i++)
  {
    //cout << "v1[i]" << v1[i]<<  endl; 
    find_iter = find( v_temp.begin(), v_temp.end(), v1[i]); 
			
    if ( find_iter  == v_temp.end() )
    {
      // not  found
      v_temp.push_back(v1[i]);
      //cout << "v1[i]" << v1[i]<<  endl; 
    }
  }
		
  v1.clear();
  v1 = v_temp;
	
	
}



//  makes  the  map  map <unsigned int tiber_BC_region, vector<unsigned int> region_nodes 
void  ReadISEGrid::set_BC_data()
{
  //	unsigned int elem_dim = 0;
  //	 vector<unsigned int> elem_nodes_list;
  //   elem_nodes_list.clear();
   
  // *****************************************************************************
   
  vector<ISE_Element*> BC_elements ;  // elements belonging to  a  single  phys reg
   
  unsigned int  phys_reg_id;
  unsigned int  BC_elements_size ;
  vector<unsigned int> BC_elem_nodes_list;
   
  vector<unsigned int> region_nodes;
   
  map <unsigned int , vector<ISE_Element*> >   :: iterator  p_reg_elem;
    
  unsigned int  regions_1D_size = regions_1D.size();
  unsigned int  regions_2D_size = regions_2D.size();
  
  unsigned int  tiber_BC_region =0;   // user  BC region ID  (should be consistent with input file !!!)
  
  if (dimension == 2)
  {
  
  
    for (unsigned int i =0; i<regions_1D_size  ;i++)
    {
	      	
      phys_reg_id = regions_1D[i];
	      	
      tiber_BC_region++;
	      	
	      	 
      p_reg_elem = map_1D_region_elements.find(phys_reg_id);
	      	 
	      	 
      if (p_reg_elem != map_1D_region_elements.end() )
      {
	      	 	
        BC_elements = p_reg_elem->second  ;
	      	 	
	      	 	
        //	cout << "phys_reg_id  , BC_elements_nodes   " << phys_reg_id << "  " << 
        //	((dynamic_cast<ISE_Element_1D*> (BC_elements[0] )  )->get_nodes() )[0]<< endl;
      }
      else
        cout  <<  "error in map_1D_region_elements"<< endl  ;
	    
      BC_elements_size = BC_elements.size();
      for (unsigned int j =0; j< BC_elements_size  ;j++)
      {
        BC_elem_nodes_list = (dynamic_cast<ISE_Element_1D*> (  BC_elements[j] ) )-> get_nodes();
        unique_nodes(BC_elem_nodes_list);  
	     		 	
        for (unsigned int k =0; k< BC_elem_nodes_list.size()  ;k++)
        {
          region_nodes.push_back(BC_elem_nodes_list[k]);
        }
	     		 	
        BC_elem_nodes_list.clear();
        //  next BC element		 	
      }
	    
      unique_nodes(region_nodes);
	   		 
      // makes map <unsigned int tiber_BC_region, vector<unsigned int> region_nodes )
	  	  
      // cout << "tiber_BC_region, region_nodes" << tiber_BC_region << "   "<<region_nodes[0] << endl;
	  	  
      map_BC_region_nodes.insert(make_pair(tiber_BC_region, region_nodes) );
	  
      region_nodes.clear();
      // next  region 
	      	
    }
	   
    //  for  all 1D_regions
    // BC_region++
    // from  map <phys_reg_id, 1D_elements>  ->  1D elements
    // elem_nodes_list = (dynamic_cast<ISE_Element_1D*> (  elements_list[j] ) )-> get_nodes();
    // region_nodes.push_back(elem_nodes_list )
    //  make  map <BC_region, region_nodes>
	  
	   
    //	for (unsigned int j =0; j< nb_elements ;++j)
    //      {//for j
    //
    //		elem_dim  = elements_list[j]->get_dimension();
    //		if (elem_dim == (dimension-1) )  //  (n-1) D elements ->  boundary  regions  !!!
    //		{
    //	
    //			elem_nodes_list = (dynamic_cast<ISE_Element_1D*> (  elements_list[j] ) )-> get_nodes();
    //		}
    //		
    //      }
	
  } //  end  dimension  = 2
  
  else if (dimension == 3)
  {
  	
    for (unsigned int i =0; i<regions_2D_size  ;i++)
    {
	      	
      phys_reg_id = regions_2D[i];
	      	
      tiber_BC_region++;
	      	
	      	 
      p_reg_elem = map_2D_region_elements.find(phys_reg_id);
	      	 
	      	 
      if (p_reg_elem != map_2D_region_elements.end() )
      {
	      	 	
        BC_elements = p_reg_elem->second  ;
      }
      else
        cout  <<  "error in map_1D_region_elements"<< endl  ;
	    
      BC_elements_size = BC_elements.size();
      for (unsigned int j =0; j< BC_elements_size  ;j++)
      {
        BC_elem_nodes_list = (dynamic_cast<ISE_Element_2D*> (  BC_elements[j] ) )-> get_nodes();
        unique_nodes(BC_elem_nodes_list);  
	     		 	
        for (unsigned int k =0; k< BC_elem_nodes_list.size()  ;k++)
        {
          region_nodes.push_back(BC_elem_nodes_list[k]);
        }
	     		 	
        BC_elem_nodes_list.clear();
        //  next BC element		 	
      }
	    
      unique_nodes(region_nodes);
	   		 
      // makes map <unsigned int tiber_BC_region, vector<unsigned int> region_nodes )
      map_BC_region_nodes.insert(make_pair(tiber_BC_region, region_nodes) );
	  
      region_nodes.clear();
      // next  region 
	      	
    }
  	
  	
  	
  	
  	
  }//  end  dimension  =3
  
  else cout  <<  "NOT IMPLEMENTED !" <<  endl;
		
	
	
}



//map <unsigned int , vector<unsigned int>  > 

void  ReadISEGrid::get_BC_data (map<unsigned int, vector<unsigned int> >& BoundCond_map )
{

  BoundCond_map = map_BC_region_nodes;
 
}



// write  mesh  and  meshdata from  .xda and  .xta  files

void   ReadISEGrid::read_mesh_and_data(Mesh& mesh, MeshData_elements&  mesh_data )

{

  //Mesh mesh (2);

  //MeshData_elements  mesh_data(mesh);
  //mesh_data.enable_compatibility_mode();

  // string  mesh_file_inp = "in.xda";

  string  mesh_file_data = "elem_data.xta";

  string  mesh_file_inp = "in.xda";

  //  string  mesh_file_inp = "in_pippo.xda";




  mesh.read (mesh_file_inp,&mesh_data  ); 
  
  mesh_data.read(mesh_file_data);


}




















