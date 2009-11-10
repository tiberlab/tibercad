#include "ReadISEGrid.h"


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
  
  int elements_list_size = elements_list.size();
  int vertices_size = vertices.size();
  int edges_size = edges.size();
  int faces_size = faces.size();



  for (int i = 0; i < vertices_size; i++)
  {
    delete vertices[i];
  }

  for (int i = 0; i < edges_size; i++)
  {
    delete edges[i];
  }

  for (int i = 0; i < faces_size; i++)
  {
    delete faces[i];
  }

  //delete    pointers to  elements_list !!!
  for (int i = 0; i < elements_list_size; i++)
  {
    delete elements_list[i];
  }

}



void ReadISEGrid::integrity_check(int ver, string tp)
{
  if (ver == 1.0)
  {
    if (tp == "grid" || tp == "boundary")
    {return;}
    else
    {
      cerr << "Version or Type mismatch! Aborting program. " << endl;
      exit (1);
    }
  }

}

void ReadISEGrid::scan_grid_file()
{
	 
  //  local  pointers to ISE entities

  ISE_Vertex* vertex_point;
  ISE_Edge* edge_point;
  ISE_Face* face_point;
  ISE_Element* elements_list_point;


  //  
  //   lists of nD Elements with nD region    //   LOCAL   !!!!
  //
  vector<ISE_Element*>  region_elements_0D;
  vector<ISE_Element*>  region_elements_1D;
  vector<ISE_Element*>  region_elements_2D;
  vector<ISE_Element*>  region_elements_3D;



  ifstream ISE_INPUT(ISE_file_name);

  if ( !(ISE_INPUT.good()) )
  {
    //cout << "The path entered is incorrect or does not represent a valid input file, aborting program" 
    //  << endl << endl;
    exit(1);
  };

	
  //cout << "Beginning ISE_FILE Reading process: " << endl;	 


  //cout << "Opening file streaming... Done" << endl;

  string dummy;
     
      
  do
  {
    ISE_INPUT >> dummy;
  }
  while (dummy != "Info");
    


  //cout << "Reading Info Block Variables... ";
  ////     READ  INFO BLOCK VARIABLES
  //    
  string qualifier, equal, type;
  unsigned int version;
    
  ISE_INPUT >> dummy >> qualifier >> equal >> version;
  ISE_INPUT >> qualifier >> equal >> type;

  //cout << "Done" << endl;


  //cout << "Integrity check... ";

  integrity_check(version, type);   

  //cout << "Done" << endl;

     
  ISE_INPUT >> qualifier >> equal >> dimension;
  ISE_INPUT >> qualifier >> equal >> nb_vertices;
  ISE_INPUT >> qualifier >> equal >> nb_edges;
  ISE_INPUT >> qualifier >> equal >> nb_faces;
  ISE_INPUT >> qualifier >> equal >> nb_elements;
  ISE_INPUT >> qualifier >> equal >> nb_regions;
    

  //cout << endl << endl << "Variables: " << endl << endl;
  //cout << "Dimension: " << dimension << endl;
  //cout << "Vertices: " << nb_vertices << endl;
  //cout << "Edges: " << nb_edges << endl;
  //cout << "Faces: " << nb_faces << endl;
  //cout << "Elements: " << nb_elements << endl;
  //cout << "Regions: " << nb_regions << endl;

  //cout << endl;
    
    
  //cout << endl << "Reading Region and Material Information: " << endl;
  // array  of  strings  for  "region" and  "material" labels
  string  regions[nb_regions], materials[nb_regions];
    
  ISE_INPUT >> qualifier >> equal >> dummy;
    
  for (int i = 0; i < nb_regions; i++)
  {
    ISE_INPUT >> regions[i];
    //cout << "Region ID: " << regions[i] << endl;
  }
    
  ISE_INPUT >> dummy >> qualifier >> equal >> dummy;
    
  for (int i = 0; i < nb_regions; i++)
  {
    ISE_INPUT >> materials[i];
    //cout << "Material type: " << materials[i] << endl;
  }
    
  //    
  ////END OF INFO BLOCK VARIABLES
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
    
    vector<double> node_coord;
    double coord_id;
	
    vertices.clear();
    //cout << endl << endl << endl << "Vertices Section: " << endl << endl;
	
    for (unsigned int i=0; i < nb_vertices ; i++)
    {
		
      for (unsigned int j = 0; j < dimension; j++)
      {
        ISE_INPUT >> coord_id ;
        //cout << coord_id << "   ";
        node_coord.push_back( coord_id);
            
      }
		
      //cout << endl;
      vertex_point = new  ISE_Vertex( node_coord , i);	
		
      vertices.push_back(vertex_point);
			
      node_coord.clear();

		
    } 
    
    
    
    //
    // *****************************
    ////EDGES   section
    // *****************************
    //******************************

        edges.clear();
    
    if (dimension > 1)
    {
      do
      {
        ISE_INPUT >> dummy;  //  da graffa  che  chiude Vertices
      }
      while (dummy != "{");
                     
      unsigned int id_1;
      unsigned int id_2;
                   
      //cout << endl << endl << "Edges Section: " << endl;
    
      for (unsigned int i = 0; i < nb_edges; i++)
      {
        ISE_INPUT >> id_1; 
        ISE_INPUT >> id_2; 
        //cout << id_1 << "    " << id_2 << endl;
        				 
        edge_point = new  ISE_Edge(vertices[id_1] , vertices[id_2] );	
        edges.push_back(edge_point);

       					 
                         
      }
    }
    
    
    
    
    //
    // *****************************
    ////FACES   section
    // *****************************
    //******************************

       

        if (dimension == 3)
      {
        do
        {
          ISE_INPUT >> dummy;  //  da graffa  che  chiude Edges
        }
        while (dummy != "{");
                     
        unsigned int n_edg;
        int edge_id;
        vector<ISE_Edge*> face_edgs;
        vector<bool> neg_edges;
        face_edgs.clear();
        neg_edges.clear();

        //cout << endl << endl << "Faces Section: " << endl;

	for (unsigned int f = 0; f < nb_faces; f++)
	{
          ISE_INPUT >> n_edg;
          //cout << n_edg << "       ";

          for (unsigned int e = 0; e < n_edg; e++)
          {
            ISE_INPUT >> edge_id;
            //cout << edge_id << "  ";

            if (edge_id < 0)
            {
              edge_id = -edge_id-1;
              neg_edges.push_back(true);
            }
            else
            {
              neg_edges.push_back(false);
            }

            face_edgs.push_back(edges[edge_id]);
          }

          //cout << endl;
          face_point = new  ISE_Face(face_edgs, neg_edges);	
		
          faces.push_back(face_point);

          face_edgs.clear();
          neg_edges.clear();
	} 
      }
                   
                  
    
  
    // *****************************************************************************
    //                                    ELEMENTS
    // *****************************************************************************  
    
    elements_list.clear();
  
    int  vertex_0, vertex_1;
    
    vector<ISE_Edge*> edge_list;  // local  vector
    edge_list.clear();
    
    vector<ISE_Face*> face_list;  // local  vector
    face_list.clear();
    

    vector<bool> negative_edges; // local  vector
    negative_edges.clear();

    vector<bool> negative_faces; // local  vector
    negative_faces.clear();
	
	
    vector<unsigned int> :: iterator p;  //  iterator for  ISE_element_type_list vector 
    ISE_element_type_list.clear();
  	
    unsigned int  elem_type, elem_edges, elem_faces; //local variables
    int id; //local variable   
	

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
                     
    //cout << endl << endl << endl << "Element Section: " << endl << endl;
                     
    for (unsigned int i = 0; i < nb_elements; i++)
    {
      ISE_INPUT >> elem_type ; //read element type;

      //cout << "Element type: " << elem_type << "  " << endl;                     	  
                     	                  	  
      switch(elem_type) 
      {
      case 0:
	 
        //  // NO: never in XDA list !  !!!!update list of  ISE element types:
        //           p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(), elem_type);

        //           if (p == ISE_element_type_list.end())   //  not  found
        //           {
        //             ISE_element_type_list.push_back(elem_type);  
						 		 
        //           }


        //  0D element = Point 
        //cout << "0D element = Point: " << endl;
        ISE_INPUT >> vertex_0;
        //cout << "vertex: " << vertex_0 << endl;

        elements_list_point = new  ISE_Element_0D(vertices[vertex_0]);
        elements_list_point->set_type(elem_type);
        elements_list_point->set_dimension(0);
	     						
        elements_list.push_back(elements_list_point);


        //cout << "0D element inserted into list..." << endl << endl;

        break;

      case 1:
     			
        if (dimension == 1)  //  if  dim = 2D or  3D, 1D elements don't go in xda file !!
        {		 
          // update list of  ISE element types:
          p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(), elem_type);

          if (p == ISE_element_type_list.end())   //  not  found
          {
            ISE_element_type_list.push_back(elem_type);  
						 		 
          }
        }

        //  1D element = Segment (<vertex0,vertex1>)
        //cout << "1D element = Segment (<vertex0,vertex1>): " << endl;
     						
        ISE_INPUT >> vertex_0;
        //cout  <<  "vertex_0: " << vertex_0 << "   " ;
        ISE_INPUT >> vertex_1;
        //cout  <<  "vertex_1: " << vertex_1 << "   " ;

        elements_list_point = new  ISE_Element_1D(vertices[vertex_0] , vertices[vertex_1]);
        elements_list_point->set_type(elem_type);
        elements_list_point->set_dimension(1);
	     						
        elements_list.push_back(elements_list_point);


        //cout << endl << "1D element inserted into list..." << endl << endl;

     						 
        break;
     						 	
      case 2: // Triangle
      case 3: //  Rectangle
        // ***** 2D elements *******
     						 
        if (dimension == 2)  //  if  dim = 3D, 2D elements don't go in xda file, only BC !!
        {
          // update list of  ISE element types:
          p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(), elem_type);

          if (p == ISE_element_type_list.end())   //  not  found
          {
            ISE_element_type_list.push_back(elem_type);  
						 		 
          }

        }


       							 
          
        //cout << "2D element = Triangle/Rectangle: " << endl;
     						
       							 
        //cout  << "edge id: " << endl;
     						 
        // for  all  edges  of  the  element -> elem_type+1
        for (unsigned int j = 0; j<(elem_type+1); j++)
        {
          ISE_INPUT >> id;   // edge id
          //cout  << id << "   " ;
	     						
          if (id < 0)  //  negative  edge  id 
          {
            id = (-id-1);   //  ISE code  for inverted edge
            //  !!!!! INVERT THE  ORDER OF  THE  NODES OF  THIS  (NEGATIVE) EDGE :
            //  (1,2) ->   (2,1)
            negative_edges.push_back(true);  //  flag  for ISE_Element_2D::set_element_nodes()
            //cout << "NEGATIVE EDGE ********** ";
          }
          else {negative_edges.push_back(false);}
          //cout << endl;
					
          edge_list.push_back(edges[id]); 
       							
        } //   edges cycle	
	     					 
        //cout << endl;
 						      
        //  creates new  Element_2D

        elements_list_point = new  ISE_Element_2D( edge_list, negative_edges );
   						     
        elements_list_point->set_type(elem_type);
        elements_list_point->set_dimension(2);
   						 
        elements_list.push_back(elements_list_point);


        //cout << "2D element type inserted into list..." << endl << endl; 
   						   
        edge_list.clear();
        negative_edges.clear();
 
        break;





      case 4: //  Polygon
        // ***** 2D element *******
     					
        //cout << "2D element = Polygon: " << endl;
        if (dimension == 2)  //  if  dim = 3D, 2D elements don't go in xda file, only BC !!
        {					 
          // update list of  ISE element types:
          p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(), elem_type);

          if (p == ISE_element_type_list.end())   //  not  found
          {
            ISE_element_type_list.push_back(elem_type);  
						 		 
          }
        }

       							 
       				
        ISE_INPUT >> elem_edges;
        //cout << "Number of Edges:  " << elem_edges << endl;
			 
        //cout  << "edge id: " << endl;
     						 
        // for  all  edges  of  the  element
        for (unsigned int j = 0; j<(elem_edges); j++)
        {
          ISE_INPUT >> id;   // edge id
          //cout  << id << "   " ;
	     						
          if (id < 0)  //  negative  edge  id 
          {
            id = (-id-1);   //  ISE code  for inverted edge
            //  !!!!! INVERT THE  ORDER OF  THE  NODES OF  THIS  (NEGATIVE) EDGE :
            //  (1,2) ->   (2,1)
            negative_edges.push_back(true);  //  flag  for ISE_Element_2D::set_element_nodes()
            //cout << "NEGATIVE EDGE ********** ";
          }
          else {negative_edges.push_back(false);}
          //cout << endl;
					
          edge_list.push_back(edges[id]); 
       							
        } //   edges cycle	
	     					 
        //cout << endl;
 						      
        //  creates new  Element_2D

        elements_list_point = new  ISE_Element_2D( edge_list, negative_edges );
   						     
        elements_list_point->set_type(elem_type);
        elements_list_point->set_dimension(2);
   						 
        elements_list.push_back(elements_list_point); 


        //cout << "2D element type inserted into list..." << endl << endl;
   						   
        edge_list.clear();
        negative_edges.clear();
          

        break;




      case 5: //  Tetrahedron
        // ***** 3D element *******
     						
        //cout << "3D element = Tetrahedron: " << endl;
     						 
        // update list of  ISE element types:
        p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(), elem_type);

        if (p == ISE_element_type_list.end())   //  not  found
        {
          ISE_element_type_list.push_back(elem_type);  	 
        }
       							 
       	
        //cout  << "face id: " << endl;
     						 
        // for  all  faces  of  the  element
        for (unsigned int j = 0; j<4; j++)
        {
          ISE_INPUT >> id;   // face id
          //cout  << id << "   " ;
	     						
          if (id < 0)  //  negative  face  id 
          {
            id = (-id-1);   //  ISE code  for inverted face
            //  !!!!! INVERT THE  ORDER (AND ORIENTATION) OF  THE  EDGES OF  THIS  (NEGATIVE) FACE :
            //  (1,2,3,4) ->   (-4,-3,-2,-1)
            negative_faces.push_back(true);  //  flag  for ISE_Element_3D::set_element_nodes()
            //     cout << "NEGATIVE FACE ********** ";
          }
          else {negative_faces.push_back(false);}
          //cout << endl;

          //    negative_faces.push_back(false);  //  TEST  !!!
					
          face_list.push_back(faces[id]); 
       							
        } //   faces cycle	
	     					 
        //cout << endl;
 						      
        //  creates new  Element_3D

        //    elements_list_point = new  ISE_Element_3D( face_list, negative_faces );
        elements_list_point = new  ISE_Element_3D( face_list, negative_faces,elem_type );

   						     
        //       elements_list_point->set_type(elem_type);
        elements_list_point->set_dimension(3);
   						 
        elements_list.push_back(elements_list_point);

 
        //cout << "3D Tetrahedron element type inserted into list..." << endl << endl;
						 	
   						   
        face_list.clear();
        negative_faces.clear();
          
	  

        break;




      case 6:
      case 7: //  Pyramid/Prism
        // ***** 3D element *******
     		
        //cout << "3D element = Pyramid/Prism: " << endl;
     										 
        // update list of  ISE element types:
        p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(), elem_type);

        if (p == ISE_element_type_list.end())   //  not  found
        {
          ISE_element_type_list.push_back(elem_type);    
        }
       							 
       	
        //cout  << "face id: " << endl;
     						 
        // for  all  faces  of  the  element
        for (unsigned int j = 0; j<5; j++)
        {
          ISE_INPUT >> id;   // face id
          //cout  << id << "   " ;
	     						
          if (id < 0)  //  negative  face  id 
          {
            id = (-id-1);   //  ISE code  for inverted face
            //  !!!!! INVERT THE  ORDER (AND ORIENTATION) OF  THE  EDGES OF  THIS  (NEGATIVE) FACE :
            //  (1,2,3,4,5) ->   (-5,-4,-3,-2,-1)
            negative_faces.push_back(true);  //  flag  for ISE_Element_3D::set_element_nodes()
            //cout << "NEGATIVE FACE ********** ";
          }
          else {negative_faces.push_back(false);}
          //cout << endl;
					
          face_list.push_back(faces[id]); 
       							
        } //   faces cycle	
	     					 
        //cout << endl;
 						      
        //  creates new  Element_3D

        elements_list_point = new  ISE_Element_3D( face_list, negative_faces,elem_type );
   						     
        //       elements_list_point->set_type(elem_type);
        elements_list_point->set_dimension(3);
   						 
        elements_list.push_back(elements_list_point);


        //cout << "3D Pyramid/Prism element type inserted into list..." << endl << endl;
						 		
   						   
        face_list.clear();
        negative_faces.clear();



          
        break;



      case 8: //  Brick
        // ***** 3D element *******
     						 
        //cout << "3D element = Brick: " << endl;
     						
        // update list of  ISE element types:
        p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(), elem_type);

        if (p == ISE_element_type_list.end())   //  not  found
        {
          ISE_element_type_list.push_back(elem_type);   
        }
       							 
       	
        //cout  << "face id: " << endl;
     						 
        // for  all  faces  of  the  element
        for (unsigned int j = 0; j<6; j++)
        {
          ISE_INPUT >> id;   // face id
          //cout  << id << "   " ;
	     						
          if (id < 0)  //  negative  face  id 
          {
            id = (-id-1);   //  ISE code  for inverted face
            //  !!!!! INVERT THE  ORDER (AND ORIENTATION) OF  THE  EDGES OF  THIS  (NEGATIVE) FACE :
            //  (1,2,3,4,5,6) ->   (-6,-5,-4,-3,-2,-1)
            negative_faces.push_back(true);  //  flag  for ISE_Element_3D::set_element_nodes()
            //cout << "NEGATIVE FACE ********** ";
          }
          else {negative_faces.push_back(false);}
          //cout << endl;
					
          face_list.push_back(faces[id]); 
       							
        } //   faces cycle	
	     					 
        //cout << endl;
 						      
        //  creates new  Element_3D

        elements_list_point = new  ISE_Element_3D( face_list, negative_faces,elem_type );
   						     
        //       elements_list_point->set_type(elem_type);
        elements_list_point->set_dimension(3);
   						 
        elements_list.push_back(elements_list_point);

   						   
        //cout << "3D Brick element type inserted into list..." << endl << endl;
						 		 
        face_list.clear();
        negative_faces.clear();

          
        break;




      case 9: //  Tetrabrick
        // ***** 3D element *******
     					
        //cout << "3D element = Tetrabrick: " << endl;
     							 
        // update list of  ISE element types:
        p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(), elem_type);

        if (p == ISE_element_type_list.end())   //  not  found
        {
          ISE_element_type_list.push_back(elem_type);   
        }
       							 
       	
        //cout  << "face id: " << endl;
     						 
        // for  all  faces  of  the  element
        for (unsigned int j = 0; j<7; j++)
        {
          ISE_INPUT >> id;   // face id
          //cout  << id << "   " ;
	     						
          if (id < 0)  //  negative  face  id 
          {
            id = (-id-1);   //  ISE code  for inverted face
            //  !!!!! INVERT THE  ORDER (AND ORIENTATION) OF  THE  EDGES OF  THIS  (NEGATIVE) FACE :
            //  (1,2,3,4,5,6,7) ->   (-7,-6,-5,-4,-3,-2,-1)
            negative_faces.push_back(true);  //  flag  for ISE_Element_3D::set_element_nodes()
            //cout << "NEGATIVE FACE ********** ";
          }
          else {negative_faces.push_back(false);}
          //cout << endl;
					
          face_list.push_back(faces[id]); 
       							
        } //   faces cycle	
	     					 
        //cout << endl;
 						      
        //  creates new  Element_3D

        elements_list_point = new  ISE_Element_3D( face_list, negative_faces,elem_type );
   						     
        //       elements_list_point->set_type(elem_type);
        elements_list_point->set_dimension(3);
   						 
        elements_list.push_back(elements_list_point);

   						   
        //cout << "3D Tetrabrick element type inserted into list..." << endl << endl;
						 		 
        face_list.clear();
        negative_faces.clear();

	  

        break;




      case 10: //  Polyhedron
        // ***** 3D element *******
     					
        //cout << "3D element = Polyhedron: " << endl;
     							 
        // update list of  ISE element types:
        p = find(ISE_element_type_list.begin(), ISE_element_type_list.end(), elem_type);

        if (p == ISE_element_type_list.end())   //  not  found
        {
          ISE_element_type_list.push_back(elem_type);   	 
        }
       							 
       	
        ISE_INPUT >> elem_faces;
        //cout  << "Number of Element Faces: " << elem_faces << endl;
        //cout  << "face id: " << endl;
     						 
        // for  all  faces  of  the  element
        for (unsigned int j = 0; j<(elem_faces); j++)
        {
          ISE_INPUT >> id;   // face id
          //cout  << id << "   " ;
	     						
          if (id < 0)  //  negative  face  id 
          {
            id = (-id-1);   //  ISE code  for inverted face
            //  !!!!! INVERT THE  ORDER (AND ORIENTATION) OF  THE  EDGES OF  THIS  (NEGATIVE) FACE :
            //  (1,2,3,...,n) ->   (-n,...,-3,-2,-1)
            negative_faces.push_back(true);  //  flag  for ISE_Element_3D::set_element_nodes()
            //cout << "NEGATIVE FACE ********** ";
          }
          else {negative_faces.push_back(false);}
          //cout << endl;
					
          face_list.push_back(faces[id]); 
       							
        } //   faces cycle	
	     					 
        //cout << endl;
 						      
        //  creates new  Element_3D

        elements_list_point = new  ISE_Element_3D( face_list, negative_faces,elem_type  );
   						     
        //       elements_list_point->set_type(elem_type);
        elements_list_point->set_dimension(3);
   						 
        elements_list.push_back(elements_list_point);


        //cout << "3D Polyhedron element type inserted into list..." << endl << endl;
						 	
   						   
        face_list.clear();
        negative_faces.clear();


        break;

      default:
        //cout << "Error : NOT implemented yet... is it possible? :)" << endl<<endl; 
        break;
      }		//  end  switch
                     	   
    } //  next  element
   						
    

	    
    // *******************************
    //     END  ELEMENT  SECTION 
    // *******************************
        



  
    // *****************************************************************************
    //                           PHYSICAL REGIONS
    // *****************************************************************************  
    
    //cout << endl << "Begin of Physical Region Section: " << endl << endl;

    int numb_region_elements;
    unsigned int phys_reg_id,tiber_phys_reg_id, tiber_BC_reg_id  ;
    bool  tib_reg_incremented,tib_BC_reg_incremented  ;
    string region_name;

    phys_reg_id = 0;

    tiber_phys_reg_id = 0;
    tiber_BC_reg_id= 0;

    regions_0D.clear();
    regions_1D.clear();
    regions_2D.clear();
    regions_3D.clear();

    region_names_ID_map.clear();
    BC_region_names_ID_map.clear();


    //    increment_tiber_region = false;
    tib_reg_incremented = false;
    tib_BC_reg_incremented= false;


    ISE_INPUT >> dummy; //  read closing  bracket } of  element section

    for  (unsigned int i = 0; i < nb_regions; i++)
    { // for nb_regions

      // regions_0D.clear();
      //       regions_1D.clear();
      //       regions_2D.clear();
      //       regions_3D.clear();

      region_elements_0D.clear();
      region_elements_1D.clear();
      region_elements_2D.clear();
      region_elements_3D.clear();
	
      phys_reg_id = i+1;  //  phys  reg =  {1,...,nb_regions}



      //cout << "Region[" << phys_reg_id << "]: " << endl;
	
      //********************************************************
          // read  region  name

 
          ISE_INPUT >> dummy;  //  read "Region"  keyword  

      ISE_INPUT >> region_name;   // region name in the  form: ("region_1")
      //   cerr<< "region_name *********" << region_name << endl;

      string::size_type loc = region_name.find_first_of( '(', 0 );
      if( loc != string::npos )
      {
        //     cout << "Found  at " << loc << endl;
      }
      else
        //     cout << "Didn't find  " << endl; 
        throw InitFailedException("ERROR in  reading ISE grid file (physical regions section)."); 

      region_name.erase(loc, loc+2);
      //     cout  <<  "  NEW region_name:   "<< region_name <<  endl << endl;

      string::size_type loc2 = region_name.find( '"', 0 );
      if( loc2 != string::npos )
      {
        //     cout << "Found # at " << loc2 << endl;
      }
      else
        //    cout << "Didn't find # " << endl; 
        throw InitFailedException("ERROR in  reading ISE grid file (physical regions section)."); 

      region_name.erase(loc2);
      //   cout  <<  "  NEW region_name:   "<< region_name <<  endl << endl;

      // exit(1);

      //********************************************************



          //   do
          //       {
          //     	ISE_INPUT >> dummy;
          //       }
          //       while (dummy != "{");  
          //       //  begin  region data



    
        do
        {
          ISE_INPUT >> dummy;  //  read  "material =..." and "Elements" keyword
        }
        while (dummy != "(");  
        //  find  numb  elements of  regions
        ISE_INPUT >> numb_region_elements;

        //cout << "Region Elements Number: " << numb_region_elements << endl;
        do
        {
          ISE_INPUT >> dummy;
        }
        while (dummy != "{"); 



        //list of  elements id

        //cout << "Elements id: " << endl;
    
        for  (unsigned int j = 0; j < numb_region_elements; j++)
        {
          ISE_INPUT >> id;  //  id  is  element id ,  that  is  its  position in the  array elements_list
    	
          //cout <<  id <<  "   ";
    	  	
          // ******************************************
          // associate current physical  region (phys_reg_id)  to  element elements_list[id]
          //	elements_list[id]->set_physical_region(phys_reg_id); 
          //!!see below, only  for elem_dim = sim_dim
    	
          // makes a vector of all elements nD in the  ISE phys. reg.


          switch ( (elements_list[id]->get_dimension()) )
          {

          case 0:  //  element_dimension = 0!
            {
              region_elements_0D.push_back(elements_list[id]);

              //  regions_0D.push_back(phys_reg_id) ;

              //             // makes map <phys reg, elements>
    		     		 
              //             map_0D_region_elements.insert(make_pair(phys_reg_id, region_elements_0D) );


    		 
              //cout << "0D Element inserted into appropriate map" << endl;
    		 
    	 	
            };

            break;


          case 1: //  element_dimension = 1!
            {
              region_elements_1D.push_back(elements_list[id]);



              if ( (dimension ==2) ||  (dimension ==3) )
              {

                if  (tib_BC_reg_incremented ==  false)
                {
                  tib_BC_reg_incremented =  true;
                  tiber_BC_reg_id++;
                }
              }


              //  regions_1D.push_back(phys_reg_id) ;
              //             cerr << "regions_1D = " << phys_reg_id <<  endl;

              //             // makes map <phys reg, elements>
    		     		 
              //             map_1D_region_elements.insert(make_pair(phys_reg_id, region_elements_1D) );

              //             //cout << "1D Element inserted into appropriate map" << endl;
    		 
    	 	
            };

            break;


          case 2://  element_dimension = 2!
            {

              // element_dimension = simulation dimension ->  to tiber_phys_region !!
              if  (dimension == 2)
              {// 

                if  (tib_reg_incremented ==  false)
                {
                  tib_reg_incremented =  true;
                  tiber_phys_reg_id++;
                }
              
                elements_list[id]->set_physical_region(tiber_phys_reg_id);
              }

              if  (dimension ==3) 
              {

                if  (tib_BC_reg_incremented ==  false)
                {
                  tib_BC_reg_incremented =  true;
                  tiber_BC_reg_id++;
                }
              }




              region_elements_2D.push_back(elements_list[id]);

 




              //  regions_2D.push_back(phys_reg_id) ;

              //             cerr << "regions_2D = " << phys_reg_id <<  endl;

              //             // makes map <phys reg, elements>
    		     		 
              //             map_2D_region_elements.insert(make_pair(phys_reg_id, region_elements_2D) );


              //             //cout << "2D Element inserted into appropriate map" << endl;
    		 
    	 	
            };

            break;


          case 3://  element_dimension = 3!
            {
        

              // element_dimension = simulation dimension ->  to tiber_phys_region !!
              if  (dimension == 3)
              {
                if  (tib_reg_incremented ==  false)
                {
                  tib_reg_incremented =  true;
                  tiber_phys_reg_id++;
                }
              
                elements_list[id]->set_physical_region(tiber_phys_reg_id);

              }

              region_elements_3D.push_back(elements_list[id]);


              //  regions_3D.push_back(phys_reg_id) ;

              //             // makes map <phys reg, elements>
    		     		 
              //             map_3D_region_elements.insert(make_pair(phys_reg_id, region_elements_3D) );


              //             //cout << "3D Element inserted into appropriate map" << endl;
    		 
    	 	
            };

            break;

          default:
            {


              cerr << " ERROR: Could not insert element " << id 
                   << " into appropriate region/element vector" << endl;
              throw InitFailedException(""); 

            };
          } //switch end
     	    	
        	

        } // end of  present region elements

     
        //  if  tib_reg_incremented =  true ->   it's a  physical region -> put it  in the  map
        // put region name in  map with  tiber_phys_reg_id
        if   (tib_reg_incremented == true)
        {
          //   cerr<< " put in map physical region :  " <<tiber_phys_reg_id<< "  "<<  region_name<< endl;
          region_names_ID_map.insert(make_pair(tiber_phys_reg_id, region_name)  ); 
        }
        //  else put  in  BC_region_names_ID_map
        //     else BC_region_names_ID_map.insert(make_pair(tiber_phys_reg_id, region_name)  );
        if   (tib_BC_reg_incremented == true)
        {
          //    cerr<< " put in map BC region :  " << tiber_BC_reg_id << "  "<<region_name<< endl;
          BC_region_names_ID_map.insert(make_pair(tiber_BC_reg_id, region_name)  ); 
        }  


        tib_reg_incremented =  false; // reset flag  for  next region
        tib_BC_reg_incremented =  false; // reset flag  for  next region

        //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //  here   regions_nD.push_back  and  map_nD_region_elements.insert

        switch ( (elements_list[id]->get_dimension()) )
        {

        case 0:
          {
            

            regions_0D.push_back(phys_reg_id) ;

            // makes map <phys reg, elements>
    		     		 
            map_0D_region_elements.insert(make_pair(phys_reg_id, region_elements_0D) );


    		 
            //cout << "0D Element inserted into appropriate map" << endl;
    		 
    	 	
          };

          break;


        case 1:
          {
           

            regions_1D.push_back(phys_reg_id) ;
            //    cerr << "regions_1D = " << phys_reg_id <<  endl;

            // makes map <phys reg, elements>
    		     		 
            map_1D_region_elements.insert(make_pair(phys_reg_id, region_elements_1D) );

            //cout << "1D Element inserted into appropriate map" << endl;
    		 
    	 	
          };

          break;


        case 2:
          {
           

            regions_2D.push_back(phys_reg_id) ;

            //   cerr << "regions_2D = " << phys_reg_id <<  endl;

            // makes map <phys reg, elements>
    		     		 
            map_2D_region_elements.insert(make_pair(phys_reg_id, region_elements_2D) );


            //cout << "2D Element inserted into appropriate map" << endl;
    		 
    	 	
          };

          break;


        case 3:
          {
            

            regions_3D.push_back(phys_reg_id) ;

            // makes map <phys reg, elements>
    		     		 
            map_3D_region_elements.insert(make_pair(phys_reg_id, region_elements_3D) );


            //cout << "3D Element inserted into appropriate map" << endl;
    		 
    	 	
          };

          break;

        default:
          {
            cerr << " ERROR: Could not insert element " << id 
                 << " into appropriate region/element vector" << endl;
            exit (1);
          };
        } //switch end





    
        //cout << endl;
    
        phys_reg_id = 0;

 
        ISE_INPUT >> dummy; //  read first closing  bracket } of  physical region section
        ISE_INPUT >> dummy; //  read second closing  bracket } of  physical region section


    }  //  next  phys. region
    
    
    //   cerr  <<  "regions_2D.size   "  <<  regions_2D.size() <<  endl ;

    //  END OF PHYSICAL REGIONS       

    //***************************************
        //        END OF ISE  GRID  FILE
        //***************************************
            //  close  input file
            ISE_INPUT.close(); 


    // *******************deallocation  local pointers
    // delete vertex_point;
    // delete edge_point;
    // delete edge_point;
    //  delete elements_list_point;




    //cout << "File reading... Done" << endl << "Closing file stream... Done" << endl << endl;
 	 	
} // end  scan_grid_file

















void ReadISEGrid::write_xda() 
{

  //cout << "Writing xda file: In progress..." << endl << endl;
  //cout << "Converting ISE element types into Libmesh ones..." << endl;
	
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
 
  vector<unsigned int> xda_elem_type; // Local Vector
   
  xda_elem_type.clear();
  
  const unsigned int ISE_element_type_list_size = ISE_element_type_list.size();
  

  //cout << ISE_element_type_list_size << " element types to convert: " << endl;

  
  for (unsigned int i =0; i< ISE_element_type_list_size ;i++)
  {

    //cout << "ISE Element type " << ISE_element_type_list[i] << " -> ";
    switch(ISE_element_type_list[i]) 
    {

    case 0:
      // 0D :  not  implemented in Libmesh
      //cout << "Element not implemented in Libmesh!! Nothing to be done (Aborting Program)." << endl;
      exit(1);

      break;

    case 1:
      // 1D:  SEGMENT  -> line EDGE2 = 0 (2 nodes)
      //cout << "Libmesh Element type 0 " << endl;

      xda_elem_type.push_back(0);

      break;
	
    case 2:
      //   TRIANGLE 3nodes       // 3

      //cout << "Libmesh Element type 3 " << endl;

      xda_elem_type.push_back(3);
      break;

    case 3:
      // Rectangle 4nodes,      // 5
      //cout << "Libmesh Element type 5 " << endl;

      xda_elem_type.push_back(5);
      break;

    case 4:
      //    Polygon n nodes
      //cout << "Element not implemented in Libmesh!! Nothing to be done (Aborting Program)." << endl;
      exit(1);


      break;

    case 5:
      // Tetrahedron  4faces   
      //cout << "Libmesh Element type 8 " << endl;

      xda_elem_type.push_back(8);
      break;

    case 6:
      //  Pyramid 5faces      
      //cout << "Libmesh Element type 16 " << endl;

      xda_elem_type.push_back(16);
      break;

    case 7:
      //  Prism 5faces 
      //cout << "Libmesh Element type 13 " << endl;

      xda_elem_type.push_back(13);
      break;

    case 8:
      // Brick 6 faces 
      //cout << "Libmesh Element type 10 " << endl; 

      xda_elem_type.push_back(10);
     
      break;

    case 9:
      //  Tetrabrick 
      //cout << "Element not implemented in Libmesh!! Nothing to be done (Aborting Program)." << endl;
      exit(1);


      break;


    case 10:
      //  Polyhedron n  faces 
      //cout << "Element not implemented in Libmesh!! Nothing to be done (Aborting Program)." << endl;
      exit(1);


      break;

    default:
      //cout <<  "Error : wrong  value  of ISE_elem_type  "; 
      break;

    }
  }
 

  // *****************************************************
  //   SYNTAX OF  FILE  .XDA   
  // *****************************************************
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
 
  //  fname_xda = "a.xda";
  fname_xda = "in.xda";

  // Open the output file stream
  std::ofstream out_xda (fname_xda.c_str());
 
  if ( !( out_xda.good() ) )
  {
    cerr << "The name of the output file is incorrect, aborting program." << endl;
    exit(1);
  };
    
  //0	 # Num. Boundary Conds.  (= 0 )
  unsigned int num_bc =0; //local variable
 	
 	
  //	 # Num. Element Blocks. = "number  of  mesh  blocks"  !
  //     * A mesh block by definition contains
  //     * only a single type of element.
  unsigned int num_mesh_block =  xda_elem_type.size();
    


  // # Num. of elements in each block = total number of elements 
  //            of a given type
  //* The size of each element block is
      //     * the total number of a given type of
      //     * element in the mesh.
 
      // ************************************************************
      //  calculation of  list  of  elements( per  type block) 
      //  and  of  num_elem_per_type =  num  of  elements of  each type in elem_type

 
      vector <unsigned int> num_elem_per_type; //Local Vector
  num_elem_per_type.clear();
   
  vector<unsigned int> elem_nodes_list;  //Local Vector
  elem_nodes_list.clear();
   
  unsigned int  count = 0; //Local Variable
  unsigned int el_weight = 0; //Local Variable
  unsigned int sum_of_element_nodes = 0; //Local Variable
    
    
     
  xda_list_elements.clear();
  unsigned int el_dimension = 0;

  for (unsigned int i = 0; i < ISE_element_type_list_size; i++)
  { //for i

    for (unsigned int j = 0; j < nb_elements; j++)
    {//for j


      if ( (elements_list[j]->get_type() )== ISE_element_type_list[i])

      {


        xda_list_elements.push_back( elements_list[j] );
        count++;  //  counter of  element of  current  type
         		 
        el_dimension = elements_list[j]->get_dimension();

	
        elem_nodes_list.clear();


        switch (el_dimension)
	{

	case 1:
          elem_nodes_list = (dynamic_cast<ISE_Element_1D*> (  elements_list[j] ) )-> get_nodes_id();
          break;

	case 2:
          elem_nodes_list = (dynamic_cast<ISE_Element_2D*> (  elements_list[j] ) )-> get_nodes_id();
          break;

	case 3:
          elem_nodes_list = (dynamic_cast<ISE_Element_3D*> (  elements_list[j] ) )-> get_nodes_id();
          break;

	default:
	  cerr << "ERROR: Wrong element dimension." << endl;
	  cerr << "Cannot execute switch (el_dimension) to get elem_nodes_list value." << endl;
	  cerr << "Aborting program." << endl;
	  exit(1);
		
	} // switch end


        unique_nodes(elem_nodes_list);  // nodes = list  of  nodes   belonging to  current  element, 
        // (without  repetitions !)
	
        sum_of_element_nodes = elem_nodes_list.size();
         		 
        el_weight += sum_of_element_nodes;

      }
     
		

    }

    num_elem_per_type.push_back(count);
    count = 0;
    sum_of_element_nodes = 0;

  }
    
    
    
  unsigned int xda_elem_type_size = xda_elem_type.size();  //Local Variable
  unsigned int num_of_nodes = nb_vertices;  //Local Variable
  unsigned int el_tot = 0; //Local Variable

  for (unsigned int t = 0; t < num_elem_per_type.size(); t++)
  {
    el_tot += num_elem_per_type[t];
  };



  //  ***********************************************************
  // writes  HEADER  of  .xda  file
  // ************************************************************
  //
  out_xda << "DEAL 003:003\n";
  out_xda << el_tot  <<   "	 # Num. Elements  \n";
  out_xda << num_of_nodes  <<   "	 # Num. Nodes \n";
  out_xda <<  el_weight <<   " 	 # Sum of Element Weights \n";
  out_xda <<  num_bc  <<   "	 # Num. Boundary Conds. \n";
  out_xda << "65536	 # String Size (ignore) \n";
  out_xda <<  num_mesh_block  <<   "	 # Num. Element Blocks.	\n";


  for (unsigned int i = 0; i < xda_elem_type_size; i++)
  {
    out_xda << xda_elem_type[i] << "   ";
  };

  out_xda << " 	 # Element types in each block. \n";

  for (unsigned int i = 0; i < xda_elem_type_size ; i++)
  {
    out_xda <<  num_elem_per_type[i] << "   ";
  }

  out_xda << "  	 # Num. of elements in each block. \n"; 

  out_xda << "Id String  \n";
  out_xda <<"Title   String  \n";
 
 
  // ***********************************
  // write   elements  section :  node list for  each  element
  // ***********************************

  vector<unsigned int> nodes_list; //Local Vector
  unsigned int nodes_size, xda_el_dim; // Local Variables

  for (unsigned int i = 0; i < xda_list_elements.size(); i++)
  {

    xda_el_dim = xda_list_elements[i]->get_dimension();

    switch (xda_el_dim)
    {

    case 1:
      nodes_list = (dynamic_cast<ISE_Element_1D*> (  xda_list_elements[i] ) )-> get_nodes_id();
      break;

    case 2:
      nodes_list = (dynamic_cast<ISE_Element_2D*> (  xda_list_elements[i] ) )-> get_nodes_id();
      break;

    case 3:
      nodes_list = (dynamic_cast<ISE_Element_3D*> (  xda_list_elements[i] ) )-> get_nodes_id();
      break;

    default:
      cerr << "ERROR: Wrong xda element dimension." << endl;
      cerr << "Cannot execute switch (xda_el_dim) to get nodes_list value." << endl;
      cerr << "Aborting program." << endl;
      exit(1);
		
    } // switch end	
    unique_nodes(nodes_list);  // nodes = list  of  nodes   belonging to  current  element,  (without  repetitions !)
    nodes_size = nodes_list.size();

    for (unsigned int j = 0; j < nodes_size; j++)

    {
      out_xda << nodes_list[j]  << "  " ;
    }
    out_xda << "\n";
		
		
  }
  
  
 

 
  // ***********************************
  // write  nodes  section 
  // ***********************************
 
  vector<double>  node_coordinates; //Local Vector

  for (unsigned int i = 0; i < nb_vertices ; i++)
  {
    node_coordinates = vertices[i]->get_coord();
		
    for (unsigned int j = 0; j < dimension ; j++)
    {
      out_xda << std::setprecision(10) << node_coordinates[j]   << "  ";	 
    }

    if (dimension < 3 )
    {
      for (unsigned int k=0; k < (3-dimension) ; k++)
      {
        out_xda << std::setprecision(10) << 0.00 << "  ";
      }
				 	
    };
			
		
    out_xda << "\n";
		
  }
		
  out_xda.close();		
 

  //cout << endl << "Writing *.xda file... Done" << endl << endl;

}

//  END  write  xda





















void ReadISEGrid::write_xta() 
{

  //cout << "Writing xta file: In progress" << endl;

  // //cout << endl << "Insert the full name (with no path will be in the same directory) of the *.xta output file: " << endl;


  // *************************************
  //  write xta  file
  // *************************************
	
	


  //  fname_xta = "a.xta";
  fname_xta ="elem_data.xta";

  // Open the output file stream
  std::ofstream out (fname_xta.c_str());
 
  if ( !( out.good() ) )
  {
    cerr << "The name of the output file is incorrect, aborting program." << endl;
    exit(1);
  };


  double  reg_id= 0.0;  
  unsigned int el_id =0;
  
  unsigned int elem_dim = 0;
 
  // WRITES  HEADER  OF  .XTA FILE
  
  out << " # Data description \n";
  out << "REAL	 # type of values     \n";
  out << "0	 # No. of nodes for which data is stored     \n";
  out << xda_list_elements.size()  << "	 # No. of elements for which data is stored    \n";
  


  for (unsigned int i =0; i< xda_list_elements.size(); i++)
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
		

  out.close();


  //cout << endl << "Writing *.xta file... Done" << endl << endl;
}























//! Utility    to  eliminate  repetitions  in  node  list

void  ReadISEGrid::unique_nodes(vector<unsigned int>& v1)
{
	
  vector<unsigned int> v_temp;
  unsigned int v_size = v1.size();
  vector<unsigned int>::iterator find_iter;
		
  v_temp.clear();
		
  for (unsigned int i=0; i < v_size ; i++)
  {

    find_iter = find( v_temp.begin(), v_temp.end(), v1[i]); 
			
    if ( find_iter  == v_temp.end() )
    {
      // not  found
      v_temp.push_back(v1[i]);

    }
  }
		
  v1.clear();
  v1 = v_temp;
	
	
}




//!  Writes map which associates  nodes belonging to Boundary regions
// [ (dimension-1) elements] and relative regions

//  makes  the  map  map <unsigned int tiber_BC_region, vector<unsigned int> region_nodes 
void  ReadISEGrid::set_BC_data()
{
   
  vector<ISE_Element*> BC_elements ;  // elements belonging to  a  single  phys reg
   
  unsigned int  phys_reg_id;
  unsigned int  BC_elements_size ;
  vector<unsigned int> BC_elem_nodes_list;
   
  vector<unsigned int> region_nodes;
   
  map <unsigned int , vector<ISE_Element*> >   :: iterator  p_reg_elem;
    
  unsigned int  regions_1D_size = regions_1D.size();
  unsigned int  regions_2D_size = regions_2D.size();
  
  unsigned int  tiber_BC_region = 0;   // user  BC region ID  (should be consistent with input file !!!)
  
  //  cerr  <<  " *********** regions_1D_size  =  " <<  regions_1D_size <<  endl ;
  //   cerr  <<  " *********** regions_2D_size  =  " <<  regions_2D_size <<  endl ;


  if (dimension == 2)
  {
  
  
    for (unsigned int i =0; i<regions_1D_size  ;i++)
    {
	      	
      phys_reg_id = regions_1D[i];
	      	
      tiber_BC_region++;
      //   cout  <<  " *********** tiber_BC_region  =  " <<  tiber_BC_region <<  endl ;
	      	 
      p_reg_elem = map_1D_region_elements.find(phys_reg_id);
	      	 
	      	 
      if (p_reg_elem != map_1D_region_elements.end() )
      {
	      	 	
        BC_elements = p_reg_elem->second  ;
	      	 	
	      	 	
      }
      else
      {
        //cout  <<  "error in map_1D_region_elements"<< endl;
	exit(1);
      };
	    
      BC_elements_size = BC_elements.size();
      for (unsigned int j =0; j< BC_elements_size  ;j++)
      {
        BC_elem_nodes_list = (dynamic_cast<ISE_Element_1D*> (  BC_elements[j] ) -> get_nodes_id());
        unique_nodes(BC_elem_nodes_list);  
	     		 	
        for (unsigned int k =0; k< BC_elem_nodes_list.size()  ;k++)
        {
          region_nodes.push_back(BC_elem_nodes_list[k]);
        }
	     		 	
        BC_elem_nodes_list.clear();
        		 	
      }  //  next BC element
	    
      unique_nodes(region_nodes);
	   		 
      map_BC_region_nodes.insert(make_pair(tiber_BC_region, region_nodes) );
	  
      region_nodes.clear();
      
	      	
    }  // next  region 
	   
	
  } //  end  dimension  = 2
  
  else if (dimension == 3)
  {
  	
    for (unsigned int i =0; i<regions_2D_size  ;i++)
    {
	      	
      phys_reg_id = regions_2D[i];
	      	
      tiber_BC_region++;

      //    cout  <<  " *********** tiber_BC_region  =  " <<  tiber_BC_region <<  endl ;
	      	
	      	 
      p_reg_elem = map_2D_region_elements.find(phys_reg_id);
	      	 
	      	 
      if (p_reg_elem != map_2D_region_elements.end() )
      { 	
        BC_elements = p_reg_elem->second  ;
      }
      else
      {
        //cout  <<  "error in map_2D_region_elements"<< endl;
	exit (1);
      };
	    
      BC_elements_size = BC_elements.size();
      for (unsigned int j =0; j< BC_elements_size; j++)
      {
        BC_elem_nodes_list = (dynamic_cast<ISE_Element_2D*> (  BC_elements[j] )-> get_nodes_id());
        unique_nodes(BC_elem_nodes_list);  
	     		 	
        for (unsigned int k =0; k< BC_elem_nodes_list.size()  ;k++)
        {
          region_nodes.push_back(BC_elem_nodes_list[k]);
        }
	     		 	
        BC_elem_nodes_list.clear();
        		 	
      } //  next BC element
	    
      unique_nodes(region_nodes);
	   		 
      map_BC_region_nodes.insert(make_pair(tiber_BC_region, region_nodes) );
	  
      region_nodes.clear();
       
	      	
    } // next  region
  	
  	
  	
  	
  	
  }//  end  dimension  = 3
  
  else
    cout  <<  "NOT IMPLEMENTED !" <<  endl;
		
	
}


//  returns BC_region_map
void  ReadISEGrid::get_BC_data (map<unsigned int, vector<unsigned int> >& BoundCond_map )
//map<unsigned int, vector<unsigned int> >&  ReadISEGrid::get_BC_data()
{

  BoundCond_map = map_BC_region_nodes;
  //return map_BC_region_nodes;
 
}


// write  mesh  and  meshdata from  .xda and  .xta  files

void   ReadISEGrid::read_mesh_and_data(Mesh& mesh, MeshData_elements&  mesh_data )

{

  //Mesh mesh (2);

  //MeshData_elements  mesh_data(mesh);
  //mesh_data.enable_compatibility_mode();

  // string  mesh_file_inp = "in.xda";

  //  string  mesh_file_data = "elem_data.xta";

  //   string  mesh_file_inp = "in.xda";

  //  string  mesh_file_inp = "in_pippo.xda";


  //fname_xta ="elem_data.xta";
  // fname_xda = "in.xda";

  //  mesh.read (mesh_file_inp,&mesh_data  ); 

  mesh.read (fname_xda,&mesh_data  ); 
  

  //  mesh_data.read(mesh_file_data);

  mesh_data.read(fname_xta);



}





void   ReadISEGrid::get_region_names_map (map<ID, string >& region_names_map )

{

  region_names_map = region_names_ID_map;


}




void   ReadISEGrid::get_BC_region_names_map (map<ID, string >& BC_region_names_map )

{

  BC_region_names_map = BC_region_names_ID_map;

}

