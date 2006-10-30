#include "ISE_Element_3D.h"

#include <vector>
#include <string>
#include <stdio.h>
#include <math.h>
#include <cassert>

using namespace std;

ISE_Element_3D::~ISE_Element_3D()
{
};


void ISE_Element_3D::set_element_nodes()
{
cout << endl << "Setting 3D Element nodes: " << endl;

  unsigned int e_size;
  ISE_Vertex* current_node;
  ISE_Edge* current_edge;
  bool edge_state = false;

	
  for (unsigned int i=0; i < element_faces.size(); i++)
  {

   e_size = element_faces[i]->get_face_edges_size();

cout <<  "  face " <<  i << " : " << endl;
		
		
    if (negative_faces[i] == true)
    {
      //swap face edges order and their orientation !!!!!!!;
cout << "Negative Face: Swapping edges order and orientation" << endl;

        

	for (int s = e_size-1; s >= 0; s--)
	{

cout << "edge " << s << ": " << endl;

	 current_edge = element_faces[i]->get_edge(s);
	 edge_state = element_faces[i]->get_edge_state(s);

	 if (!edge_state)
	 {

cout << "Negative edge, swapping vertices..." << endl;

	  	for (unsigned int j=2; j > 0; j--)
          	{
          	 current_node = (current_edge->get_vertex(j) ) ;
cout << " current_node =    " <<  current_node << endl;
          	 element_nodes.push_back(current_node);		
          	}
	 }

	 else
	 {	
     		for (unsigned int j=0; j < 2; j++)
     		{
      		 current_node = (current_edge->get_vertex(j+1) ) ;
cout << " current_node =    " <<  current_node<< endl;
       		 element_nodes.push_back(current_node);		
     		}
    	 }
	}


     
    }

    else
    {

	for (unsigned int e = 0; e < e_size; e++)
	{

cout << "edge " << e << ": " << endl;

	 current_edge = element_faces[i]->get_edge(e);
	 edge_state = element_faces[i]->get_edge_state(e);

	 if (edge_state)
	 {

cout << "Negative edge, swapping vertices..." << endl;

	  	for (unsigned int j=2; j > 0; j--)
          	{
          	 current_node = (current_edge->get_vertex(j) ) ;
cout << " current_node =    " <<  current_node << endl;
          	 element_nodes.push_back(current_node);		
          	}
	 }

	 else
	 {	
     		for (unsigned int j=0; j < 2; j++)
     		{
      		 current_node = (current_edge->get_vertex(j+1) ) ;
cout << " current_node =    " <<  current_node<< endl;
       		 element_nodes.push_back(current_node);		
     		}
    	 }
	}
		
    }
		
cout << endl;

  }	
	
	
	
	
};


void ISE_Element_3D::unique_nodes_point()
{
	
  vector<ISE_Vertex*>  v_temp;
	
  unsigned int 	element_nodes_size = element_nodes.size();
	
  vector<ISE_Vertex*> ::iterator find_iter;
		
  v_temp.clear();
		
  for (unsigned int i=0; i < element_nodes_size ; i++)
  {
			
    find_iter = find( v_temp.begin(), v_temp.end(), element_nodes[i] ); 
			
    if ( find_iter  == v_temp.end() )
    {
      // not  found
      v_temp.push_back(element_nodes[i]);
    }
  }
		
		
  element_nodes.clear();
  element_nodes = v_temp;
		
};




void ISE_Element_3D::set_element_nodes_id()
{

  unsigned int current_node;
  unsigned int element_nodes_size = element_nodes.size();

  for (unsigned int i=0; i < element_nodes_size ; i++)
  {

    current_node = (element_nodes[i])-> get_node_id() ;
			
    element_nodes_id.push_back(current_node);


  }


};

