#include "ISE_Element_2D.h"




#include <vector>
#include <string>

#include <stdio.h>
#include <math.h>
#include <cassert>
using namespace  std;


ISE_Element_2D::~ISE_Element_2D()
{
}


// writes  nodes  belonging to element, both  as  ISE_Vertex*,  and  as  node_id
void ISE_Element_2D::set_element_nodes()
{
	
//cout << endl << "Setting 2D Element nodes: " << endl;
  ISE_Vertex* current_node;
	
  for (unsigned int i=0; i < element_edges.size(); i++)
  {
//cout <<  "  edge " <<  i << " : " << endl;
		
		
    if (negative_edges[i] == true)
    {
      //swap edge vertices !!!!!!!
//cout << "Negative Edge: Swapping vertices" << endl;


     for (unsigned int j=2; j > 0; j--)
     {
       current_node = (element_edges[i]->get_vertex(j) ) ;
//cout << " current_node =    " <<  current_node<< endl;
       element_nodes.push_back(current_node);
			
     }
    }

    else
    {
		
     for (unsigned int j=0; j < 2; j++)
     {
       current_node = (element_edges[i]->get_vertex(j+1) ) ;
//cout << " current_node =    " <<  current_node<< endl;
       element_nodes.push_back(current_node);
			
     }
    }
		
		
		
  }
	
	
	
	
	
}





void  ISE_Element_2D::unique_nodes_point()
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
		
	
	
}




void ISE_Element_2D::set_element_nodes_id()
{

  unsigned int current_node;
  unsigned int element_nodes_size = element_nodes.size();

  for (unsigned int i=0; i < element_nodes_size ; i++)
  {

    current_node = (element_nodes[i])-> get_node_id() ;
			
    element_nodes_id.push_back(current_node);


  }


}

