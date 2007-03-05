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


// //writes  nodes  belonging to element, both  as  ISE_Vertex*,  and  as  node_id
// void ISE_Element_2D::set_element_nodes()
// {
	
// //cout << endl << "Setting 2D Element nodes: " << endl;
//   ISE_Vertex* current_node;
	
//   for (unsigned int i=0; i < element_edges.size(); i++)
//   {
// //cout <<  "  edge " <<  i << " : " << endl;
		
		
//     if (negative_edges[i] == true)
//     {
//       //swap edge vertices !!!!!!!
// //cout << "Negative Edge: Swapping vertices" << endl;


//      for (unsigned int j=2; j > 0; j--)
//      {
//        current_node = (element_edges[i]->get_vertex(j) ) ;
// //cout << " current_node =    " <<  current_node<< endl;
//        element_nodes.push_back(current_node);
			
//      }
//     }

//     else
//     {
		
//      for (unsigned int j=0; j < 2; j++)
//      {
//        current_node = (element_edges[i]->get_vertex(j+1) ) ;
// //cout << " current_node =    " <<  current_node<< endl;
//        element_nodes.push_back(current_node);
			
//      }
//     }
		
		
		
//   }
	


// // temp  !!!!
// //*********************	
// //  check_orientation_2D();	
	
	
// }


// **********************************-----------------------------


// OLD  VERSION !!!
//writes  nodes  belonging to element, both  as  ISE_Vertex*,  and  as  node_id
void ISE_Element_2D::set_element_nodes()
{
	
  unsigned int  last,  last2;
  // cout <<  " element_edges.size() =   " << element_edges.size()<<  endl;
  //unsigned int current_node =0 ;
	
  ISE_Vertex* current_node;
	
  for (unsigned int i=0; i < element_edges.size(); i++)
  {
    //cout <<  "  edge " <<  i << endl;
		
		
		
    for (unsigned int j=0; j < 2; j++)
    {
      current_node = (element_edges[i]->get_vertex(j+1) ) ;
      //	cout << " current_node =    " <<  current_node<< endl;
      element_nodes.push_back(current_node);
			
      //element_nodes_id.push_back(current_node->get_node_id() );
			
    }
		
    if (negative_edges[i] == true)
    {
      //swap last  2  nodes !!!!!!!;
			
      //	cout << "element_nodes[last] =  " << element_nodes_id[last2]<< endl;
      //	cout << "element_nodes[last-1] =  "<<  element_nodes_id[last2-1]<<  endl;
			
      last= (element_nodes.size() ) - 1;
      swap(element_nodes[last], element_nodes[last-1]);
			
      //last2= (element_nodes_id.size() ) - 1;
      //swap(element_nodes_id[last2], element_nodes_id[last2-1]);
			
			
      //	cout << "element_nodes[last] =  "<<  element_nodes_id[last2]<< endl;
      //	cout << "element_nodes[last-1] =  "<<  element_nodes_id[last2-1]<<  endl;
			
			
    }
			
		
		
  }
	
	
  //  unique_nodes (element_nodes)
	
  // unique_nodes_point(element_nodes)
	
	
  check_orientation_2D();

	
	
	
	
	
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



// ******************************************************
//   temp  :   check  of  nodes  orientation in  2D  element

void  ISE_Element_2D::check_orientation_2D()
{
 	
  ISE_Vertex*  temp;
 	
  bool swap;
  swap = false;
	 	 
  double det  ;    //    determinant
  double x0,x1,z0, y0,y1,z1,  x2,y2, z2, x3, y3, z3 ,x4,y4,z4, 
    a11,a21, a12, a22, a31, a32, a13, a23, a33;
	
  // 	 unsigned int  size_node_list ;
	
  // 	size_node_list =  element_nodes.size();
	  	
  unique_nodes_point();
 	

  // OK  FOR  FIRST  ORDER ELEMENTS
  //   SECOND ORDER ...............TO  DO  !!?
  
  
  
 
  //  if  (dim == 2)
  //  ok for  TRI4  AND  QUAD4 

  // {

  vector<double>  coord;
  coord.clear();

  coord = element_nodes[0]->get_coord();

  //inline vector<double>
  //ISE_Vertex::get_coord()

  x0 = coord[0];
  // cout << "x0 " << x0 << endl;
  y0 = coord[1];
  // cout << "y0 " << y0 << endl;
  coord.clear();
  
  coord = element_nodes[1]->get_coord();
  x1 = coord[0];
  // cout << "x1 " << x1 << endl;
  y1 = coord[1];
  // cout << "y1 " << y1 << endl;
  coord.clear();
  
  coord = element_nodes[2]->get_coord();
  x2 = coord[0];
  //  cout << "x2 " << x2 << endl;
  y2 = coord[1];
  // cout << "y2 " << y2 << endl;
  


  //    x0 = node_coord[( temp_nodes_list[0]) ][0];
  //    x1 = node_coord[( temp_nodes_list[1]) ][0];
  //
  //
  //
  //    y0 = node_coord[( temp_nodes_list[0]) ][1];
  //    y1 = node_coord[( temp_nodes_list[1]) ][1];
  //
  //    x2 = node_coord[( temp_nodes_list[2]) ][0];
  //    y2  = node_coord[( temp_nodes_list[2]) ][1];

  a11 = x1-x0;
  a21 =   y1-y0;

  a12 = x2 - x0;
  a22 = y2-y0;


  det = a11 * a22  - (a21*a12) ;

  //cout << fabs(det) ;

  assert(fabs(det)> 1e-12);

  if  (det > 0.0)
  {swap = false;}
  else swap = true;

  if (swap)   // if  swap = true  (det <= 0 )   then  reorder nodes of  element
  {

      
    // swap node_id_list[0] and  node_id_list[2]
      
    temp = element_nodes[0];
    element_nodes[0]= element_nodes[2];
    element_nodes[2] = temp;
      
    //      temp = node_id_list[0];
    //      node_id_list[0] = node_id_list[2];
    //      node_id_list[2] = temp;
      
  //   cout << " change orientation  !!!  " << endl;
   //   cout << " det =     " << det << endl;
    //	  count_ok++;

  }
  else  // swap = false  (det >  0 )  
  {
      
    // cout << " det =     " << det << endl;
    //  cout  << " SWAP NODES  !!!  " << endl;
    //  count_swap++;

    // swap node_id_list[0] and  node_id_list[2]
    //      temp = node_id_list[0];
    //      node_id_list[0] = node_id_list[2];
    //      node_id_list[2] = temp;


  }

	
}
 
