/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file ISE_Element_3D.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */

#include "ISE_Element_3D.h"

#include <vector>
#include <string>
#include <cstdio>
#include <cmath>

#include <cassert>

using namespace std;

ISE_Element_3D::~ISE_Element_3D()
{
}


void ISE_Element_3D::set_element_nodes()
{
  //cout << endl << "Setting 3D Element nodes: " << endl;

  unsigned int e_size;
  ISE_Vertex* current_node;
  ISE_Edge* current_edge;
  bool edge_state = false;

	
  for (unsigned int i=0; i < element_faces.size(); i++)
  {

    e_size = element_faces[i]->get_face_edges_size();

    //cout <<  "  face " <<  i << " : " << endl;
		
		
    if (negative_faces[i] == true)
    {
      //swap face edges order and their orientation !!!!!!!;
      //cout << "Negative Face: Swapping edges order and orientation" << endl;

        

      for (int s = e_size-1; s >= 0; s--)
      {

        //cout << "edge " << s << ": " << endl;

        current_edge = element_faces[i]->get_edge(s);
        edge_state = element_faces[i]->get_edge_state(s);

        if (!edge_state)
        {

          //cout << "Negative edge, swapping vertices..." << endl;

          for (unsigned int j=2; j > 0; j--)
          {
            current_node = (current_edge->get_vertex(j) ) ;
            //cout << " current_node =    " <<  current_node << endl;
            element_nodes.push_back(current_node);		
          }
        }

        else
        {	
          for (unsigned int j=0; j < 2; j++)
          {
            current_node = (current_edge->get_vertex(j+1) ) ;
            //cout << " current_node =    " <<  current_node<< endl;
            element_nodes.push_back(current_node);		
          }
        }
      }


     
    }

    else
    {

      for (unsigned int e = 0; e < e_size; e++)
      {

        //cout << "edge " << e << ": " << endl;

        current_edge = element_faces[i]->get_edge(e);
        edge_state = element_faces[i]->get_edge_state(e);

        if (edge_state)
        {

          //cout << "Negative edge, swapping vertices..." << endl;

          for (unsigned int j=2; j > 0; j--)
          {
            current_node = (current_edge->get_vertex(j) ) ;
            //cout << " current_node =    " <<  current_node << endl;
            element_nodes.push_back(current_node);		
          }
        }

        else
        {	
          for (unsigned int j=0; j < 2; j++)
          {
            current_node = (current_edge->get_vertex(j+1) ) ;
            //cout << " current_node =    " <<  current_node<< endl;
            element_nodes.push_back(current_node);		
          }
        }
      }
		
    }
		
    //cout << endl;

  }	
	
	
  check_orientation_3D();  //  necessary  for  consistency of nodes order,  normals, etc.	
	
}


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
		
}




void ISE_Element_3D::set_element_nodes_id()
{

  unsigned int current_node;
  unsigned int element_nodes_size = element_nodes.size();

  for (unsigned int i=0; i < element_nodes_size ; i++)
  {

    current_node = (element_nodes[i])-> get_node_id() ;
			
    element_nodes_id.push_back(current_node);


  }


}




void  ISE_Element_3D::check_orientation_3D()
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
 	
  vector<double>  coord;
  coord.clear();


  coord = element_nodes[0]->get_coord();

  x0 = coord[0];
  y0 = coord[1]; 
  z0 = coord[2];
  coord.clear();

  coord = element_nodes[1]->get_coord();
  x1 = coord[0];
  // cout << "x1 " << x1 << endl;
  y1 = coord[1];
  // cout << "y1 " << y1 << endl;
  z1 = coord[2];

  coord.clear();
  
  coord = element_nodes[2]->get_coord();
  x2 = coord[0];
  
  y2 = coord[1];
  z2 = coord[2];

  coord.clear();

  coord = element_nodes[3]->get_coord();
  x3 = coord[0];
 
  y3 = coord[1];
  z3 = coord[2];
  coord.clear();





  // unsigned int  elem_type;

  //  HEX8   ISE 8  Brick

  if (elem_type == 8)
  {

    coord = element_nodes[4]->get_coord();
    x4 = coord[0];
 
    y4 = coord[1];
    z4 = coord[2];
    coord.clear();


 
    a11 = x1-x0;
    a21 =   y1-y0;
    a31 = z1-z0;

    a12 = x3-x0;
    a22 = y3 - y0;
    a32 = z3- z0;

    a13 =x4 - x0;
    a23 = y4 - y0;
    a33 = z4 - z0;

  }

  ////  TET4  ISE 5  Tetrahedron 
  // PRISM6   ISE 7  Prism 

  else  if ( (elem_type == 5) || (elem_type == 7)) //    TET4 , PRISM6

 

  {
    a11 = x1-x0;
    a21 =   y1-y0;
    a31 = z1-z0;

    a12 = x2 - x0;
    a22 = y2-y0;
    a32 = z2-z0;

    a13 = x3-x0 ;
    a23 = y3-y0;
    a33= z3-z0;
  }

  // det 

  det =  a11 * (a22*a33 - a23*a32 ) - a12 * (a21*a33 - a23*a31)+ 
    a13 * (a21*a32 - a22*a31);

  //|a_1 a_2 a_3; b_1 b_2 b_3; c_1 c_2 c_3| ==
  // a_1b_2c_3-a_1b_3c_2-a_2b_1c_3+a_2b_3c_1+a_3b_1c_2-a_3b_2c_1

  assert(abs(det)> 1e-12);

  if  (det > 0.0)
  {swap = false;}
  else 

  {
    swap = true;

    // cout <<  "  SWAP !!!  " ; //<< endl ;
  }

  if (swap)     
  {// if  swap = true  (det <= 0 )   then  reorder nodes of  element

    //   cerr << " ************elem_type " << elem_type <<  endl;
   
    switch(elem_type) {
      //   switch(ISE_elem_type) 
    case 5: //  TET4

      // swap node_id_list[1] and  node_id_list[3]
      temp = 0;
      temp = element_nodes[1];

      element_nodes[1] =element_nodes[3];
      element_nodes[3]= temp;
       

      break;

    case 7: // ISE  PRISM6
    {

      ISE_Vertex*  temp_list[3];

      // swap   (0 1 2) -> (3 4 5 )


      for (unsigned int i=0; i<3; i++)
      {

        temp_list[i] = element_nodes[i];

      }

      for (unsigned int i=0; i<3; i++)

      {
        element_nodes[i] =element_nodes[i+3] ;
      }

      for (unsigned int i=0; i<3; i++)

      {
        element_nodes[i+3] = temp_list[i]  ;

      }

      break;
    }

    case 8 : // ISE HEX8
    {

      ISE_Vertex*  temp_list[4];
      // swap   (0 1 2 3) -> ( 4 5  6  7  )

      for (unsigned int i=0; i<4; i++)
      {

        temp_list[i] = element_nodes[i];

      }

      for (unsigned int i=0; i<4; i++)

      {
        element_nodes[i] =element_nodes[i+4] ;
      }

      for (unsigned int i=0; i<4; i++)

      {
        element_nodes[i+4] = temp_list[i]  ;

      }


      break;
    }
	   
      // 

    default:
      cerr <<  "Error : ISE_elem_type NOT YET IMPLEMENTED "; 

    }  //  end  switch

  } // endif (swap)     

}
