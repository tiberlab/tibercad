/*=============================================================================
    Copyright (c) 2002-2003 Joel de Guzman
    http://spirit.sourceforge.net/

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

///////////////////////////////////////////////////////////////////////////////
#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/dynamic.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

#include <vector>
#include <string>


#include "Read_MSH.h"


///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace boost::spirit;

///////////////////////////////////////////////////////////////////////////////
//
//  Class  for  parser of  GMSH  .MSH file (v. 1) :   get  nodes of  BC regions (1D  physical regions)
//   -------------------------------------------
//  MUST   include  SPIRIT  libraries !
// -----------------------------------------------
//
///////////////////////////////////////////////////////////////////////////////


//  constructor:  scan  .msh  file,  parse it, get info on BC,  writes .xda and  .xta  file for  Libmesh
//  needs physical region ID,  bound. conditions region  ID (2D or 1D), dim. of  simulation 
Read_MSH::Read_MSH(string filename, vector<unsigned int>& phys_reg_ID, vector<unsigned int>& BC_reg_ID, unsigned int sim_dim, Mesh& mesh, MeshData_elements&  mesh_data)


  //Read_MSH::Read_MSH(string filename, vector<unsigned int>& phys_reg_ID, vector<unsigned int>& BC_reg_ID, int sim_dim  )

{

  dim =  sim_dim; 

  initialize_vectors(); 


  scan_input(filename);

  //  get_BC_info();

  get_BC_info(BC_reg_ID); //  ->  after   write_xda(); ?????

  //  get_physical_elem(phys_reg_ID);  //  ->  after   write_xda();  !!!!!


  //  get_elem_nodes();

  get_nodes_coord();
  get_elem_nodes();



  write_xda();

  get_physical_elem(phys_reg_ID);

  write_xta();

  read_mesh_and_data(mesh,mesh_data );


}


Read_MSH::~Read_MSH()
{
}


void Read_MSH::initialize_vectors()
{

  elem_values.clear();
  
}



//  parser of  ELEMENT  section of  .msh  file  (version 1)
void Read_MSH::parse_elem_section(ifstream& in_stream )
{
     
  unsigned int n_elem;
  string  str;
  vector<unsigned int> v;
  v.clear();
     
  //  rule  for element  line  (list  of  integers uint)
  rule<>list_of_numbers_space_sep = uint_p[push_back_a(v)] >> *( *(space_p) >> uint_p[push_back_a(v)]);
     
  //rule<> r = *(space_p)>> uint_p[assign_a(id)] >> *(space_p) >>list_of_numbers_space_sep  >> *(anychar_p);

  rule<> r = list_of_numbers_space_sep; 


  //  reads  elements' total number (one  integer)
  getline(in_stream, str);
  if (parse(str.c_str(), uint_p[assign_a(n_elem)] >> *(space_p)   , space_p).full) 
    {
      // //cout << n_elem<<endl ;
    } 
      
      
  // reads  element   line : elem number, type of  elem, phys reg,  geom reg, number "n" of nodes in elem, "n" nodes
  //  and  put it in vector elem_values 
  while (getline(in_stream, str))  //   
    {
      
      if (  parse(str.c_str(), r) .full )
           
	{
	  
	  elem_values.push_back(v);
	                          
	
          v.clear();
	     
                             
                             
	} 


    } 
      
     
}



//  check if  $ELM  header  is  found 
void Read_MSH::find_elem_section(char const* str, ifstream& in_stream)
{
         
  string name;
    
  if  (parse(str, if_p("$")[(+alpha_p)[assign_a(name)]].else_p[nothing_p]  , space_p).full)
    {
      // //cout << name<< endl;
      if (name == "ELM")

    	{
    	  parse_elem_section(in_stream);
    	}
         
         
    }
}                                                                     







// scan  .msh  file  name,  line  by  line 
void Read_MSH::scan_input(string file_name)

{

  std::string str;

  std::ifstream in_stream (file_name.c_str());



  if ( !in_stream.good() )
    {
      std::cerr << "ERROR: Input file not good." 
		<< std::endl;
      //   error();
    }


  //cout << "/////////////////////////////////////////////////////////\n\n";
  //cout << "\t\tA space  separated list parser for Spirit...\n\n";
  //cout << "\t\t read  from .msh  file ELEM and NODES section \n\n";
  //
  //cout << "/////////////////////////////////////////////////////////\n\n";


  //getline(in_stream, str);
  //  while (str.empty() )
  //    {
  //      getline(in_stream, str);
  //    }
 


  while (getline(in_stream, str))  //   
    {
      read_data_section(str.c_str(),in_stream);
	   
    }





  //   while (getline(in_stream, str))  //   
  //     {
      

  //       find_node_section(str.c_str(),in_stream);
     
  //       //  find_elem_section(str.c_str(),in_stream);
  

  //     }

  //  //cout << "fine  node section"<< endl;

  //    while (getline(in_stream, str))  // 
  //     {

  //       find_elem_section(str.c_str(),in_stream);
  

  //     }



}





// utility method  to  read  element  lines ...
void  Read_MSH::get_data ( vector< vector<unsigned int> >& glob_elem_values )
{
  glob_elem_values = elem_values;
}



//   processes    elements lines :  extracts   (dim-1)D elements and associates   nodes of  each (dim-1)D element 
// to its own physical region number 

void Read_MSH::get_BC_info( vector<unsigned int>& BC_reg_ID )
{ // begin

  //  works  only  for  .msh file  format  v.1  !!!
  //  elem_id , el_type, phys_reg, geom_reg, number of  nodes, list of  nodes

  //  dim = n -> extract (n-1)D elements
  //  ( 1D elements for  2D problem,  2D  elements for 3D  problem)

  // put (n-1)D elements in  separate  vector
  //  dim_BC = dim - 1

  vector< vector<unsigned int> > BC_elem_line;  // 

  //  for (unsigned int i =0; i< elem_values.size();++i)
  //     {

  //       //  //cout <<  global_elem_values[i][1] << "   " ; 
  //       if (elem_values[i][1] == 1) 
  // 	{
  // 	  BC_elem_line.push_back(elem_values[i]);
  // 	}

  //     }


  int el_type;
 
  //
  //
  //  makes vector with lines of BC elements   
  for (unsigned int i =0; i< elem_values.size();++i)
    {

      el_type = elem_values[i][1];

      switch(dim) {
      case 1:
	// 1D  sim ................
	break;
      case 2:
	if ( (elem_values[i][1] == 1) || (elem_values[i][1] == 8) )
	  {
	    BC_elem_line.push_back(elem_values[i]);
	  }

	break;
      case 3:
	if ( (elem_values[i][1] == 2) || (elem_values[i][1] == 3) || (elem_values[i][1] == 9) || (elem_values[i][1] == 10) )
	  {
	    BC_elem_line.push_back(elem_values[i]);
	  }

	break;
      default:
	cout <<  "Error : wrong  value  of  dim"; 

      }

    }




  //  ***************************************************************************************

  //  BEWARE  !!!!  passed   from  calling  routine (BC  reg. ID  are read  from input file)

  //  ***************************************************************************************
  //   // creates vector of  BC regions

  //   vector<int> BC_reg_ID ;
  //   vector<int> :: iterator p;
  //   int id;
  //   id = 0;


  //   for (int i =0; i< BC_elem_line.size();++i)
  //     {

  //       id = BC_elem_line[i][2];
  //       p = find(BC_reg_ID.begin(), BC_reg_ID.end(), id );  // find  BC ID

  //       if (p == BC_reg_ID.end())   //  not  found
  // 	{
  // 	  BC_reg_ID.push_back(id);
  // 	  //  //cout<< BC_reg_ID[0]  << endl ;

  // 	}

  //       id = 0;
  //     }



  //  p_2 = BC_reg_ID.begin();

  //  //cout<< endl  << " BC  regions ....   "  << endl ;


  //  while (p_2 != BC_reg_ID.end())  //display 
  //    //cout<<*p_2++ << endl ;


  // **********************************
  //creates vector of BC_nodes vectors (one for each  BC  region)
  // **********************************

  vector< vector<unsigned int> > BC_nodes;
  vector< vector<unsigned int> > :: iterator p ;


  // initialization of BC_nodes vector
  vector<unsigned int> null_vector;
  for (unsigned int i =0; i< BC_reg_ID.size(); ++i)
    {
      BC_nodes.push_back(null_vector);
    }



  //  int  node1,  node2,  pos, id;
  unsigned int  pos, id, number_of_nodes ;

  id = 0;
  //  node1 = 0;
  //  node2 =0;

  vector<unsigned int> node_id;

  //  null_vector.push_back(0);

 
 
  // (BC_nodes[0]).push_back(node1);


  for (unsigned int i =0; i< BC_elem_line.size(); ++i)
    {

      id = BC_elem_line[i][2]  ;

      number_of_nodes = BC_elem_line[i][4];
      //   //cout<< " number_of_nodes  = " << number_of_nodes<< endl;

      for (unsigned int j =0; j<  number_of_nodes; ++j)
	{

	  //  //cout<< " (BC_elem_line[i][5+j]  = " << BC_elem_line[i][5+j]<< endl ;

	  node_id.push_back(BC_elem_line[i][5+j]) ;



	  //  node1 = BC_elem_line[i][5]  ;
	  // //cout<<  node1  << endl ;

	  //   node2 = BC_elem_line[i][6]  ;
	  // //cout<<  node2  << endl ;

	}
   
      pos = find_pos(id,BC_reg_ID);
 

      //   if (BC_nodes[pos].size() == 1)
      //        	{
      //        	  BC_nodes[pos].pop_back();  //  remove  initialization null
      //        	}

      //   //cout<< " node_id.size()  = " << node_id.size() << endl ;
      for (unsigned int i =0; i< node_id.size(); ++i)
	{ 

	  //  BC_nodes[pos].push_back(node1);
	  //  BC_nodes[pos].push_back(node2);

	  //	  //cout<< " node_id[i]  = " << node_id[i] << endl ;

	  BC_nodes[pos].push_back(node_id[i]);
	}



      //   p = BC_nodes.begin();
      //   p += pos+1;

      //   BC_nodes.insert(p,1,node_id);




      //   BC_nodes[pos].push_back(node_id);
      node_id.clear();

      //  node1 = 0;
      //  node2 =0;


    }

  //  //cout<< " BC_nodes.size()  = " << BC_nodes.size()  << endl ;
  //  //cout<< " BC_nodes[0].size()  = " << BC_nodes[0].size()  << endl ;




  //   for (int i =0; i< BC_nodes.size(); ++i)
  //     { 
  //       for (int j =0; j< BC_nodes[i].size(); ++j)
  // 	{

  // 	  //cout<< " BC_nodes[i][j]  = " <<BC_nodes[i][j] << endl;
  // 	}
  //     }



  //  unique  eliminates  repetitions !!!

  vector<unsigned int>::iterator new_end_BC;

  // = unique(test.begin(), test.end());

  // test.erase(new_end, test.end());



  for (unsigned int i =0; i< BC_nodes.size();++i)
    {
      new_end_BC = unique(BC_nodes[i].begin(), BC_nodes[i].end());

      BC_nodes[i].erase(new_end_BC, BC_nodes[i].end());

    }


  //   END creates vector of BC_nodes vectors (one for each  BC  region) 



  //  

  // ****************************************
  //  ONLY  IF   DIM = 3  !!!!!!!!
  // creates also  vector of BC_elements
  // ****************************************


  if  (dim == 3) 
    {  //if

      unsigned int el_id=0;
      id = 0;

      vector< vector<unsigned int> > BC_elements;

      // initialization of BC_nodes vector
      //  vector<int> null_vector;
      for (unsigned int i =0; i< BC_reg_ID.size(); ++i)
	{
	  BC_elements.push_back(null_vector);
	}


      for (unsigned int i =0; i< BC_elem_line.size(); ++i)
	{

	  id = BC_elem_line[i][2]  ;


	  el_id = BC_elem_line[i][0]  ;

	  pos = find_pos(id,BC_reg_ID);

	  BC_elements[pos].push_back(el_id);
	  el_id = 0;

	}



      //  creates  map  <BC_region, BC_elements >      //   for   dim = 3

      //  map<int, vector<int> >  // BoundCond_elem   in  .h  !!!

      for (unsigned int i =0; i< BC_elements.size();++i)
	{
	  BoundCond_elements.insert(make_pair(BC_reg_ID[i], BC_elements[i]) );
	  //   //cout << BC_reg_ID[i]<< endl ;
	  // //cout << BC_nodes[i][0]<< endl ;

	}

    }

  //   END  BC_elements

  // ****************************************





  //  creates  map  <BC_region, BC_nodes>

  //  map<int, vector<int> >  BoundCond   in  .h  !!!

  for (unsigned int i =0; i< BC_nodes.size();++i)
    {
      BoundCond.insert(make_pair(BC_reg_ID[i], BC_nodes[i]) );
      //   //cout << BC_reg_ID[i]<< endl ;
      // //cout << BC_nodes[i][0]<< endl ;

    }


}

//  END  get_BC_info() !!!



//  utility  function 
unsigned int  Read_MSH::find_pos( unsigned int  reg_id ,   vector<unsigned int>& BC_reg_ID )
{

  unsigned int  pos;
  unsigned int count = 0;
  bool found = false;


  while (found == false)  {
     
    if  (  (BC_reg_ID[count]) != reg_id)
      {
	count++; 
      }
    else 
      {
	found = true;
      }


  }  //  end  while

  pos =  count;
  return  pos;

}





// public  method  to  get  map    <BC_region, BC_nodes> 
void  Read_MSH::get_BC_data (    map<unsigned  int, vector<unsigned int> >&  BoundCond_map    )
{

  BoundCond_map = BoundCond;
 
}


// public  method  to  get  map  <BC_region, BC_nodes> and map   <BC_region, BC_elements> 
void  Read_MSH::get_BC_data (map<unsigned int, vector<unsigned int> >&  BoundCond_map, map<unsigned int, vector<unsigned int> >&  BoundCond_el_map      )
{

  BoundCond_map = BoundCond;

  //   if  dim  =  3  : BoundCond_elem 

  BoundCond_el_map = BoundCond_elements;

 
}


// public  method  to  get  map  <phys_region, elements>  
void  Read_MSH::get_elem_data (map<unsigned int, vector<unsigned int> >& PhysReg_elements_map)
{ 
      
PhysReg_elements_map = PhysReg_elements;  
}





// get  physical  regions  for  each  element
void Read_MSH::get_physical_elem(vector<unsigned int>& phys_reg_ID)
{ //get_phys

	
  //  get   elements  associated  to  each  physical  region

  //  works  only  for  .msh file  format  v.1  !!!
  //  elem_id , el_type, phys_reg, geom_reg, number of  nodes, list of  nodes

 

  vector< vector<unsigned int> > phys_elem_line;  //  ->  rename  ok!   BC_elem_line !

  //  for (int i =0; i< elem_values.size();++i)
  //     {

  //       //  //cout <<  global_elem_values[i][1] << "   " ; 
  //       if (elem_values[i][1] == 1) 
  // 	{
  // 	  BC_elem_line.push_back(elem_values[i]);
  // 	}

  //     }

 
  unsigned int el_type;

 
  // **********************************************

  //  idea  !:  work  on    vector<Element>  list_elements   instead   of phys_elem_lines
  //   list_elements  has  already  the  elements with dim  =  sim_dim !!
  // **********************************************

  for (unsigned int i =0; i< elem_values.size();++i)
    {

      el_type = elem_values[i][1];

      switch(dim) {
      case 1:
	// 1D  sim ................elem_values[i][1] == 1) || (elem_values[i][1] == 8)
	break;
	
      case 2:
	if ( (el_type == 2) || (el_type == 3) || (el_type == 9) || (el_type == 10) )
	  {
	    phys_elem_line.push_back(elem_values[i]);
	  }

	break;
      case 3:
	if ( (el_type == 4) || (el_type == 5) || (el_type == 6) || (el_type == 7)|| (el_type == 11)|| 
	     (el_type == 12)|| (el_type == 13) || (el_type == 14))
	  {
	    phys_elem_line.push_back(elem_values[i]);
	  }

	break;
	
      default:
	//cout <<  "Error : wrong  value  of  dim"; 
	break;

      }

    }


  // ****************************************
 
  // creates   vector of BC_elements
  // ****************************************


  unsigned int id;  // physical reg id
  unsigned int  pos;
  unsigned int el_id=0;
  id = 0;

  vector<unsigned int> null_vector;

  vector< vector<unsigned int> > phys_elements;  //get_physical_elem

  // initialization of phys_elements vector
  //  vector<int> null_vector;
  for (unsigned int i =0; i< phys_reg_ID.size(); ++i)
    {
      phys_elements.push_back(null_vector);
    }


 //  for (unsigned int i =0; i< phys_elem_line.size(); ++i)
//     {

//       id = phys_elem_line[i][2]  ;


//       //	  if  (id != 100)     //  OLD   FOR  PYTHON SCRIPT  !!!
//       //	    {
//       //   el_id = phys_elem_line[i][0]  ;  //  !!!!!!!!!!!!!!!!!    ATTENZIONE  !!!!  SE SALTA  NUM ORDINE  ELEMENTI !!??
//       el_id = (i +1);  // renumber elements (BC_elements excluded !)  (i  starts  from zero !!)
	   
//       //   if  (el_id != (i+1))
//       // 		{ //cout <<  "******** WARNING !!!! : i-th element id != i (break in  element list of  GMSH)"<< el_id <<endl<< endl; }  


//       pos = find_pos(id,phys_reg_ID);

//       phys_elements[pos].push_back(el_id);
//       //	    }
//       el_id = 0;

//     }

  // *******************************************************
  //***************************************************
  //  with  vector<Element>  list_elements;

  unsigned int  physic_id;
  physic_id = 0;

 for (unsigned int i =0; i< list_elements.size(); ++i)
    {
      physic_id = list_elements[i].phys_id;

      pos = find_pos(physic_id,phys_reg_ID);
      el_id = (i +1);
      phys_elements[pos].push_back(el_id);
      //	    
      el_id = 0;
    }

  // *******************************************************
  //************************************************




  //  ??????????????????????????
  //  creates  map  <phys_region, phys_elements >      //  

  //  map<int, vector<int> >  // PhysReg_elem   in  .h  !!!

  for (unsigned int i =0; i< phys_elements.size();++i)
    {
      PhysReg_elements.insert(make_pair(phys_reg_ID[i], phys_elements[i]) );
      //   //cout << BC_reg_ID[i]<< endl ;
      // //cout << BC_nodes[i][0]<< endl ;

    }
	
  //  ??????????????????????????


	



  //   *********************************************************************
  //    map   < element number , phys_reg >  
 	 
  //	map<int, int >   elem_region_map;  in  .h  !!!

  //	//cout << endl <<  "phys_elements.size()" << phys_elements.size()  <<endl<< endl ;
	
  for (unsigned int i =0; i< phys_elements.size();++i)
    {
      for (unsigned int j =0; j< phys_elements[i].size();++j)
	{
	  
	  elem_region_map.insert(make_pair( phys_elements[i][j],phys_reg_ID[i] ) );
	}
    }
        
        
  ////cout << endl <<  "content of  map elem / PhysReg " << endl<< endl ;
       


  unsigned int  temp;
  map <unsigned int, unsigned int>   :: iterator  p_el_data;



  //  *********************************************
  //    DELETE  !!!
      
  //  display  map
    
  for (unsigned int i =0; i< phys_elements.size();++i)
    { // display

      for (unsigned int j =0; j< phys_elements[i].size();++j)
        {
        
	  p_el_data = elem_region_map.find(phys_elements[i][j]);
    
	  if  (p_el_data != elem_region_map.end() )
    
	    {  
	      temp  =  (p_el_data -> second) ;
	      ////cout << endl <<  "BC #"<< phys_reg_numbers[i] << "   " << endl; 
	      //for (int i =0; i< temp_3.size();++i)
	      //    	      {
	      //    
	      //    		if ( (i % 3) ==  0)     //cout << endl; 
    
	      // //cout << phys_elements[i][j]<< "   "  << temp ;  //<< endl; 
    
	      //   }
	      ////cout << endl; 
    
    
	    }
    
	  else
	    cout  <<  "error elem_region_map" << endl ;
    	  	  
	}
  
      //   }    
        
      //  
					
    }

  // ***************************************


  //  
 //    DELETE  !!!
  // ****************************************


  //  END  get_physical_elem  () !!!
	
}




//  to  make  private !!!
// public  method  to  get  map  <phys_region, elements>  
void  Read_MSH::get_elem_phys_map (map<unsigned int,unsigned int> &elem_phys  )
{ 
      
  elem_phys =   elem_region_map ;
   
}



// write  .xda  file
void  Read_MSH::write_xda ( )
{ // xda
  // 
 
  //  vector<Element>  list_elements;
  list_elements.clear();
  num_elem_per_type.clear();
  unsigned int	count = 0;
  count =0;



  // ************************************************************************

  //  write  .xda file  
  // ************************************************************************

  { // write .xda



    // elem_type_conversion from GMSH_type to  xda (libmesh) type 

    // ************************************************************
    //  xda (libmesh)  ELEM TYPES 
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

    //**************************************************************************
    // GMSH  ELEM TYPES
    //**************************************************************************


    // 1 Line (2 nodes).
    //  2 Triangle (3 nodes).
    //  3 Quadrangle (4 nodes). 
    // 4 Tetrahedron (4 nodes). 
    // 5 Hexahedron (8 nodes). 
    // 6 Prism (6 nodes). 
    // 7 Pyramid (5 nodes).
    //  8 Second order line (3 nodes: 2 associated with the vertices and 1 with the edge). 
    // 9 Second order triangle (6 nodes: 3 associated with the vertices and 3 with the edges). 
    // 10 Second order quadrangle (9 nodes: 4 associated with the vertices, 4 with the edges and 1 with the face).
    //  11 Second order tetrahedron (10 nodes: 4 associated with the vertices and 6 with the edges). 
    // 12 Second order hexahedron (27 nodes: 8 associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume). 
    // 13 Second order prism (18 nodes: 6 associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
    //  14 Second order pyramid (14 nodes: 5 associated with the vertices, 8 with the edges and 1 with the quadrangular face). 
    // 15 Pounsigned int (1 node).




	
    // elem_type_conversion from GMSH_type to  xda (libmesh) type 
    //  ??
    //  gmsh_elem_type  is  a  vector
// 	for (int i =0; i< gmsh_elem_type.size();++i)
// 	  { switch(gmsh_elem_type[i]) {
// case 
//  elem_type.push_back[...] 

    elem_type.clear();

    for (unsigned int i =0; i< gmsh_elem_type.size();++i)
      { //for
	switch(gmsh_elem_type[i]) {
	  //   switch(gmsh_elem_type) {

	case 1:
	  // 1D  sim ................  line
	  break;
	
	case 2:
	  //   TRI3       // 3
	  //    elem_type = 3; 
	  elem_type.push_back(3);
	  break;

	case 3:
	  //QUAD4,      // 5
	  //   elem_type = 5; 
	  elem_type.push_back(5);
	  break;

	case 4:
	  //    TET4,       // 8 
	  //     elem_type = 8;  
	  elem_type.push_back(8);
	  break;

	case 5:
	  // HEX8,       // 10
	  //  elem_type = 10 ; 
	  elem_type.push_back(10);
	  break;

	case 6:
	  //  PRISM6,     // 13
	  //   elem_type = 13;    
	  elem_type.push_back(13);
	  break;

	case 7:
	  //  PYRAMID5,   // 16
	  //  elem_type = 16 ;
	  elem_type.push_back(16);
	  break;

	case 8:
	  // 1D   second order line  ........
	  //  elem_type = ; 
     
	  break;

	case 9:
	  //  TRI6,       // 4
	  //   elem_type = 4 ; 
	  elem_type.push_back(4);
	  break;


	case 10:
	  //    QUAD9,      // 7
	  //  elem_type = 7; 
	  elem_type.push_back(7);
	  break;

	case 11:
	  // TET10,      // 9
	  //  elem_type = 9  ; 
	  elem_type.push_back(9);
	  break;


	case 12:
	  //      HEX27,      // 12
	  //     elem_type = 12; 
	  elem_type.push_back(12);
	  break;


	case 13:
	  //  PRISM18,    // 15
	  //  elem_type = 15 ; 
	  elem_type.push_back(15);
	  break;

	case 14:
	  // Second order pyramid (14 nodes)  not  present  in  xda format !!!!
	  //  elem_type = ;      
	  break;

	case 15:
	  // 15 Point (1 node).
	  //   elem_type = ;    ........................  0 D  
	  break;

	
	default:
	  cout <<  "Error : wrong  value  of gmsh_elem_type  "; 

	}


      }


    string fname;
    fname = "in.xda";

 

    // Open the output file stream
    std::ofstream out (fname.c_str());
 
    assert (out.good());



    // Syntax of  file  .xda

    //    see 

    //DEAL 003:003
    //30	 # Num. Elements  =  num of  elements
    //20	 # Num. Nodes  =  num of  nodes
    //90	 # Sum of Element Weights= sum of nodes of  all the  elements = Sigma(i=Nel)[ nodes(i)]
    //0	 # Num. Boundary Conds.
    //65536	 # String Size (ignore)
    //1	 # Num. Element Blocks. = "number  of  mesh  blocks" =  1 !
    //3 	 # Element types in each block. = vector of  element types (TRI = 3) (see  elem_type.h)
    //30 	 # Num. of elements in each block = total number of elements of a given type
    //Id String
    //Title String
    //  for  (all elements)
    // { list  of  nodes}
    //  for  (all  nodes) 
    // { list of  node coordinates}
    //


    // declared in .h 	 
    num_bc =0;
    //  num_mesh_block = 1 ;

    num_mesh_block =  elem_type.size();

 /**
     * A mesh block by definition contains
     * only a single type of element.
     *
     * @return The number of mesh blocks.
     */

    //  elem_type = 3;  //4 ; //3;
    // num_elem_per_type  must  be  a  vector  !!

/**
     * The size of each element block is
     * the total number of a given type of
     * element in the mesh.
     *
     * @return The vector of block sizes
     */




//   num_elem_per_type =  num_of_elem ;
    //  num_elem_per_type =  to  be  calculated  !!


    // ************************************************************
    //  calculation of  list  of  elements( per  type block) 
    //  and  of  num_elem_per_type =  num  of  elements of  each of the  types in elem_type

   
 for (unsigned int i =0; i< gmsh_elem_type.size();++i)
   { //for i

     for (unsigned int j =0; j< All_elements.size();++j)
       {//for j


	 if (All_elements[j].type == gmsh_elem_type[i])
	   {
	     // //cout << "All_elements[j].type " << All_elements[j].type << endl;

	     list_elements.push_back(All_elements[j]);
	     count++;

	   }

       }

     num_elem_per_type.push_back(count);
     count = 0;


   }






 // ************************************************************







    //  int el_id = 0 ;
    //   double  reg_id= 0.0;  // = 101.0; 

    //    int n_elem = 30 ; //33313;
    //cout << " num_of_nodes" <<  num_of_nodes;
    


    //  ***********************************************************
    // HEADER  of  .xda  file
    // ************************************************************
    //
    out << "DEAL 003:003\n";
    out << num_of_elem  <<   "	 # Num. Elements  \n";
    out << num_of_nodes  <<   "	 # Num. Nodes \n";
    out <<  el_weight <<   " 	 # Sum of Element Weights \n";
    out <<  num_bc  <<   "	 # Num. Boundary Conds. \n";
    out << "65536	 # String Size (ignore) \n";
    out <<  num_mesh_block  <<   "	 # Num. Element Blocks.	\n";
    // elem_type = vector of  elem types present in  mesh
    for (unsigned int i =0; i< elem_type.size();++i)
      {

	out <<  elem_type[i]<< "   ";

      }
    //   out <<  elem_type << " 	 # Element types in each block. \n";
    out << " 	 # Element types in each block. \n";

    //   out << num_elem_per_type << "  	 # Num. of elements in each block. \n"; 

    for (unsigned int i =0; i< elem_type.size();++i)
      {

	out <<  num_elem_per_type[i]<< "   ";

      }
    out << "  	 # Num. of elements in each block. \n"; 

    out << "Id String  \n";
    out <<"Title   String  \n";
    //
    //



//  element  section 
    // ******************************
    unsigned int size_element_list;
    unsigned int nodes_size;
    unsigned int  node;

    size_element_list=  list_elements.size();


    //   for (unsigned int i =0; i<list_elements.size() ; i++)
    for (unsigned int i =0; i<size_element_list  ; i++)
      {
	nodes_size=list_elements[i].nodes.size();

      //	for (unsigned int j =0; j<list_elements[i].nodes.size(); j++)
	for (unsigned int j =0; j< nodes_size; j++)

	  {

	    node = list_elements[i].nodes[j] - 1;
	    //   out << ( list_elements[i].nodes[j] - 1)   << "  " ;
	    out << node  << "  " ;

	  }

	out << "\n";

      }

 // ******************************



	
//     //  element  section 
//     for (unsigned int el_id =0; el_id<num_of_elem ; el_id++)
//       {
// 	for (unsigned int i =0; i<elem_nodes[el_id].size() ; i++)
// 	  {



// 	    out << (elem_nodes[el_id][i]-1)    << "  ";
	   
// 	    //    out << std::setprecision(10) << reg_id   << "  	   # Values   ";
	  
// 	  }

// 	out << "\n";

//       }




    // nodes  section 

    for (unsigned int node_id =0; node_id<num_of_nodes ; node_id++)
      {
	for (unsigned int i =0; i<node_coord[node_id].size() ; i++)
	  {
	     
	    out << std::setprecision(10) << node_coord[node_id][i]   << "  ";
	    //  out <<  node_coord[node_id][i]   << "  ";

	   
	    //    out << std::setprecision(10) << reg_id   << "  	   # Values   ";
	  
	  }

	out << "\n";

      }




  }


  // ************************************************************************


  // END  

   
}



// get  element data  for  each  element
void Read_MSH::get_elem_nodes()
{ //get_elem_nodes

  Element  current_element;  //   struct *************
  All_elements.clear();

  unsigned int el_type;

  vector< vector<unsigned int> > elem_line;
  vector<unsigned int> :: iterator p;
  gmsh_elem_type.clear();

  for (unsigned int i =0; i< elem_values.size();++i)
    { //for

      el_type = elem_values[i][1];

    
      //    gmsh_elem_type = el_type;



      // 

      //   vector<int> BC_reg_ID ;
      //  vector<int> :: iterator p;
      // if (! find (el_type) in  gmsh_elem_type_vector
      //      gmsh_elem_type_vector.push_back(el_type)


      //  ******  put in  switch *****************

   //    p = find(gmsh_elem_type.begin(), gmsh_elem_type.end(),el_type  );  // 

//       if (p == gmsh_elem_type.end())   //  not  found
//  	{
//  	  gmsh_elem_type.push_back(el_type);
// 	  // 	  //  //cout<< BC_reg_ID[0]  << endl ;

// 	}


      //.......................


      switch(dim) {
      case 1:
	// 1D  sim ................elem_values[i][1] == 1) || (elem_values[i][1] == 8)
	break;
	
      case 2:
	if ( (el_type == 2) || (el_type == 3) || (el_type == 9) || (el_type == 10) )
	  {
	    elem_line.push_back(elem_values[i]);
	    p = find(gmsh_elem_type.begin(), gmsh_elem_type.end(),el_type  );

	    if (p == gmsh_elem_type.end())   //  not  found
	      {
		gmsh_elem_type.push_back(el_type);

		////cout<< "el_type*************"<<  el_type   << endl ;

	      }


	  }

	break;
      case 3:
	if ( (el_type == 4) || (el_type == 5) || (el_type == 6) || (el_type == 7)|| (el_type == 11)|| 
	     (el_type == 12)|| (el_type == 13) || (el_type == 14))
	  {
	    elem_line.push_back(elem_values[i]);
	    p = find(gmsh_elem_type.begin(), gmsh_elem_type.end(),el_type  );

	    if (p == gmsh_elem_type.end())   //  not  found
	      {
		gmsh_elem_type.push_back(el_type);

		////cout<< "el_type*************"<<  el_type   << endl ;

	      }



	  }

	break;
	
      default:
	cout <<  "Error :    wrong  value  of  dim"; 

      }

    //   if (p == gmsh_elem_type.end())   //  not  found
//  	{
//  	  gmsh_elem_type.push_back(el_type);

// 	  //cout<< "el_type*************"<<  el_type   << endl ;

// 	}



    }



  unsigned int new_node_id;
  double node_id =0.0 ; 
  unsigned int count = 0;
  unsigned int  number_of_nodes_in_elem = 0;
  vector<unsigned int> node_id_list;

  for (unsigned int i =0; i< elem_line.size(); ++i)
    {  

      number_of_nodes_in_elem = elem_line[i][4];

    

      for (unsigned int j =0; j<  number_of_nodes_in_elem; ++j)
	{


	  node_id  = (double)(elem_line[i][5+j]);    //  ??????   change  node_label and  node_entry to  int !!!
	  for (unsigned int k =0; k<num_of_nodes ; ++k)
	    { 
	      //  //cout << " node_label[k]" <<  node_label[k]<< endl;
	      if  (node_label[k] ==  node_id)
		{ new_node_id = k+1; }
	    }
	  node_id_list.push_back( new_node_id) ;
	  //	 //cout << "new_node_id " << new_node_id<< endl ;

	  //	 node_id_list.push_back(elem_line[i][5+j]) ;

	  count ++;
	  new_node_id = 0;
	  node_id  =0;

	}

      //*********************

      current_element.nodes = node_id_list;
      current_element.type =  elem_line[i][1];

      current_element.phys_id =  elem_line[i][2];

      All_elements.push_back(current_element);
      // //cout << "******************" << current_element.type<< endl ;
      current_element.nodes.clear();
      current_element.type = 0;


      //****************



      elem_nodes.push_back(node_id_list);

      node_id_list.clear();

    }



  el_weight = count ;
  num_of_elem = elem_nodes.size();




} //  end  get_elem_nodes




//  parser of  NODES  section of  .msh  file  (version 1)
void Read_MSH::parse_node_section(ifstream& in_stream )
{ // parse nodes

  bool end_nodes = false;
 

  string  str;
  vector<double> v;
  v.clear();
     
  //  rule  for node  line  (list  of  real)
  rule<>list_of_numbers_space_sep = real_p[push_back_a(v)] >> *( *(space_p) >> real_p[push_back_a(v)]);
     
  //rule<> r = *(space_p)>> uint_p[assign_a(id)] >> *(space_p) >>list_of_numbers_space_sep  >> *(anychar_p);

  rule<> r = list_of_numbers_space_sep; 



  //  reads  nodes  total number (one  integer)
  getline(in_stream, str);
  if (parse(str.c_str(), uint_p[assign_a(num_of_nodes)] >> *(space_p)   , space_p).full) 
    {
      //cout << " NUMERO  NODI *********"<< num_of_nodes  <<endl ;
    } 
      
      
  // reads  NODE   line : NODE  number, 3  coordinates
  //  and  put it in vector node_entry
 
  while  (getline(in_stream, str))      //   
    {
      
      if (  parse(str.c_str(), r) .full )
           
	{
	  
	  node_entry.push_back(v);
	  //     //cout << v[0] <<endl ;                
	
          v.clear();
	     
                             
                             
	} 

      if (parse(str.c_str(), str_p("$ENDNOD") ).full) 

	//    if (  parse(str.c_str(), r) .full )
           
	{

	  //cout << " BREAK   "  <<endl ;  
	  end_nodes = true; 
  	  break;
	}

      // if     trovato $ENDNOD 
      //   {exit
      //      }


    } 
      
     
} // end  parse nodes




//  check if  $NOD  header  is  found 
void Read_MSH::find_node_section(char const* str, ifstream& in_stream)
{ //
         
  string name;
    
  if  (parse(str, if_p("$")[(+alpha_p)[assign_a(name)]].else_p[nothing_p]  , space_p).full)
    {
      //    //cout << name<< endl;
      if (name == "NOD")

    	{
    	  parse_node_section(in_stream);
    	}
         
         
    }
}                                                                     



// get  node  data  for  each  node
void Read_MSH::get_nodes_coord()
{ //get_nodes_coord

  vector<double> temp;

  node_label.clear();
  for (unsigned int i =0; i<num_of_nodes ; ++i)
    {  

      //     //cout << "(node_entry[i][0] " <<node_entry[i][0]<< endl ;
      node_label.push_back(node_entry[i][0]) ; // ????


      //   //cout << " node_label[i]" <<  node_label[i];

      for (unsigned int j =0; j<3 ; ++j)
	{ 

	  temp.push_back(node_entry[i][j+1]);
	}

      //  node_coord[i].push_back(temp);
      node_coord.push_back(temp);


      temp.clear();


    }


} // end get_nodes_coord






void Read_MSH::read_data_section(char const* str,ifstream& in_stream )
  // void InputParser::read_data_section(char const* str)


  // void read_data_section(char const* str,ifstream& in_stream )

{

  string name;

  if  (parse(str, if_p("$")[(+alpha_p)[assign_a(name)]].else_p[nothing_p]  , space_p).full)
    {
      //cout << name<< endl;
      if (name == "ELM")
	//	if (name == name_data )
	{
	  parse_elem_section(in_stream);
	}

      else if  (name == "NOD")
	{
	  //  parse_2(in_stream);
	  parse_node_section(in_stream);

	}



    }



}




// ************************************************************************

//  write  .xta file (data file for  meshdata (elem_data.xta )  
// 
// ************************************************************************



void  Read_MSH::write_xta()
{ // xta
  //



  // write .xta


  map <unsigned int, unsigned int>   :: iterator  p_el_data;


  string fname;
  fname = "elem_data.xta";

 

  // Open the output file stream
  std::ofstream out (fname.c_str());
 
  assert (out.good());



  // Syntax of  file  .xta
  //
  //  # Data description
  // REAL	 # type of values
  // 0	 # No. of nodes for which data is stored
  // n_elem (excluded  reg 100 !!!) 	 # No. of elements for which data is stored


  // for  all_elements
  //
  // el_id -1	 # Foreign element id
  //  1	 # vector length
  //   reg_ID	 # Values
  // 

  //  int el_id = 0 ;
  double  reg_id= 0.0;  // = 101.0; 

  // int n_elem = 150963; // 524;  // 248 ; //30 ; //33313;

  //  num_of_elem is  a  member data
  
  out << " # Data description \n";
  out << "REAL	 # type of values     \n";
  out << "0	 # No. of nodes for which data is stored     \n";
  out << num_of_elem  << "	 # No. of elements for which data is stored    \n";

  for (unsigned  int el_id =0; el_id<num_of_elem; el_id++)
    {

      p_el_data =elem_region_map.find((el_id+1)  );

      //   p_el_data = elem_phys_map.find((el_id+1)  );

      //   if  (p_el_data != elem_phys_map.end() )

      if  (p_el_data != elem_region_map.end() )
    
	{ 
	  reg_id	  =  (double) (p_el_data -> second) ;
	}
      else
	cout  <<  "error  write_xta "<< endl  ;

      out << el_id   << "     # Foreign element id ";

      out << "\n";
      out << "1	 # vector length    ";
      out << "\n";
      out << std::setprecision(10) << reg_id   << "  	   # Values   ";
      out << "\n";


    }

}


// ************************************************************************
// ************************************************************************

// ************************************************************************

// END   write  .xta file  (elem_data.xta )
// ************************************************************************








// write  mesh  and  meshdata from  .xda and  .xta  files

void   Read_MSH::read_mesh_and_data(Mesh& mesh, MeshData_elements&  mesh_data )

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





