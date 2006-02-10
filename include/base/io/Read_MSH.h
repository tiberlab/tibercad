/*=============================================================================
    Copyright (c) 2002-2003 Joel de Guzman
    http://spirit.sourceforge.net/

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
///////////////////////////////////////////////////////////////////////////////
//
//  This sample demontrates a parser for a comma separated list of numbers
//  This is discussed in the "Quick Start" chapter in the Spirit User's Guide.
//
//  [ JDG 5/10/2002 ]
//
///////////////////////////////////////////////////////////////////////////////
#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/dynamic.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include "mesh_data_elements.h"
#include "mesh_data.h"
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"


///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace boost::spirit;

///////////////////////////////////////////////////////////////////////////////
//
//  Our comma separated list parser
//
///////////////////////////////////////////////////////////////////////////////

class Read_MSH{

  // private:

  unsigned int dim;

  map<unsigned int, vector<unsigned int> >  BoundCond;

  map<unsigned int, vector<unsigned int> >  BoundCond_elements;
  
  map<unsigned int, vector<unsigned int> > PhysReg_elements;

  map<unsigned int, unsigned int >   elem_region_map;



  vector< vector<unsigned int> > elem_values  ;

  vector< vector<double> > node_entry  ;

 vector< vector<double> >  node_coord ;



 vector< vector<unsigned int> > elem_nodes  ;

 vector<double> node_label;


 unsigned int  el_weight, num_of_elem, num_of_nodes,  num_bc,  num_mesh_block ; //, num_elem_per_type;

 vector<unsigned int> gmsh_elem_type;

 vector<unsigned int> elem_type;

 vector<unsigned int> num_elem_per_type ;



 struct Element {

   vector<unsigned int> nodes;
   unsigned int type;
   unsigned int phys_id  ;

 };


 vector<Element> All_elements;

 vector<Element>  list_elements;  //  list for   .xta  file  





  void Read_MSH::initialize_vectors();

  
  void Read_MSH::parse_elem_section(ifstream& in_stream );


  void Read_MSH::find_elem_section(char const* str, ifstream& in_stream);


  void Read_MSH::scan_input(string file_name);



  unsigned int  Read_MSH::find_pos( unsigned int  reg_id ,   vector<unsigned int>& BC_reg_ID );

  void Read_MSH::get_BC_info( vector<unsigned int>& BC_reg_ID );

// get  physical  regions
   void Read_MSH::get_physical_elem(vector<unsigned int>& phys_reg_ID);

// write  .xda  file
   void  Read_MSH::write_xda ( );

   
// get  element data  for  each  element
   void Read_MSH::get_elem_nodes();


//  parser of  NODES  section of  .msh  file  (version 1)
   void Read_MSH::parse_node_section(ifstream& in_stream );

//  check if  $NOD  header  is  found 
   void Read_MSH::find_node_section(char const* str, ifstream& in_stream);

// get  node  data  for  each  node
   void Read_MSH::get_nodes_coord();


   void Read_MSH::read_data_section(char const* str,ifstream& in_stream );

//  write  .xta file (data file for  meshdata (elem_data.xta )  
   void  Read_MSH::write_xta();

//  write  mesh  and  meshdata from  .xda and  .xta  files
   void   Read_MSH::read_mesh_and_data(Mesh& mesh, MeshData_elements&  mesh_data );



 public:

//  MSH_parser(string filename , vector<int>& BC_reg_ID, int sim_dim );
  
//  Read_MSH(string filename, vector<int>& phys_reg_ID, vector<int>& BC_reg_ID, int sim_dim  );

   Read_MSH(string filename, vector<unsigned int>& phys_reg_ID, vector<unsigned int>& BC_reg_ID, unsigned int sim_dim, Mesh& mesh, MeshData_elements&  mesh_data);


  ~Read_MSH();




  // void  MSH_parser::get_data ( vector< vector<double> >& glob_reg_values,  vector<int>& glob_reg_id,  vector<string>& glob_mat);

  void  Read_MSH::get_data ( vector< vector<unsigned int> >& glob_elem_values);

  //  void  Read_MSH::get_BC_data (   map<int, vector<int> >&  BoundCond_map);

  void  Read_MSH::get_BC_data (   map<unsigned int, vector<unsigned int> >&  BoundCond_map);

  void  Read_MSH::get_BC_data (   map<unsigned int, vector<unsigned int> >&  BoundCond_map, map<unsigned int, vector<unsigned int> >&  BoundCond_el_map   );

  void  Read_MSH::get_elem_data (map<unsigned int, vector<unsigned int> >& PhysReg_elements_map);

  void  Read_MSH::get_elem_phys_map (map<unsigned int,unsigned int> &elem_phys  );


};



