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
//  Class  for  parsing of  GMSH  .MSH file (v. 1 and 2 ) :   get  nodes of  BC regions 
// (1D  physical regions for 2D simulations and 2D physical regions for 3D simulations)
//
//   and  creates .xta and  .xda  files for  mesh and  meshdata in Libmesh 
//
// -------------------------------------------
//  MUST   include  SPIRIT  libraries !
// -----------------------------------------------
//
///////////////////////////////////////////////////////////////////////////////


//  constructor:  scan  .msh  file,  parse it, get info on BC,  writes .xda and  
// .xta  file for  Libmesh
//  needs physical region ID,  bound. conditions region  ID (2D or 1D), dim. of  
// simulation 
//Read_MSH::Read_MSH(string filename, vector<unsigned int>& phys_reg_ID, 
//                   vector<unsigned int>& BC_reg_ID, unsigned int sim_dim, 
//                   Mesh& mesh, MeshData_elements&  mesh_data)

Read_MSH::Read_MSH(string filename, unsigned int sim_dim, 
                   Mesh& mesh, MeshData_elements&  mesh_data)

{

  dim =  sim_dim; 
  version = 1;  //  by default it's version 1 of .msh file  format
  // in scan_input it is override if there is a MeshFormat section

  initialize_vectors();   //  to  be  continued...


  // call parser of  .msh  file
  scan_input(filename);

  //  get_BC_info(BC_reg_ID); //  ->  after   write_xda(); 


  // put  nodes and  elements in data structures
  get_nodes_coord();
  get_elem_nodes();


  // write .xda file of mesh (elements and  nodes :  
  // only elements with correct dim are considered)
  write_xda();


  // extract nodes which  belong to  bound cond reg. (DIM-1 phys.reg.) and 
  // associate them with 
  // user defined BC regions
  //  get_BC_info(BC_reg_ID);

  get_BC_info();

  // associate elements to corrispondent physical region (as defined in user's 
  // phys region settings)
  ////  get_physical_elem(phys_reg_ID);
  get_physical_elem();


  // write  .xta file with  meshdata info for  Libmesh (here ID of  phys regions)  
  write_xta();

  // read .xta and  .xda files  in  Libmesh  data structures
  read_mesh_and_data(mesh,mesh_data );


}


Read_MSH::~Read_MSH()
{
}


void Read_MSH::initialize_vectors()
{

  elem_values.clear();
  // ......................
  
}



//  parser of  ELEMENT  section of  .msh  file  (version 1 and 2)
//
//  ****  VERSION  2.0 ***
// element   line : elem number, type of  elem, number of tags (default=3) ,
// tag1=physical_region, tag2=geom_region, tag3= 0 ("mesh_partition"), list of the nodes
//
void Read_MSH::parse_elem_section(ifstream& in_stream )
{
     
  unsigned int n_elem;
  string  str;
  vector<unsigned int> v;
  v.clear();
     
  //  rule  for element  line  (list  of  integers uint)
  //  rule<>list_of_numbers_space_sep = uint_p[push_back_a(v)] >> 
  //    *( *(space_p) >> uint_p[push_back_a(v)]);

  // 24.4.08
  rule<>list_of_numbers_space_sep = uint_p[push_back_a(v)] >> 
    *( *(space_p) >> uint_p[push_back_a(v)])  >> *(anychar_p)  ;


     
  //rule<> r = *(space_p)>> uint_p[assign_a(id)] >> *(space_p) >>
  // list_of_numbers_space_sep  >> *(anychar_p);

  rule<> r = list_of_numbers_space_sep; 


  //  reads  elements' total number (one  integer)
  getline(in_stream, str);
  if (parse(str.c_str(), uint_p[assign_a(n_elem)] >> *(space_p), space_p).full) 
  {
    // cout << n_elem<<endl ;
  } 
      
      
  // reads  element   line : elem number, type of  elem, phys reg,  geom reg, 
  // number "n" of nodes in elem, "n" nodes
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




// NOT  USED  !!??
//  check if  $ELM  header  is  found 
void Read_MSH::find_elem_section(char const* str, ifstream& in_stream)
{
         
  string name;
    
  if  (parse(str, if_p("$")[(+alpha_p)[assign_a(name)]].else_p[nothing_p], space_p).full)
  {
    // cout << name<< endl;
    if (name == "ELM")

    {
      parse_elem_section(in_stream);
    }
         
         
  }
}                                                                     



// **************************************  
// scan  .msh  file  name
void Read_MSH::scan_input(string file_name)

{

  std::string  Nodes_label, Elements_label ,label ;

  std::ifstream in_stream (file_name.c_str());

  int file_type,data_size;
  double version_number;


  if ( !in_stream.good() )
  {
    std::cerr << "ERROR: Input file not good." 
              << std::endl;
    //   error();
  }

  in_stream >> label;

  // **********************  version 2
  if (label == "$MeshFormat")
  {
    parse_meshformat_section(in_stream);
    in_stream >> label;
 

  }

  //  in_stream >> Nodes_label;
  parse_node_section(in_stream);

 

  in_stream >> Elements_label;
  parse_elem_section(in_stream);

}







//**********************************************



  // utility method  to  read  element  lines ...
void  Read_MSH::get_data ( vector< vector<unsigned int> >& glob_elem_values )
{
  glob_elem_values = elem_values;
}



//   processes    elements lines :  extracts   (dim-1)D elements and associates 
//  nodes of  each (dim-1)D element 
// to its own physical region number 

// void Read_MSH::get_BC_info( vector<unsigned int>& BC_reg_ID )
void Read_MSH::get_BC_info()

{ // begin

  //  works  only  for  .msh file  format  v.1  !!!
  //  elem_id , el_type, phys_reg, geom_reg, number of  nodes, list of  nodes

  //  dim = n -> extract (n-1)D elements
  //  ( 1D elements for  2D problem,  2D  elements for 3D  problem)

  // put (n-1)D elements in  separate  vector
  //  dim_BC = dim - 1

  vector< vector<unsigned int> > BC_elem_line;  // 

  //  for (int i =0; i< elem_values.size();++i)
  //     {

  //       //  cout <<  global_elem_values[i][1] << "   " ; 
  //       if (elem_values[i][1] == 1) 
  // 	{
  // 	  BC_elem_line.push_back(elem_values[i]);
  // 	}

  //     }


  unsigned int el_type;
 
  //
  //
  //  makes vector with lines of BC elements   
  for (unsigned int i =0; i< elem_values.size();++i)
  {

    el_type = elem_values[i][1];   //  OK  VERSION 2

    switch(dim) {
    case 1:
      // 1D  sim ................
      if  (elem_values[i][1] == 15)   //   BC = 0D  point (1 node)
      {
        BC_elem_line.push_back(elem_values[i]);
      }

      break;
    case 2:
      if ( (elem_values[i][1] == 1) || (elem_values[i][1] == 8) )
      {
        BC_elem_line.push_back(elem_values[i]);
      }

      break;
    case 3:
      if ( (elem_values[i][1] == 2) || (elem_values[i][1] == 3) || 
           (elem_values[i][1] == 9) || (elem_values[i][1] == 10) )
      {
        BC_elem_line.push_back(elem_values[i]);
      }

      break;
    default:
      cout <<  "Error : wrong  value  of  dim"; 

    }

  }


  unsigned int   BC_id;
  BC_id = 0;
  vector<unsigned int> BC_id_vec;
  vector<unsigned int> :: iterator p1;

  // makes  a  vector  of  all different BC  id  from  gmsh  file,  
  // to  compare with user's BC_reg_ID  ->  THIS  CHECK  WILL  BE  DONE  SOMEWHERE ELSE !

  for (unsigned int i =0; i< BC_elem_line.size(); ++i)
  {


    if (version == 1)
    {
      BC_id = BC_elem_line[i][2]  ;
    }

    else if (version == 2)
    {
      BC_id =  BC_elem_line[i][3]  ;
    }
   

    // **** VERSION 2 ***
    //  BC_id =  BC_elem_line[i][3]  ; !!!!!!!!!!!!!
    //

    p1 = find(BC_id_vec.begin(), BC_id_vec.end(), BC_id );  
    // find if physic_id has  been already taken

    if (p1 ==  BC_id_vec.end())   //  not  found
    {
      BC_id_vec.push_back(BC_id);
      //   cout << " BC_id_vec " << BC_id << endl ;
    }

  }

  // cout << " BC_id_vec " << BC_id_vec[0]<< "   " << BC_id_vec[1]<< "   " 
  // <<BC_id_vec[2]<<endl;

 
  // ************************************************************
  //  cross check between BC_id_vec and BC_reg_ID ->  THIS  CHECK  WILL  BE  DONE  SOMEWHERE ELSE !
  //  cross_check_regions( BC_reg_ID,BC_id_vec, "BC");
  //  ***********************************************************

  //  p_2 = BC_reg_ID.begin();

  //  cout<< endl  << " BC  regions ....   "  << endl ;


  //  while (p_2 != BC_reg_ID.end())  //display 
  //    cout<<*p_2++ << endl ;


  // **********************************
  //creates vector of BC_nodes vectors (one for each  BC  region)
  // **********************************

  vector< vector<unsigned int> > BC_nodes;
  vector< vector<unsigned int> > :: iterator p ;


  // initialization of BC_nodes vector
  vector<unsigned int> null_vector;
  //  for (unsigned int i =0; i< BC_reg_ID.size(); ++i)  
  for (unsigned int i =0; i< BC_id_vec.size(); ++i)
  {
    BC_nodes.push_back(null_vector);
  }



  //  int  node1,  node2,  pos, id;
  unsigned int  pos, id, number_of_nodes ;
  unsigned int new_node_id;
  double old_node_id = 0.0;

  unsigned int number_of_tags;  // for  version 2
  number_of_tags = 0;

  vector<unsigned int> node_id;

  id = 0;
  //  node1 = 0;
  //  node2 =0;

  //  vector<int> node_id;

  //  null_vector.push_back(0);

 
 
  // (BC_nodes[0]).push_back(node1);


  for (unsigned int i =0; i< BC_elem_line.size(); ++i)
  {

    if (version == 1)
    {
      id = BC_elem_line[i][2]  ;
    }
    else if (version == 2)
    {
      id = BC_elem_line[i][3]  ;
    }

    // **** VERSION 2 ***
    // id = BC_elem_line[i][3]  ;
    //


    if (version == 1)
    {
      number_of_nodes = BC_elem_line[i][4];  //  number  of  nodes  
      // belonging  to current  element
    }
    else if (version == 2)
    {

      number_of_tags = BC_elem_line[i][2];
      number_of_nodes = (BC_elem_line[i].size() ) - (3 + number_of_tags);

    }

    //**************  VERSION 2
        // number_of_nodes = BC_elem_line[i].size() - (3 + BC_elem_line[i][2])   (see  get_elem_nodes() )


        //   cout<< " number_of_nodes  = " << number_of_nodes<< endl;

        for (unsigned int j =0; j<  number_of_nodes; ++j)
        {

          //  cout<< " (BC_elem_line[i][5+j]  = " << BC_elem_line[i][5+j]<< endl ;


          // ************************************
          // NEW  CHANGE :  13.12.05
          //  convert  to  new  node  numeration :
          // ************************************

          //  num_of_nodes = total  number  of  nodes in the  mesh



          if (version == 1)
          {
            old_node_id = (double)(BC_elem_line[i][5+j]);
          }
          else if (version == 2)
          {
            old_node_id = (double)(BC_elem_line[i][3 + number_of_tags+j]);
          }
          // **************  VERSION 2
          // number_of_tags = BC_elem_line[i][2]
          //  BC_elem_line[i][3 + number_of_tags+j])


          for (unsigned int k =0; k<num_of_nodes ; ++k)
          { 
	     
            if  (node_label[k] ==  old_node_id)
            { new_node_id = k; }

            //  	{ new_node_id = k+1; }
            // +1  -1  -> +0 ,  to  get  nodes  starting  from  zero [ see write_xda ( )]
          }


          //	  node_id.push_back(BC_elem_line[i][5+j]) ;
          node_id.push_back( new_node_id);

          old_node_id = 0;
          new_node_id =0;



          //  node1 = BC_elem_line[i][5]  ;
          // cout<<  node1  << endl ;

          //   node2 = BC_elem_line[i][6]  ;
          // cout<<  node2  << endl ;

        }
   
    //  pos = find_pos(id,BC_reg_ID);
    pos = find_pos(id,BC_id_vec);

    //   if (BC_nodes[pos].size() == 1)
    //        	{
    //        	  BC_nodes[pos].pop_back();  //  remove  initialization null
    //        	}

    //   cout<< " node_id.size()  = " << node_id.size() << endl ;
    for (unsigned int i =0; i< node_id.size(); ++i)
    { 

      //  BC_nodes[pos].push_back(node1);
      //  BC_nodes[pos].push_back(node2);

      //	  cout<< " node_id[i]  = " << node_id[i] << endl ;

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

  //  cout<< " BC_nodes.size()  = " << BC_nodes.size()  << endl ;
  //  cout<< " BC_nodes[0].size()  = " << BC_nodes[0].size()  << endl ;




  //   for (int i =0; i< BC_nodes.size(); ++i)
  //     { 
  //       for (int j =0; j< BC_nodes[i].size(); ++j)
  // 	{

  // 	  cout<< " BC_nodes[i][j]  = " <<BC_nodes[i][j] << endl;
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
    //   for (unsigned int i =0; i< BC_reg_ID.size(); ++i)
    for (unsigned int i =0; i<BC_id_vec.size(); ++i)
    {
      BC_elements.push_back(null_vector);
    }


    for (unsigned int i =0; i< BC_elem_line.size(); ++i)
    {


      if (version == 1)
      {
        id = BC_elem_line[i][2]  ;
      }
      else if (version == 2)
      {
        id = BC_elem_line[i][3]  ;
      }
      // **** VERSION 2 ***
      //  id = BC_elem_line[i][3]  ;
      //

      el_id = BC_elem_line[i][0]  ;

      //   pos = find_pos(id,BC_reg_ID);
      pos = find_pos(id,BC_id_vec);

      BC_elements[pos].push_back(el_id);
      el_id = 0;

    }




    //  creates  map  <BC_region, BC_elements >      //   for   dim = 3

    //  map<int, vector<int> >  // BoundCond_elem   in  .h  !!!

    for (unsigned int i =0; i< BC_elements.size();++i)
    {
      //     BoundCond_elements.insert(make_pair(BC_reg_ID[i], BC_elements[i]) );
      BoundCond_elements.insert(make_pair(BC_id_vec[i], BC_elements[i]) );
      //   cout << BC_reg_ID[i]<< endl ;
      // cout << BC_nodes[i][0]<< endl ;

    }

  }

  //   END  BC_elements (only if   dim=3 )

  // ****************************************



  // ++++++++++++++++++++++++++++++++++++++++
  // FOR  2D  !!!!!!!!!!
  //+++++++++++++++++++++++++++++++++
  //  creates  map  <BC_region, BC_nodes>

  //  map<int, vector<int> >  BoundCond   in  .h  !!!

  for (unsigned int i =0; i< BC_nodes.size();++i)
  {
    //  BoundCond.insert(make_pair(BC_reg_ID[i], BC_nodes[i]) );
    BoundCond.insert(make_pair(BC_id_vec[i], BC_nodes[i]) );
    //    cout << BC_reg_ID[i]<< endl ;
    //   cout << "BC_id_vec[i]" << BC_id_vec[i]<<  endl ;

    //    for (unsigned int j =0; j< BC_nodes[i].size();++j)
    //     {
    //       cout << "BC_nodes[i][j]" << BC_nodes[i][j]<<  endl ;
    //     }


  }

}

//  END  get_BC_info() !!!





//  utility  function 
//unsigned int  Read_MSH::find_pos( unsigned int  reg_id , 
//  vector<unsigned int>& BC_reg_ID )
unsigned int  Read_MSH::find_pos( unsigned int  reg_id ,
                                  vector<unsigned int>& user_reg_ID )

{

  unsigned int  pos;
  unsigned int count = 0;
  bool found = false;


  while (found == false)  {
     
    if  (  (user_reg_ID[count]) != reg_id)
    {
      count++; 
    }
    else 
    {
      found = true;
    }


  }  //  end  while

  if  (found == true)

  {

    pos =  count;
    return  pos;

  }

  else  cerr  <<  "Error :  inconsistent  Physical or  BC  region definition !! " 
              << endl ;
  // physic ID  not  found in  user's  list 


}





// public  method  to  get  map    <BC_region, BC_nodes> 
void  Read_MSH::get_BC_data (map<unsigned int, vector<unsigned int> >& BoundCond_map )
{

  std::map< unsigned int,std::vector<unsigned int> >::iterator it(BoundCond.begin());
  const std::map< unsigned int,std::vector<unsigned int> >::iterator end(BoundCond.end());
  for ( ; it != end; ++it)
  {

    //   cerr << "BC_map ID  =  " << (it->first) <<  endl;
    for ( int i =0; i < ((it->second).size()); ++i)
    {
      //    cerr <<  " node " << i << "  " << (it->second)[i] << endl;
    }

  }


  BoundCond_map = BoundCond;
 
}


// public  method  to  get  map  <BC_region, BC_nodes> and map <BC_region, BC_elements> 
void  Read_MSH::get_BC_data (map<unsigned int, vector<unsigned int> >& 
                             BoundCond_map, map<unsigned int, vector<unsigned int> >&
                             BoundCond_el_map      )
{

  BoundCond_map = BoundCond;

  //   if  dim  =  3  : BoundCond_elem 

  BoundCond_el_map = BoundCond_elements;

 
}


// public  method  to  get  map  <phys_region, elements>  
void  Read_MSH::get_elem_data (map<unsigned int, vector<unsigned int> >& 
                               PhysReg_elements_map)
{ 
      
  PhysReg_elements_map = PhysReg_elements;  
}





// get  physical  regions  for  each  element
//void Read_MSH::get_physical_elem(vector<unsigned int>& phys_reg_ID)
void Read_MSH::get_physical_elem()
{ //get_phys

	
  //  get   elements  associated  to  each  physical  region

  // ******************************************************************
  //  works  only  for  .msh file  format  v.1  !!!
  // *********************************************************
  //  elem_id , el_type, phys_reg, geom_reg, number of  nodes, list of  nodes
 
  unsigned int el_type;

 
  // **********************************************

  //  idea  !:  work  on    vector<Element>  list_elements   instead   
  // of phys_elem_lines
  //   list_elements  has been  already created [in  write_xda() ] with  
  // the  elements with dim  =  sim_dim !!
  // **********************************************

 


  unsigned int id;  // physical reg id
  unsigned int  pos;
  unsigned int el_id=0;
  id = 0;

  vector<unsigned int> null_vector;

  vector< vector<unsigned int> > phys_elements;  //get_physical_elem



  unsigned int  physic_id;

  vector<unsigned int> physic_id_vec;   //  NEW :   use   only  this  vector for  physical regions  ID !!
  physic_id = 0;
  // int id = 0;
  vector<unsigned int> :: iterator p;
  vector<unsigned int> :: iterator p1;
  vector<unsigned int> :: iterator p2;





  // makes  a  vector  of  all different phisic  id  from  gmsh  file,  
  // to  compare with user's phys_reg_ID
  for (unsigned int i =0; i< list_elements.size(); ++i)
  {
    physic_id = list_elements[i].phys_id;

    // IDEA :  put physic_id in  vector  
    p = find(physic_id_vec.begin(), physic_id_vec.end(), physic_id );  
    // find if physic_id has  been already taken

    if (p ==  physic_id_vec.end())   //  not  found
    {
      physic_id_vec.push_back(physic_id);
   	 
    }

  }

  // cout << " physic_id_vec " << physic_id_vec[0]<< "   " << physic_id_vec[1]<< 
  // "   " <<physic_id_vec[2]<<endl;

  // NEW :   ELIMINATE   phys_reg_ID !!!!
  // ---------------------------------------
  // ************************************************************
  //  cross check between physic_id_vec and phys_reg_ID

  //   if  (phys_reg_ID.empty() ) 
  //   {
  //     cerr  <<  " ERROR : No Physical  region has   been    defined !! " <<  endl;
  //     exit(1);
  //   }

  //   cross_check_regions( phys_reg_ID,physic_id_vec, "phys");
  //   THIS  CHECK   WILL  BE   SOMEWHERE  ELSE !




  // **************************************************************
  // creates   vector of  physical elements
  // ****************************************

  // initialization of phys_elements vector
  //  vector<int> null_vector;
  //  for (unsigned int i =0; i< phys_reg_ID.size(); ++i)
  for (unsigned int i =0; i< physic_id_vec.size(); ++i)
  {
    phys_elements.push_back(null_vector);
  }




  // *******************************************************
  //***************************************************
      //  with  vector<Element>  list_elements;

 

      for (unsigned int i =0; i< list_elements.size(); ++i)
    {
      physic_id = list_elements[i].phys_id;


      //    pos = find_pos(physic_id,phys_reg_ID);
      pos = find_pos(physic_id,physic_id_vec);

      // if  not  found ->  error  inconsistency

      el_id = (i +1);
      phys_elements[pos].push_back(el_id);
      //	    
      el_id = 0;
    }


  // **************************************************************


  // *******************************************************
  //************************************************




      //  ??????????????????????????
      //  creates  map  <phys_region, phys_elements >      //  

      //  map<int, vector<int> >  // PhysReg_elem   in  .h  !!!

      for (unsigned int i =0; i< phys_elements.size();++i)
    {
      //    PhysReg_elements.insert(make_pair(phys_reg_ID[i], phys_elements[i]) );

      PhysReg_elements.insert(make_pair(physic_id_vec[i], phys_elements[i]) );


      //   cout << BC_reg_ID[i]<< endl ;
      // cout << BC_nodes[i][0]<< endl ;

    }
	
  //  ??????????????????????????


	



  //   *********************************************************************
  //    map   < element number , phys_reg >  
 	 
  //	map<int, int >   elem_region_map;  in  .h  !!!

  //	cout << endl <<  "phys_elements.size()" << phys_elements.size()  <<endl<< endl ;
	
  for (unsigned int i =0; i< phys_elements.size();++i)
  {
    for (unsigned int j =0; j< phys_elements[i].size();++j)
    {
	  
      //    elem_region_map.insert(make_pair( phys_elements[i][j],phys_reg_ID[i] ) );

      elem_region_map.insert(make_pair( phys_elements[i][j],physic_id_vec[i] ) );



    }
  }
        
        
  // cout << endl <<  "content of  map elem / PhysReg " << endl<< endl ;
       
  //  END  get_physical_elem  () !!!
	
}



//  not  used ???
//  
// public  method  to  get  map  <phys_region, elements>  
void  Read_MSH::get_elem_phys_map (map<unsigned int,unsigned int> &elem_phys  )
{ 
      
  elem_phys =   elem_region_map ;
   
}



// write Libmesh  .xda  file
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

    //                  PYRAMID5,   // 16    BEWARE :   Element NOT correctly handled in  LIBMESH !! (23.11.08)
    //
    //
    // ************************************************************

    //**************************************************************************
      // GMSH  ELEM TYPES
      //**************************************************************************


      // 1 Line (2 nodes).  [ -> EDGE2 = 0] 
      // 2 Triangle (3 nodes).
      // 3 Quadrangle (4 nodes). 
      // 4 Tetrahedron (4 nodes). 
      // 5 Hexahedron (8 nodes). 
      // 6 Prism (6 nodes). 
      // 7 Pyramid (5 nodes).
      //  8 Second order line (3 nodes: 2 associated with 
      // the vertices and 1 with the edge).  [ -> EDGE3 = 1]
      //  
      // 9 Second order triangle (6 nodes: 3 associated with 
      // the vertices and 3 with the edges). 
      //
      // 10 Second order quadrangle (9 nodes: 4 associated with 
      // the vertices, 4 with the edges and 1 with the face).
      //
      //  11 Second order tetrahedron (10 nodes: 4 associated with 
      // the vertices and 6 with the edges). 
      //
      // 12 Second order hexahedron (27 nodes: 8 associated with 
      // the vertices, 12 with the edges, 6 with the faces and 1 with the volume). 
      //
      // 13 Second order prism (18 nodes: 6 associated with the vertices, 
      // 9 with the edges and 3 with the quadrangular faces).
      //
      //  14 Second order pyramid (14 nodes: 5 associated with the vertices, 
      // 8 with the edges and 1 with the quadrangular face). 
      //
      // 15 Point (1 node).  [no  correspondance in  Libmesh !]




	
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
        // 1D  sim ................  line EDGE2 = 0 (2 nodes)
        elem_type.push_back(0);

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
        // BEWARE !   PYRAMID  element  not  correctly  handled in  LIBMESH
        //
        throw InitFailedException("ERROR: PYRAMID5 element not implemented...");
        //   elem_type.push_back(16);
        break;

      case 8:
        // 1D   second order line (3 nodes) ........ EDGE3 = 1
        //  elem_type =1 ; 
        elem_type.push_back(1);
     
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
    //90	 # Sum of Element Weights= sum of nodes of  all the  
    //             elements = Sigma(i=Nel)[ nodes(i)]
    //0	 # Num. Boundary Conds.
    //65536	 # String Size (ignore)
    //1	 # Num. Element Blocks. = "number  of  mesh  blocks" =  1 !
    //3 	 # Element types in each block. = vector of  element types 
    //             (TRI = 3) (see  elem_type.h)
    //30 	 # Num. of elements in each block = total number of elements 
    //            of a given type
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
          //    cout << "All_elements[j].type " << All_elements[j].type << endl;

          list_elements.push_back(All_elements[j]);
          count++;

        }

      }

      num_elem_per_type.push_back(count);
      count = 0;


    }






    // ************************************************************







  
    // cout << " num_of_nodes" <<  num_of_nodes<< endl;
    


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

        node = list_elements[i].nodes[j] - 1; //  -1  ->   .xda  starts  from  zero
        //   out << ( list_elements[i].nodes[j] - 1)   << "  " ;
        out << node  << "  " ;

      }

      out << "\n";

    }

    // ******************************



    //    OLD
    //	
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

    for (unsigned  int node_id =0; node_id<num_of_nodes ; node_id++)
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


  // END  write  .xda  file

   
}

//





// get  element data  for  each  element
// BEWARE ! we  select ONLY  the  elements relevant to the  
// simulation dimension (dim):
// e.g. if  dim = 2 ->  here we take  only  2D elements,  
// which correspond to  phys. regions;
//   1D elements,  if  present, are processed  separately to  get bound. cond. reg.
//

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

    el_type = elem_values[i][1];   // OK FOR VERSION 2

    
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
    // 	  // 	  //  cout<< BC_reg_ID[0]  << endl ;

    // 	}


    //.......................


    switch(dim) {
    case 1:
      // 1D  sim ................elem_values[i][1] == 1) || (elem_values[i][1] == 8)
      if ( (el_type == 1) || (el_type == 8)  )
      {
        elem_line.push_back(elem_values[i]);
        p = find(gmsh_elem_type.begin(), gmsh_elem_type.end(),el_type  );

        if (p == gmsh_elem_type.end())   //  not  found
        {
          gmsh_elem_type.push_back(el_type);

          //	cout<< "el_type*************"<<  el_type   << endl ;

        }


      }





      break;
	
    case 2:
      if ( (el_type == 2) || (el_type == 3) || (el_type == 9) || (el_type == 10) )
      {
        elem_line.push_back(elem_values[i]);
        p = find(gmsh_elem_type.begin(), gmsh_elem_type.end(),el_type  );

        if (p == gmsh_elem_type.end())   //  not  found
        {
          gmsh_elem_type.push_back(el_type);

          //	cout<< "el_type*************"<<  el_type   << endl ;

        }


      }

      break;
    case 3:
      if ( (el_type == 4) || (el_type == 5) || (el_type == 6) || 
           (el_type == 7)|| (el_type == 11)|| 
           (el_type == 12)|| (el_type == 13) || (el_type == 14))
      {
        elem_line.push_back(elem_values[i]);
        p = find(gmsh_elem_type.begin(), gmsh_elem_type.end(),el_type  );

        if (p == gmsh_elem_type.end())   //  not  found
        {
          gmsh_elem_type.push_back(el_type);

          //	cout<< "el_type*************"<<  el_type   << endl ;

        }



      }

      break;
	
    default:
      cout <<  "Error :    wrong  value  of  dim"; 

    }

    //   if (p == gmsh_elem_type.end())   //  not  found
    //  	{
    //  	  gmsh_elem_type.push_back(el_type);

    // 	  cout<< "el_type*************"<<  el_type   << endl ;

    // 	}



  }

  //
  //NOW IN elem_line THERE ARE ONLY THE ELEMENT ENTRIES CONSISTENT WITH DIM OF PROBLEM !!!
  //


  unsigned int new_node_id;
  double node_id =0.0 ; 
  unsigned int count = 0;
  unsigned int  number_of_nodes_in_elem = 0;
  vector<unsigned int> node_id_list;

  unsigned int number_of_tags;
  number_of_tags = 0;

  //  int count_ok = 0;
  //   int count_swap= 0;
  //  count_ok = 0;
  //  count_swap= 0;

  //  unsigned int temp,type ;
  unsigned int type;

   

  

  for (unsigned int i =0; i< elem_line.size(); ++i)
  {  // list of lines


   

    if (version == 1)
    {
      number_of_nodes_in_elem = elem_line[i][4];
    }
    else if (version == 2)
    {
      number_of_tags = elem_line[i][2];
      number_of_nodes_in_elem = ( elem_line[i].size() ) - (3 + number_of_tags);

    }

    //**************  VERSION 2
        //elem_line[i].size() -  (1 + 1+ + 1 + number_of_tags),  that is
        //    number_of_nodes_in_elem = elem_line[i].size() - (3 + elem_line[i][2])

                  
        type =  elem_line[i][1];   //  ok VERSION 2 




    for (unsigned int j =0; j<  number_of_nodes_in_elem; ++j)
    {



      if (version == 1)
      {
        node_id  = (double)(elem_line[i][5+j]);// ?????? TO DO: change  node_label 
        //                                      and  node_entry to  int !!!
      }
      else if (version == 2)
      {      
        node_id  = (double)(elem_line[i][3 + number_of_tags+j]);
      }

      // **************  VERSION 2
      //  elem_line[i][3 + number_of_tags+j])
      //

      for (unsigned int k =0; k<num_of_nodes ; ++k)
      { 
        //  cout << " node_label[k]" <<  node_label[k]<< endl;
        if  (node_label[k] ==  node_id)
        { new_node_id = k+1; }  //   +1  ->  starts from  1
      }
      node_id_list.push_back( new_node_id) ;
      //	 cout << "new_node_id " << new_node_id<< endl ;

      //	 node_id_list.push_back(elem_line[i][5+j]) ;

      count ++;
      new_node_id = 0;
      node_id  =0;

    }

    //*********************
        //  NEW   10.1.06

        //  HERE CHECK  ORIENTATION OF  ELEMENT :

        // for each element: calculate determinant delta on node = (node_id_list[i]-1)
        // coord are given  by node_coord[node][j]
        // if  delta > 0 -> ok,   else swap nodes in  node_id_list

        //	  cout << "**********"<< endl;
        //   cout << node_id_list[0]-1 << "  " << node_id_list[1]-1 << 
        // "  " << node_id_list[2]-1 << endl ;
        //   cout <<   node_coord[1][0]<< endl ;

        // **********************************************
        check_orientation(node_id_list,type);
        
    //  TEST !!!!!!!!!!!
    //  cout <<  " Second  time  det  should  be  > 0 !!!   "  << endl;
    //         check_orientation(node_id_list,type);
        
    //*****************************************************

        //    cout << "+++++++++++++++++++    CHECK orient +++++++++++++" << endl;


        //    if ( check_orientation(node_id_list) ) 
        // 	{
        // 	  cout << " check_orientation OK  !!!  " << endl;
        // 	  count_ok++;

        // 	}
        //       else 
        // 	{
        // 	  cout  << " SWAP NODES  !!!  " << endl;
        // 	  count_swap++;

        // 	  if (dim == 2)
        // 	    {

        // 	      // swap node_id_list[0] and  node_id_list[2]
        // 	      temp = node_id_list[0];
        // 	      node_id_list[0] = node_id_list[2];
        // 	      node_id_list[2] = temp;
        // 	    }
        // 	  else if (dim == 3)
        // 	    {

        // 	      if (type == 6)   // gmsh  PRISM6

        // 		{
        // 		  // swap
        // 		}
        // 	  //  




        // 	    }



        // 	}


  




        // *******************************************************


        current_element.nodes = node_id_list;
    current_element.type =  elem_line[i][1]; //  ok VERSION 2 

    if (version == 1)
    {
      current_element.phys_id =  elem_line[i][2];
    }
    else if (version == 2)
    {
      current_element.phys_id = elem_line[i][3];
    }

    // **** VERSION 2 ***
    // current_element.phys_id = elem_line[i][3]; !!!!!!!!!!!!!!!


    All_elements.push_back(current_element);
    // cout << "******************" << current_element.type<< endl ;
    current_element.nodes.clear();
    current_element.type = 0;


    //****************



        elem_nodes.push_back(node_id_list);

    node_id_list.clear();

  }

  // end list of  (elem)lines



  el_weight = count ;
  num_of_elem = elem_nodes.size();  //  NECESSARY FOR .XDA FILE (num_of_elem = NUMBER OF PROPER ELEMENTS 
  //  TO  BE  PUT  IN  LIBMESH  FILE,  (DIFFERENT FROM THE  TOTAL NUMB. OF  GMSH ELEM.  !!!!)

  //   cout << "   count_ok = " << count_ok << endl;
  //   cout << "   count_swap = " <<  count_swap << endl;

}  //  end  get_elem_nodes










// ************************* NEW
//  parser of  NODES  section of  .msh  file  (version 1)
void Read_MSH::parse_node_section(ifstream& in_stream )
{ // parse nodes

  //  bool end_nodes = false;
 
  double node_id, x_coord, y_coord, z_coord; 

  string  str, end_nodes_label   ;
  vector<double> v;
  v.clear();

 
  in_stream >> num_of_nodes;
 
  for (int i =0; i< num_of_nodes;++i)
  {

    v.clear();

    in_stream >>  node_id;
    v.push_back( node_id);

    in_stream >>  x_coord;
    v.push_back( x_coord);

    in_stream >>  y_coord ;
    v.push_back(y_coord );

    in_stream >>  z_coord;
    v.push_back(z_coord );

    node_entry.push_back(v);


  }


  in_stream >> end_nodes_label;

}


// ************************* END NEW ***********************




//  check if  $NOD  header  is  found 
void Read_MSH::find_node_section(char const* str, ifstream& in_stream)
{ //
         
  string name;
    
  if  (parse(str, if_p("$")[(+alpha_p)[assign_a(name)]].else_p[nothing_p],
             space_p).full)
  {
    //    cout << name<< endl;
    if (name == "NOD")

    {
      parse_node_section(in_stream);
    }
         
         
  }
}                                                                     




// OK FOR  VERSION 2!!
// get  node  data  for  each  node
void Read_MSH::get_nodes_coord()
{ //get_nodes_coord

  vector<double> temp;

  node_label.clear();
  for (unsigned int i =0; i<num_of_nodes ; ++i)
  {  

    //     cout << "(node_entry[i][0] " <<node_entry[i][0]<< endl ;
    node_label.push_back(node_entry[i][0]) ; // ????


    //   cout << " node_label[i]" <<  node_label[i];

    for (unsigned int j =0; j<3 ; ++j)
    { 

      temp.push_back(node_entry[i][j+1]);
    }

    //  node_coord[i].push_back(temp);
    node_coord.push_back(temp);


    temp.clear();


  }


} // end get_nodes_coord






// void Read_MSH::read_data_section(char const* str,ifstream& in_stream )
// // void InputParser::read_data_section(char const* str)


// // void read_data_section(char const* str,ifstream& in_stream )

// {

//   string name;

//   if  (parse(str, if_p("$")[(+alpha_p)[assign_a(name)]].else_p[nothing_p],
//              space_p).full)
//   {
//     //  cout << name<< endl;
//     //    if (name == "ELM")


//     if ( (name == "ELM") || (name == "Elements") )
//       //	if (name == name_data )
//     {
//       parse_elem_section(in_stream);

//     }

//     //   else if  (name == "NOD")
//     else if  ( (name == "NOD") ||  (name == "Nodes") )
//     {
//       //  parse_2(in_stream);
//       parse_node_section(in_stream);

//     }

//     // **********************  version 2
//     else if (name == "MeshFormat")
//     {
//       parse_meshformat_section(in_stream);

//     }

//     //*********************



//         }



// }



void Read_MSH::parse_meshformat_section(ifstream&  in_stream)
{

  std::string  End_mesh_format  ;
  int file_type,data_size;
  double version_number;

  in_stream >>  version_number ;  // int(version_number ) = 2

  version = int(version_number );  //  global variable private  member

  if (version != 2)
    throw InitFailedException("ERROR: this GMSH mesh file version is not yet implemented"); 


  in_stream >>  file_type ;
  in_stream >>  data_size ;

  in_stream >>  End_mesh_format ;

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

 

  //  num_of_elem is  a  member data
  
  out << " # Data description \n";
  out << "REAL	 # type of values     \n";
  out << "0	 # No. of nodes for which data is stored     \n";
  out << num_of_elem  << "	 # No. of elements for which data is stored    \n";

  for (unsigned int el_id =0; el_id<num_of_elem; el_id++)
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

  string  mesh_file_data = "elem_data.xta";
  string  mesh_file_inp = "in.xda";

  mesh.read (mesh_file_inp,&mesh_data  ); 
  
  mesh_data.read(mesh_file_data);

}




//  check if orientation  of  nodes  is  positive,  otherwise swap  nodes  of  element
void   Read_MSH::check_orientation( vector<unsigned int>&  node_id_list , unsigned  int type )

{


  vector<unsigned int> temp_nodes_list;
  bool swap;
  unsigned int temp  ;
  unsigned int temp_list[3];

  swap = false;


  double det  ;    //    determinant
  double x0,x1,z0, y0,y1,z1,  x2,y2, z2, x3, y3, z3 ,x4,y4,z4, 
    a11,a21, a12, a22, a31, a32, a13, a23, a33;

  unsigned int  size_node_list ;

  size_node_list =  node_id_list.size();

  for (unsigned int i =0; i< size_node_list ; ++i)
  {
    temp_nodes_list.push_back(node_id_list[i]-1 );
  }




  // OK  FOR  FIRST  ORDER ELEMENTS
  //   SECOND ORDER ...............TO  DO  !!?





 
  if  (dim == 2)
    //  ok for  TRI4  AND  QUAD4 

  {

    x0 = node_coord[( temp_nodes_list[0]) ][0];
    x1 = node_coord[( temp_nodes_list[1]) ][0];



    y0 = node_coord[( temp_nodes_list[0]) ][1];
    y1 = node_coord[( temp_nodes_list[1]) ][1];

    x2 = node_coord[( temp_nodes_list[2]) ][0];
    y2  = node_coord[( temp_nodes_list[2]) ][1];

    a11 = x1-x0;
    a21 =   y1-y0;

    a12 = x2 - x0;
    a22 = y2-y0;


    det = a11 * a22  - (a21*a12) ;

    assert(abs(det)> 1e-12);

    if  (det > 0.0)
    {swap = false;}
    else swap = true;

    if (swap)   // if  swap = true  (det <= 0 )   then  reorder nodes of  element
    {
      
      // swap node_id_list[0] and  node_id_list[2]
      temp = node_id_list[0];
      node_id_list[0] = node_id_list[2];
      node_id_list[2] = temp;
      
      // cout << " change orientation  !!!  " << endl;
      //  cout << " det =     " << det << endl;
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


  if (dim == 3)

  { // dim=3


    x0 = node_coord[( temp_nodes_list[0]) ][0];
    y0 = node_coord[( temp_nodes_list[0]) ][1];
    z0 =  node_coord[( temp_nodes_list[0]) ][2];


    x1 = node_coord[( temp_nodes_list[1]) ][0];     
    y1 = node_coord[( temp_nodes_list[1]) ][1];
    z1 = node_coord[( temp_nodes_list[1]) ][2];

    x2 = node_coord[( temp_nodes_list[2]) ][0];
    y2  = node_coord[( temp_nodes_list[2]) ][1];
    z2 = node_coord[( temp_nodes_list[2]) ][2];

    x3 =  node_coord[( temp_nodes_list[3]) ][0];
    y3 =  node_coord[( temp_nodes_list[3]) ][1];
    z3 = node_coord[( temp_nodes_list[3]) ][2];

    //     x4 =  node_coord[( temp_nodes_list[4]) ][0];
    //       y4 =  node_coord[( temp_nodes_list[4]) ][1];
    //       z4 =  node_coord[( temp_nodes_list[4]) ][2];


    if (type == 5)  //  HEX8
      //     if ( (type == 5) || (type == 7) ) //  HEX8 or  PYR5

    {

      x4 =  node_coord[( temp_nodes_list[4]) ][0];
      y4 =  node_coord[( temp_nodes_list[4]) ][1];
      z4 =  node_coord[( temp_nodes_list[4]) ][2];

 
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

    else if (type == 7)  //  PYR5

    {

      x4 =  node_coord[( temp_nodes_list[4]) ][0];
      y4 =  node_coord[( temp_nodes_list[4]) ][1];
      z4 =  node_coord[( temp_nodes_list[4]) ][2];

 
      a11 = x1-x0;
      a21 = y1-y0;
      a31 = z1-z0;

      a12 = x4-x0;
      a22 = y4 - y0;
      a32 = z4- z0;

      a13 =x3 - x0;
      a23 = y3 - y0;
      a33 = z3 - z0;

    }




    else  //    TET4 , PRISM6

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

      switch(type) {
        //   switch(gmsh_elem_type) 
      case 4: //  TET4

        // swap node_id_list[1] and  node_id_list[3]
        temp = 0;
        temp = node_id_list[1];
        node_id_list[1] = node_id_list[3];
        node_id_list[3] = temp;



        break;

      case 6: // gmsh  PRISM6
        // swap   (0 1 2) -> (3 4 5 )


        for (unsigned int i=0; i<3; i++)
        {

          temp_list[i] = node_id_list[i];

        }

        for (unsigned int i=0; i<3; i++)

        {
          node_id_list[i] =node_id_list[i+3] ;
        }

        for (unsigned int i=0; i<3; i++)

        {
          node_id_list[i+3] = temp_list[i]  ;

        }

        break;

      case 5: // gmsh  HEX8
        // swap   (0 1 2 3) -> ( 4 5  6  7  )

        for (unsigned int i=0; i<4; i++)
        {

          temp_list[i] = node_id_list[i];

        }

        for (unsigned int i=0; i<4; i++)

        {
          node_id_list[i] =node_id_list[i+4] ;
        }

        for (unsigned int i=0; i<4; i++)

        {
          node_id_list[i+4] = temp_list[i]  ;

        }


      case 7: //  PYR5

        // swap node_id_list[1] and  node_id_list[3]
        temp = 0;
        temp = node_id_list[1];
        node_id_list[1] = node_id_list[3];
        node_id_list[3] = temp;



        break;

	   
        // 

      default:
        cout <<  "Error : gmsh_elem_type NOT YET IMPLEMENTED "; 

      }

    }


  }


  // return swap;






}



void   Read_MSH:: cross_check_regions( vector<unsigned int>&  user_reg_list,
                                       vector<unsigned int>&  gmsh_reg_list, 
                                       string type )

{
  unsigned int id;
  vector<unsigned int> :: iterator p1,p2;

  string error_text1_phys,  error_text2_phys,  error_text1_BC,  error_text2_BC;

  error_text1_phys = "ERROR : inconsistent Physical  regions. Physical region in  user's defined list not  found  in  GMSH  file";
  error_text2_phys = " ERROR :  inconsistent Physical  regions.  GMSH region not  found in user's defined physical region list";
  error_text1_BC = "ERROR : inconsistent BC  regions. BC region in  user's defined list not  found  in  GMSH  file";

  error_text2_BC = " ERROR :  inconsistent BC  regions.  GMSH region not  found in user's defined BC region list";

  for (int i =0; i< user_reg_list.size();++i)
  {

    id = user_reg_list[i];
    //  cout << " phys_reg_ID =  " << id << endl;
    p1 = find(gmsh_reg_list.begin(), gmsh_reg_list.end(), id );  
    // find  phys_reg_ID (user's id) in id list from file

    if (p1 == gmsh_reg_list.end())   //  not  found

    {

      if (type == "phys")
      {
        //   cerr  <<  " ERROR2 : inconsistent Physical  regions. 
        // Physical region in  user's defined list not  found  in  GMSH  file   " 
        // <<  endl;
        cerr << error_text1_phys  <<  endl;
        exit (1);
      }
      else if (type == "BC")
      {
        cerr << error_text1_BC  <<  endl;
        exit (1);

      }

    }

  }


  id = 0;
  for (int i =0; i< gmsh_reg_list.size();++i)
  {

    id = gmsh_reg_list[i];
    p2 = find(user_reg_list.begin(), user_reg_list.end(), id );  
    // find  phys_id from GMSH file in  (user's id list)

    if (p2 == user_reg_list.end())   //  not  found

    {

      if (type == "phys")
      {

        // cerr  <<  " ERROR2 :  inconsistent Physical  regions.  
        // GMSH region not  found in user's defined physical region list " <<  endl;
        cerr <<  error_text2_phys <<  endl;
        exit(1);
      }
      else if (type == "BC")
      {
        cerr << error_text2_BC  <<  endl;
        exit (1);

      }

    }

  }

  // **************************************************************


}



//  to handle CR DOS files line termination  
void Read_MSH::cut_off_CR(string& label,  ifstream& in_stream)
{

  string::size_type loc = label.find_first_of( "CR", 0 );
  if( loc != string::npos )
  {
    //    cout << "Found CR at " << loc << endl;
    label.erase(loc);

    in_stream.ignore(256,'\n');  //  if  the   read keyword is CR or  begins with CR : ignore  all  the  line

  }
  

} //  end  method



// ****************************************
