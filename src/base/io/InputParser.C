/*=============================================================================
  Copyright (c) 2002-2003 Joel de Guzman
  http://spirit.sourceforge.net/

  Use, modification and distribution is subject to the Boost Software
  License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
  http://www.boost.org/LICENSE_1_0.txt)
  =============================================================================*/

#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/dynamic.hpp>
#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include "InputParser.h"



///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace boost::spirit;



//InputParser::InputParser(string filename, string section_n)
InputParser::InputParser(string& input_file_name)

{

  //  section_name =  section_n;
  filename = input_file_name;

  initialize_vectors(); // 
  //  reg_values.clear();


  //  find_section(filename);

 

  //	cout << device_regions[0].reg_name << endl;
  //  get_device_data();

 

}

InputParser::~InputParser()
{
}


void InputParser::initialize_vectors()
{
  command_values.clear();
  reg_values.clear();
  reg_id.clear();
  mat.clear();
 

}


// *********************************
// OPTIONS (COMMAND)  SECTION
// *********************************

//void InputParser::parse_command(ifstream& in_stream )
void InputParser::parse_options(ifstream& in_stream )



{

  // *******************************************************************************************************
  //  ***** IDEA:  label = value  ->  map <string,double>  o map <string,string>  , map <string, string> 
  // variable = inp_file.read("prop", int Default)
  //   int  inp_file::read ( string,  int  Default);
  //    find string   in  map  <string,prop_value>  ->  return  value 
  //  name_property  ->  if  find in  vector< name_property>  then  property  = true
  // *******************************************************************************************************



 
  string name, str, mat_name ;
  name = "Start";
  int id;

  // ****************************
  //  assigned   values (local)
  // ****************************
  vector<double> v;              //  scalar  and   vector  values
  vector<int> v_int; 
  //  vector<double>  vect;
  vector<string>  v_string;
  v_string.clear();
  //  vect.clear();
  v.clear();
  v_int.clear();
  //******************


      //------------------------
      //  labels of  properties
      //------------------------

      vector<string> v_label;
  vector<string>vect_label;
  vector<string> v_label_string;  //  label  for  string  value
  v_label_string.clear();
  vect_label.clear();
  v_label.clear();

  //  vect.push_back(0);

   
  // *************************************************************************************************
  //  rule<phrase_scanner_t>  for  phrase level parsing (e.g. with a  separator character (comma) )
  // rule <> for  space  separated  list (character level parsing) : NEEDS  TO READ  EXPLICITLY SPACES !!!
  // *************************************************************************************************




  rule<>list_of_numbers_space_sep = ch_p('(')>> *(space_p) >>*real_p[push_back_a(v)] >>
    *( *(space_p) >> real_p[push_back_a(v)])>> *(space_p)>> ch_p(')')  ;

  // rule<>list_of_numbers_space_sep = ch_p('(')>> *(space_p) >>real_p[push_back_a(v)] >>
  //  *( *(separator) >> real_p[push_back_a(v)])>> *(separator)>> ch_p(')')  ;


  rule<>special_char =  (ch_p('_') | ch_p('-') |  ch_p('.') |  ch_p('/')   |  ch_p('+')   );

  // rule<>label  = (+alnum_p)>>   * ( ch_p('_') >> *(+alnum_p) ) ;
  rule<>label  = (+alnum_p)>>   * ( (special_char ) >> *(+alnum_p) ) ;


  // *** one  can  use $ to distinguish  string_label from label of  numerical prop,
  //  otherwise  one  line  can  contain only  numeric or only string  values !!!
  //  
  //    rule<>label_string  =  (ch_p('$'))>> (+alnum_p)>>   * ( ch_p('_') >> *(+alnum_p) ) ; 
  //   //   (+alnum_p)>> ;//  * ( (special_char ) >> *(+alnum_p) ) ;





   
  // rule<>assignement  =  (+alnum_p)[push_back_a(v_label)] >> *(space_p) >> ch_p('=')>>
  //  *(space_p) >> ( real_p[push_back_a(v)])  ; //ok !!
  // strict_real_p  for real with . !!!!
  rule<>assignement  =  (label)[push_back_a(v_label)] >> *(space_p) >> ch_p('=')>>
    *(space_p) >> ( real_p[push_back_a(v)])  ; // with _ !!

  // for  int
  //  rule<>assignement_int  =  (label)[push_back_a(v_label)] >> *(space_p) >> ch_p('=')>>
  //   *(space_p) >> ( int_p[push_back_a(v_int)])  ;



  //   rule<>assignement_string  =  (label_string)[push_back_a(v_label_bool)] >> *(space_p) >>
  //    ch_p('=')>> *(space_p) >> ((label )[push_back_a(v_bool)])   ; //  ok!!!

  rule<>assignement_string  =  (label)[push_back_a(v_label_string)] >> 
    *(space_p) >> ch_p('=')>> *(space_p) >> ((label)[push_back_a(v_string)])   ;


  rule<>assignement_vector =  (label)[push_back_a(vect_label)] >> *(space_p) >> 
    ch_p('=')>> *(space_p) >> list_of_numbers_space_sep;

 
  //  rule<>assignement = (assignement_double |  assignement_int); 


  //  rule<>list_of_assignement = (assignement | assignement_string)  >> 
  //   *( *(space_p) >> (assignement | assignement_string) ); //  ok


  // ******************************************************************* 
  //list of  only numeric  OR  list  of only  strings !!
  // 
  rule<>list_of_assignement = (assignement )  >>
    *( *(space_p) >> (assignement ) )   | (  assignement_string) >> 
    *( *(space_p) >> (assignement_string ) );

  // ************************************************************************

 
  // rule<> r_command  = *(space_p) >> (list_of_assignement | assignement_vector)    >> *(anychar_p);

  rule<> r_command  = *(space_p) >> (list_of_assignement )    >>   *(space_p) >> !(comment_p("#") >> *(anychar_p)) >>*(space_p)  ;
  // *(anychar_p);  

  //  possibly COMMENTs after list_of_assignement !


  rule<> r_command_vector  = *(space_p) >> (assignement_vector)    >> *(anychar_p);





  while ( getline(in_stream, str) )
  {

   
    if  (!(parse(str.c_str(), comment_p("#")   , space_p).full) ) 

    { // if !  comment_p("#")  

	 

      if (  parse(str.c_str(),

                  //  Begin  grammar
		      
                  r_command 

                  )

            //  ,
            //  End grammar

            //  space_p).full )
            .full )           //  not skipping spaces


      { // if parse

	     

        if ( !(v.empty()) )

        {

          if ( !(v_int.empty()) )
          {
            cout << "v_int[0-------------- "<< v_int[0]<< endl;
          }

          for (int i =0; i< v.size();++i)
          {
            prop_labels.push_back(v_label[i]);
		      
            prop_labels_map.insert(make_pair(v_label[i], v[i]) );
          }
        }


        if ( !(v_string.empty()) )
        {

          for (int i =0; i< v_label_string.size();++i)
          {
            string_prop_labels.push_back(v_label_string[i]);
            string_prop_labels_map.insert(make_pair(v_label_string[i], v_string[i]) );
          }
        }
	     


	      

        v_label_string.clear();
        vect_label.clear();
        v_label.clear();
        v_string.clear();
        v.clear();
        v_int.clear();


	    

      }

      else if(  parse(str.c_str(),  r_command_vector   )   .full ) //  not skipping spaces
        //  reads  list of  values (vector)

      {

        if ( !(v.empty()) )
        {
          vector_prop_labels.push_back(vect_label[0]);
		  
          vector_prop_labels_map.insert(make_pair(vect_label[0], v) );

        }

	     

	      
        v_label_string.clear();
        vect_label.clear();
        v_label.clear();
        v_string.clear();
        v.clear();

      }


      //	  else  if (parse(str.c_str(), if_p("$")[(+alpha_p)[assign_a(name)]]  , space_p).full)
      else  if (parse(str.c_str(), if_p("$")[(+alpha_p)[assign_a(name)] >>
                                             ch_p("$")  ] , space_p ).full)


      {
        if (name == "End")

        {  
          cout << name  << endl ;
          break;
        }


      }



      else  

      {
        cerr <<  "  SYNTAX ERROR in input  file   " <<  endl;
        cerr << " Correct syntax is : 'label' = 'value' 'label' = 'value' .......# 'comment' " 
             << endl;
        cerr << " A comment line  must be preceded by '#' "<< endl;
        cerr << " BEWARE: numerical and string values cannot  be  mixed  in  the  same line !"
             << endl;
        exit(1); 
      }


    }


    v_label_string.clear();
    vect_label.clear();
    v_label.clear();
    v_string.clear();
    v.clear();
    v_int.clear();
     

    //    if (name == "End")
    // 	{  
    // 	  cout << name  << endl ;
    // 	  break;
    // 	}



    //  OK  ***************
    //   else
    // 	{
    // 	  getline(in_stream, str);
    // 	}
    //  OK  ***************



      

  }  //  end  while


}


// *********************************
// END  OPTIONS(COMMAND)  SECTION
// *********************************


// find section = section name (from input) and  parses  just  that  section

//void InputParser::find_section(string file_name)
void InputParser::read_section(string& section_name)

{

  bool found = false;
  std::string str;

  //  std::ifstream in_stream (file_name.c_str());
  std::ifstream in_stream (filename.c_str());


  string name;

  if ( !in_stream.good() )
  {
    std::cerr << "ERROR: Input file not good." 
              << std::endl;
    //   error();
  }



  while (getline(in_stream, str))  //   
  { //while

    //  if  (parse(str.c_str(), if_p("$")[(+alpha_p)[assign_a(name)]>> ch_p("$")].else_p[nothing_p] ).full)
    if  (parse(str.c_str(), if_p("$")[(+(~ch_p('$')))[assign_a(name)]>>
                                      ch_p("$")].else_p[nothing_p] ).full)


    {
      cout << name<< endl;

      //	  if ((name == section_name) && (name == "Regions"))
      if  (name == section_name) 
      {
        found = true;

        if  (section_name == "Regions")

          //	  if (  (section_name == "Regions") && (name == section_name))

        {
          parse_device(in_stream);
        }
        else if  (section_name == "BC_Regions")

        {
          parse_device_BC(in_stream);
        }
    
        else if (section_name == "Alloy")
        
        {
          parse_alloy(in_stream);
        }
    
   
        else 
          //if (name == section_name)
      
          //   {parse_command(in_stream);}
        {
          parse_options(in_stream);
        }

      }


    }

  }

  if (found == false)

  {
    cerr << " Error: Section not  found  " << endl ;
    exit(1);

  }

  in_stream.close();



}





void  InputParser::get_data ( vector< vector<double> >& glob_reg_values,
			      vector< vector<double> >& glob_comm_values,
			      vector<int>& glob_reg_id,  vector<string>& glob_mat)
{
  glob_reg_values = reg_values;
  glob_reg_id = reg_id;
  glob_mat = mat;

  glob_comm_values = command_values;



}




double InputParser::read_input( string label , double  Default)
{
  double value;
  map <string,double>  :: iterator  p;


  //  for (int i =0; i< prop_labels.size();++i)
  //   {
  p = prop_labels_map.find( label  );

  if  (p != prop_labels_map.end() )

  {  
    value  =  (p -> second) ;
	

     
    return value ;

  }

  else 
  { 
    cout  <<  "*** Default ***  ";

    return Default;
    //  cout  <<  "error"  ;
  }
  //  }

}




int  InputParser::read_input( string label , int  Default)
{
  int value;
  map <string,double>  :: iterator  p;


  //  for (int i =0; i< prop_labels.size();++i)
  //    {
  p = prop_labels_map.find( label  );

  if  (p != prop_labels_map.end() )

  {  
    value  =  (int)(p -> second) ;
	

	 
    return value ;

  }

  else
  {
    cout  <<  "*** Default ***  ";
    return Default;
    //	  cout  <<  "error"  ;
  }
  //   }

}




string   InputParser::read_input( string label , string  Default)
{
  string  value;
  map <string,string>  :: iterator  p;


  //  for (int i =0; i< prop_labels.size();++i)
  //    {
  p = string_prop_labels_map.find( label  );

  if  (p != string_prop_labels_map.end() )

  {  
    value  =  (p -> second) ;
	

	 
    return value ;

  }

  else
  {
    cout  <<  "*** Default ***  ";
    return Default;
    //  cout  <<  "error"  ;
  }
  //   }

}



////*****  bool  return  :  overload  not  possible :  why  ???
//bool   InputParser::read_input( string label_bool , bool  default_bool)
bool   InputParser::read_input( string label_string )

{
  string  value;
  map <string,string>  :: iterator  p;
  bool flag;
  flag = 0;//default_bool ;




  //  for (int i =0; i< prop_labels.size();++i)
  //    {
  p = string_prop_labels_map.find( label_string  );

  if  (p != string_prop_labels_map.end() )

  {  
    value  =  (p -> second) ;
	 
    //	if ( strncmp (str[n],"C3**",2) == 0)


    // s.compare(s1);

    if  (  (value.compare("true")) == 0 )  

    { 
      flag = 1 ;//true;
      return flag;

    }
    else if  (  (value.compare("false")) == 0 ) 
    { 
      flag = 0;// false;
      return flag;

    }


	 
    //     return value ;

  }

  else
  {
    cout  <<  "*** Default ***  ";
    return flag ;
    //  cout  <<  "error"  ;
  }
  //   }

}



//  NEW  *********************
// overload  for   vector 

vector<double> InputParser::read_input( string label, vector<double> def_vector)
{

  map <string, vector<double> >   :: iterator  p;

  vector<double>  return_vector;

  p = vector_prop_labels_map.find( label  );

  if  (p != vector_prop_labels_map.end() )

  {  
    return_vector    =  (p -> second) ;
    return  return_vector;

  }

  else
    cout  <<  "Warning --- read_input_vector double: vector data empty or  missing  " << endl ;


}


vector<unsigned int> InputParser::read_input( string label, vector<unsigned int> def_vector)
{

  map <string, vector<double> >   :: iterator  p;

  vector<unsigned int>  return_vector;
  vector<double> temp;
  return_vector.clear();


  p = vector_prop_labels_map.find( label  );

  if  (p != vector_prop_labels_map.end() )

  {  

    temp    =  (p -> second) ;


    for (int i =0; i< temp.size();++i)
    {
      //    cout <<  "temp.size()  " <<  temp.size();
      //     return_vector[i]    = (int)(temp[i]); 
      return_vector.push_back( (int)(temp[i]) );

    }

    return  return_vector;
	
  }

  else
    cout  <<  "Warning --- read_input_vector int: vector data empty or  missing " << endl ;


}




vector<int> InputParser::read_input( string label, vector<int> def_vector)
{

  map <string, vector<double> >   :: iterator  p;

  vector<int>  return_vector;
  vector<double> temp;
  return_vector.clear();


  p = vector_prop_labels_map.find( label  );

  if  (p != vector_prop_labels_map.end() )

  {  

    temp    =  (p -> second) ;


    for (int i =0; i< temp.size();++i)
    {
      //    cout <<  "temp.size()  " <<  temp.size();
      //     return_vector[i]    = (int)(temp[i]); 
      return_vector.push_back( (int)(temp[i]) );

    }

    return  return_vector;
  
  }

  else
    cout  <<  "Warning --- read_input_vector int: vector data empty or  missing " << endl ;


}





//  END  NEW for  vectors  *********************
// ************************************************




void InputParser::read_data_maps( map <string,double>& num_map, map <string,string>&  string_map, 
				  map <string, vector<double> >&  vector_map) const

{

  num_map = prop_labels_map;
  string_map =string_prop_labels_map ;
  vector_map = vector_prop_labels_map ;


}







void InputParser::parse_device(ifstream& in_stream )

{

  // *******************************************************************************************************
  //  ***** IDEA:  List  of regions :
  //  $Regions$
  //  Reg 1 name {  reg_numb = 1  mat= Si  doping = 1e18 dop_type = donor}
  // .....................................................
  // $End$

  // $Bound_cond$
  // BC_name {  BC_numb = 1  type = dirichlet  value =  1.0 }
  //   $End$
  //
  //  vector<struct> ,   struct  ={  reg_n , reg_name, mat, dop_conc , dop type } 

  //   struct region_definition {

  //     unsigned int  reg_numb;
  //     string  reg_name;
  //     double  dop_conc  ;
  //     string dop_type;

  //   };

  //   region_definition reg_def;

  //   vector<reg_def> device_regions;

  vector<double> v_real; 
  vector<int> v_int;
  vector<string> v_string;

  string  name,str;

  RegionDefinition  current_region;  //   struct *************
  device_regions.clear();


  rule<>special_char =  (ch_p('_') | ch_p('-') |  ch_p('.') |  ch_p('/') | ch_p('+')  );
  rule<>region_name = (+alnum_p)>>   * ( (special_char ) >> *(+alnum_p) ) ;
  rule<>string_value     = (+alnum_p)>>   * (  (special_char )      >> *(+alnum_p) ) ;


  rule<>assignement_double  =   *(space_p) >> ch_p('=')>> *(space_p) >> 
    ( real_p[push_back_a(v_real)])>>*(space_p)   ; // with _ !!
  rule<>assignement_int  =   *(space_p) >> ch_p('=')>> *(space_p) >>
    ( int_p[push_back_a(v_int)])>>*(space_p)  ; // with _ !!
  rule<>assignement_string  =   *(space_p) >> ch_p('=')>> *(space_p) >>
    ( (string_value)[push_back_a(v_string)])>>*(space_p)  ; 


  //   string label1, label2, label3, label4;
  //   label1 = "reg_numb";
  //   label2 = "mat";

  // ******************************
  //  TO DO:  add label  "crystal_struct" 
  // ****************************
  //   label3 = "doping";
  //   label4 = "dop_type";

  rule<>label1_p = str_p("reg_numb");
  rule<>label2_p = str_p("mat");
  rule<>label3_p = str_p("crystal_struct");
  rule<>label4_p = str_p("doping");
  rule<>label5_p = str_p("dop_type");

  rule<> reg_prop = *(space_p) >>(label1_p) >> assignement_int >> 
    (label2_p) >>  assignement_string  >> (label3_p)>> 
    assignement_string  >> (label4_p) >> 
    assignement_double >>  (label5_p) >> assignement_string;   

  rule<> r_region = *(space_p)>>region_name[push_back_a(v_string)] >> 
    *(space_p) >> (ch_p('{'))>>  reg_prop >>*(space_p) >> (ch_p('}'))>> *(space_p) ; 


  while ( getline(in_stream, str) )
  {


    if  (!(parse(str.c_str(), comment_p("#")   , space_p).full) ) 

    { // if !  comment_p("#")  

	 

      if (  parse(str.c_str(),
                  r_region
                  )
            .full )  
      {
        //  //   cout << "v_real[0] " << v_real[0]<< endl;
        // 	      current_region.reg_name = v_string[0];
        // 	      current_region.reg_numb = v_int[0];
        // 	      current_region.mat_name = v_string[1];
        // 	      current_region.dop_conc = v_real[0];
        // 	      current_region.dop_type = v_string[2];

        // using  class  RegionDefinition 
        current_region.set_region_name(v_string[0]);
        current_region.set_region_number(v_int[0]);
        current_region.set_material_name(v_string[1]);
        current_region.set_crystal_name(v_string[2]);
        current_region.set_doping_concentration( v_real[0]);
        current_region.set_doping_type(v_string[3]);



        device_regions.push_back(current_region); // vector  of  RegionDefinition objects

        //   cout << "device_regions[0].mat_name  " <<device_regions[0].mat_name  << endl;

        //    cout << "current_region.mat_name " <<current_region.mat_name << endl;
        //    cout << v_int.size()<< endl;
        //   cout << "in parse_device: device_regions.size()  " << device_regions.size();


      }


      else  if (parse(str.c_str(), if_p("$")[(+alpha_p)[assign_a(name)] >> 
                                             ch_p("$")  ] , space_p ).full)


      {
        if (name == "End")
        {  
          cout << name  << endl ;
          break;
        }


      }

      else 
      {
        cerr <<  "  SYNTAX ERROR in input  file (section device structure)   " <<  endl;
        exit(1); 
      }





    }


    v_real.clear();
    v_int.clear();
    v_string.clear();


  }  //  end  while



}



const  vector<RegionDefinition>&   InputParser::get_device_regions()

{

  return  device_regions;


}





// void InputParser::get_device_data( vector<string>& reg_name_v,vector<unsigned int>& reg_numb_v,
// 				   vector<string>& mat_name_v, vector<double>& dop_conc_v, 
// 				   vector<string>&  dop_type_v   )

// {

//   //  cout << "device_regions.size()" << device_regions.size();

//   for (int i =0; i< device_regions.size();++i)
//     {

//       reg_name_v.push_back(device_regions[i].reg_name);
//       reg_numb_v.push_back(device_regions[i].reg_numb);
//       //  cout << "device_regions[i].mat_name" << 
//       //	device_regions[i].mat_name<< endl;
	

//       mat_name_v.push_back(device_regions[i].mat_name);
//       //  cout << " mat_name_v "<< mat_name_v[i]<< endl ; 

//       dop_conc_v.push_back(device_regions[i].dop_conc);
//       dop_type_v.push_back(device_regions[i].dop_type);

//     }


// }









void InputParser::parse_device_BC(ifstream& in_stream )

{

  // *******************************************************************************************************
  //  ***** IDEA:  List  of regions :
 

  // $Bound_cond$
  // BC_name {  BC_numb = 1  type = dirichlet  value =  1.0 }
  //   $End$
  //
 

  //   region_definition reg_def;

  //   vector<reg_def> device_regions;

  vector<double> v_real; 
  vector<int> v_int;
  vector<string> v_string;

  string  name,str;

  //  RegionDefinition  current_region;  //   struct *************
  //   device_regions.clear();


  rule<>special_char =  (ch_p('_') | ch_p('-') |  ch_p('.') |  ch_p('/') | ch_p('+')  );
  rule<>region_name = (+alnum_p)>>   * ( (special_char ) >> *(+alnum_p) ) ;
  rule<>string_value     = (+alnum_p)>>   * (  (special_char )      >> *(+alnum_p) ) ;


  rule<>assignement_double  =   *(space_p) >> ch_p('=')>> *(space_p) >> 
    ( real_p[push_back_a(v_real)])>>*(space_p)   ; // with _ !!
  rule<>assignement_int  =   *(space_p) >> ch_p('=')>> *(space_p) >>
    ( int_p[push_back_a(v_int)])>>*(space_p)  ; // with _ !!
  rule<>assignement_string  =   *(space_p) >> ch_p('=')>> *(space_p) >>
    ( (string_value)[push_back_a(v_string)])>>*(space_p)  ; 


  //   string label1, label2, label3, label4;
  //   label1 = "BC_reg_numb";
  //   label2 = "type";
  //   label3 = "value";
  //   // label4 = "dop_type";

  rule<>label1_p = str_p("BC_reg_numb");
  rule<>label2_p = str_p("type");
  rule<>label3_p = str_p("value");
  //  rule<>label4_p = str_p("dop_type");

  rule<> reg_prop = *(space_p) >>(label1_p) >> assignement_int >> 
    (label2_p) >>  assignement_string  >> (label3_p)>> 
    assignement_double; 
  // >> (label4_p) >> assignement_string;   

  rule<> r_region = *(space_p)>>region_name[push_back_a(v_string)] >> 
    *(space_p) >> (ch_p('{'))>>  reg_prop >>*(space_p) >> (ch_p('}'))>> *(space_p) ; 


  while ( getline(in_stream, str) )
  {


    if  (!(parse(str.c_str(), comment_p("#")   , space_p).full) ) 

    { // if !  comment_p("#")  

	 

      if (  parse(str.c_str(),
                  r_region
                  )
            .full )  
      {
        //   cout << "v_real[0] " << v_real[0]<< endl;

        BC_region_name_v.push_back(v_string[0]);
        BC_region_number_v.push_back(v_int[0]);
        BC_type_v.push_back(v_string[1]);
        BC_value_v.push_back(v_real[0]);



      }


      else  if (parse(str.c_str(), if_p("$")[(+alpha_p)[assign_a(name)] >> 
                                             ch_p("$")  ] , space_p ).full)


      {
        if (name == "End")
        {  
          cout << name  << endl ;
          break;
        }


      }

      else 
      {
        cerr <<  "  SYNTAX ERROR in input  file (section device boundary cond. )   " <<  endl;
        exit(1); 
      }

	 


    }


    v_real.clear();
    v_int.clear();
    v_string.clear();


  }  //  end  while



}




void InputParser::get_BC_data( vector<string>& BC_region_name_v_out,
			       vector<unsigned int>& BC_region_numb_v_out,
			       vector<string>& BC_type_v_out, vector<double>& BC_value_v_out  )


{


  BC_region_name_v_out=BC_region_name_v ;
  BC_region_numb_v_out=BC_region_number_v ;
  BC_type_v_out = BC_type_v;
  BC_value_v_out = BC_value_v;



}


void
InputParser::parse_alloy(ifstream& in_stream)
{// method  for   parsing  of  alloy model section

  // vector <AlloyModel>  alloy_model ;


  string  item, region_keyword, alloy_model_keyword, equal, model_type, start_symbol , end_symbol, x;
  double molar_fraction, start_molar_fraction, end_molar_fraction;
  unsigned int reg_numb;
  
  string dummy;
  char c;
  
  //alloy_model.clear();
  
  //AlloyModel current_alloy_model;
  
  
  
  // while (!in_stream.eof()) 
  
  //getline (in_stream, str))	  
  
  
  
  in_stream >> region_keyword;
  cout << region_keyword<< endl;
  
// if (region_keyword == "#")
// {
//  while (dummy != "\n")
//  {
//    c=getchar();
//  }
//  in_stream >> region_keyword;
// }
 
    
  while (  (region_keyword != "End") && (!in_stream.eof()) ) 
    
    
  {

    // Region n
    in_stream >>  reg_numb;
    in_stream >> start_symbol >> alloy_model_keyword >> equal >> model_type;
    alloy_model_pointer = new AlloyModel;

    //  if  (std::strncmp (model_type, "constant") = 0) 
    if (model_type == "constant")
    {

      in_stream >> x >> equal >> molar_fraction >> end_symbol;
      cout << " Region1  " << endl;
      alloy_model_pointer->set_x_constant (molar_fraction);
    }


    else
      //  if (std::strncmp (model_type, "linear") = 0)
      if (model_type =="linear")

      {
        in_stream >> x >> equal >> start_molar_fraction >>
          x >> equal >> end_molar_fraction >> end_symbol;                 

        alloy_model_pointer->set_x_min (start_molar_fraction);
        alloy_model_pointer->set_x_max (end_molar_fraction);
        cout << " Region2  " << endl;

      }

    alloy_model_pointer-> set_model (model_type);

    // alloy_model_pointer-> alloy_model_pointer = &current_alloy_model;


    // alloy_model.push_back (current_alloy_model);
    
    


    // map reg_alloy_model_map map <int reg_numb, AlloyModel&  current_alloy_model>
    cout << reg_numb << "   "  << model_type<< endl;
    reg_alloy_model_map.insert (make_pair (reg_numb, alloy_model_pointer));
    
    in_stream >> region_keyword;
    cout << region_keyword<< endl;

	      
  }// end while != "end"
      
  cout <<   "  Out of  while "<< endl;


  //  } // end while != "eof"

} //  end  method





const   map <unsigned int, AlloyModel*>&  
InputParser::get_alloy_model_map()
{
  
  return  reg_alloy_model_map;
  
}





