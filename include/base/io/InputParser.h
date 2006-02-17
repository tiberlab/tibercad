/*=============================================================================
  Copyright (c) 2002-2003 Joel de Guzman
  http://spirit.sourceforge.net/

  Use, modification and distribution is subject to the Boost Software
  License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
  http://www.boost.org/LICENSE_1_0.txt)
  =============================================================================*/


#ifndef _INPUTPARSER_H_
#define _INPUTPARSER_H_

#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/dynamic.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <iostream>
#include <fstream>

#include <map>
#include <vector>
#include <string>

#include "RegionDefinition.h"


///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace boost::spirit;

///////////////////////////////////////////////////////////////////////////////
//
// 
//
///////////////////////////////////////////////////////////////////////////////

//!  A parser  for  input  text  files . 
/*!
  Public method read_input allows  to  read  values in the  format tagname = value.
  These assignements can be placed everywhere inside the  section defined by "section_name" and should be  separated by  spaces.
  Everything following a '#' is  a  comment and  is  disregarded.
  A vector of  values can be  read in  the  same  way if  it  is  written in the  format: tagname = ( value1  value2 .... valueN ).  Only  one  vector  value  per  line  is   allowed.
  $/alpha^x = 3$

*/


class InputParser{


 public:


  //!  Constructor 
  /*!
    Defines  name  of input  file ("input_file_name").
  */

  //  InputParser(string filename , string section_name   );
  InputParser::InputParser(string& input_file_name);

  //!  Destructor 
  /*!
  
  */
  ~InputParser();


  //!  Overloaded  method  to  read input 
  /*!
    Returns   a  double if  "label"  is  found; otherwise returns Default .
  */

  double InputParser::read_input( string label , double  Default);

  //!  Overloaded  method  to  read input 
  /*!
    Returns   an  int  if  "label"  is  found; otherwise returns Default .
  */
  int  InputParser::read_input( string label , int  Default);
 
  //    bool   InputParser::read_input2( string label_bool , bool  default_bool );


  //!  Overloaded  method  to  read input 
  /*!
    Returns   a string   if  "label"  is  found; otherwise returns Default .
  */
  string   InputParser::read_input( string label , string  Default);

  //  bool   InputParser::read_input( string label_bool , bool  default_bool );

  //!  Overloaded  method  to  read input 
  /*!
    Returns   a  boolean value  if  "label_bool"  is  found .
  */
  bool   InputParser::read_input( string label_bool );



  /*   void InputParser::read_input_vector( string label , vector<double>& return_vector); */
  /*   void InputParser::read_input_vector( string label , vector<int>& return_vector); */



  //!  Overloaded  method  to  read input 
  /*!
    Returns   a  vector of  double   if  "label"  is  found; return_vector is  a  dummy vector
  */
  vector<double>  InputParser::read_input( string label, vector<double> return_vector );


  //!  Overloaded  method  to  read input 
  /*!
    Returns   a  vector of  int   if  "label"  is  found; def_vector  is  a  dummy vector
  */
  vector<unsigned int> InputParser::read_input( string label, vector<unsigned int> def_vector);





  // void  InputParser::get_data ( vector< vector<double> >& glob_reg_values,
  //  vector<int>& glob_reg_id,  vector<string>& glob_mat);

  // OLD  !!!
  void  InputParser::get_data ( vector< vector<double> >& glob_reg_values,
				vector< vector<double> >& glob_comm_values,
				vector<int>& glob_reg_id,  vector<string>& glob_mat );


 


  //  void InputParser::get_device_data( vector<RegionDefinition>& device_regions_struct );

  //!    Method  to  read device structure input 
  /*!
    Returns   five  vectors  with values of fields for  each physical region of the device: region name,
    region number, material name, doping concentration, doping  type.
  */
  void InputParser::get_device_data( vector<string>& reg_name_v,vector<unsigned int>&  reg_numb_v,
				     vector<string>&  mat_name_v, vector<double>& dop_conc_v, 
				     vector<string>&  dop_type_v   );

  void InputParser::read_section(string& section_name); 

  void InputParser::get_BC_data( vector<string>& BC_region_name_v_out,
				 vector<unsigned int>& BC_region_numb_v_out,
				 vector<string>& BC_type_v_out, vector<double>& BC_value_v_out  );

  void InputParser::read_data_maps( map <string,double>& num_map, map <string,string>&  string_map, 
				    map <string, vector<double> >&  vector_map) const  ;

  const  vector<RegionDefinition>&   InputParser::get_device_regions();


  //----------------------------------------------------------------------------------

 private:


  // string  section_name, filename;
  string  filename;

  vector< vector<double> > reg_values  ;

  vector< vector<double> > command_values  ;

 /*  struct RegionDefinition { */

/*     string reg_name ; */
/*     unsigned int  reg_numb; */
/*     string  mat_name; */
/*     double  dop_conc  ; */
/*     string dop_type; */

/*   }; */

// Use directly  class  RegionDefinition


  vector<RegionDefinition> device_regions;   //  vector of  objects  RegionDefinition

  //  BC regions vectors
  // ------------------------------------
  vector<string> BC_region_name_v;
  vector<unsigned int> BC_region_number_v;
  vector<string> BC_type_v;
  vector<double> BC_value_v;
  // ---------------------------------


  vector<int> reg_id;
  vector<string> mat;


  map <string,double>  prop_labels_map;
  map <string,string>  string_prop_labels_map;



  vector<string> prop_labels ;

  vector<string> vector_prop_labels;
  vector<string> string_prop_labels;



  map <string, vector<double> >  vector_prop_labels_map;



  //  vector<double> values;
  /*  vector< vector<double> > reg_values  ; */


  /*   vector<int> reg_id; */
  /*  vector<string> mat; */

  void InputParser::initialize_vectors();

 
  void InputParser::parse_options(ifstream& in_stream );

  void InputParser::parse_device(ifstream& in_stream );

  void InputParser::parse_device_BC(ifstream& in_stream );

};



#endif // endif define   _INPUTPARSER_H_
