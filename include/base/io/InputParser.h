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

#include <map>
#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace boost::spirit;

///////////////////////////////////////////////////////////////////////////////
//
//  Our comma separated list parser
//
///////////////////////////////////////////////////////////////////////////////

class InputParser{


 

 

  // private:

  vector< vector<double> > reg_values  ;

  vector< vector<double> > command_values  ;



  vector<int> reg_id;
 vector<string> mat;


   map <string,double>  prop_labels_map;
   map <string,string>  bool_prop_labels_map;



 vector<string> prop_labels ;

 vector<string> vector_prop_labels;
 vector<string> bool_prop_labels;



 map <string, vector<double> >  vector_prop_labels_map;



 //  vector<double> values;
 /*  vector< vector<double> > reg_values  ; */


/*   vector<int> reg_id; */
/*  vector<string> mat; */

 void InputParser::initialize_vectors();

 //  void InputParser::parse_1(ifstream& in_stream );
 void InputParser::parse_regions(ifstream& in_stream );


 void InputParser::parse_command(ifstream& in_stream );



 


  //  void InputParser::parse_2(ifstream& in_stream );
  void InputParser::parse_comma(ifstream& in_stream );


  void InputParser::parse_BC(ifstream& in_stream );





  void InputParser::read_data_section(char const* str,ifstream& in_stream );

   // void read_data_section(char const* str,ifstream& in_stream )


  //  void InputParser::scan_input(ifstream& in_stream);

  void InputParser::scan_input(string file_name);


public:

 InputParser(string filename);
  ~InputParser();


  double InputParser::read_input( string label , double  Default);
  int  InputParser::read_input( string label , int  Default);
 

  string   InputParser::read_input( string label , string  Default);




  void InputParser::read_input_vector( string label , vector<double>& return_vector);
  void InputParser::read_input_vector( string label , vector<int>& return_vector);




  // void  InputParser::get_data ( vector< vector<double> >& glob_reg_values,  vector<int>& glob_reg_id,  vector<string>& glob_mat);
void  InputParser::get_data ( vector< vector<double> >& glob_reg_values, vector< vector<double> >& glob_comm_values  , vector<int>& glob_reg_id,  vector<string>& glob_mat );

};



