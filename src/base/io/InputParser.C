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



#include "InputParser.h"


///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace boost::spirit;

///////////////////////////////////////////////////////////////////////////////
//
//  Our comma separated list parser
//
///////////////////////////////////////////////////////////////////////////////

InputParser::InputParser(string filename)
{

  initialize_vectors(); // 
//  reg_values.clear();
  scan_input(filename);
 // counter = 0;

  //  a = x;

  //std::ifstream  in_stream (file_name.c_str());

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
    // counter = 0;

}


// *********************************
// NEW: COMMAND  SECTION
// *********************************

void InputParser::parse_command(ifstream& in_stream )


{

  // *******************************************************************************************************
  //  ***** IDEA:  nome=valore  ->  map <string,double>  o map <string,string>  , map <string, bool> 
  // variable = inp_file.read("prop", int Default)
  //   int  inp_file::read ( string,  int  Default);
  //    find string   in  map  <string,prop_value>  ->  return  value 
  //  name_property  ->  if  find in  vector< name_property>  then  property  = true
  // *******************************************************************************************************



 
  string name, str, mat_name ;
  name = "Start";
  int id;

  // **********************
  //  assigned   values (local)
  // ****************************
  vector<double> v;              //  scalar  and   vector  values
  //  vector<double>  vect;
  vector<string>  v_bool;
  v_bool.clear();
  //  vect.clear();
  v.clear();
  //******************


  //------------------------
  //  labels of  properties
  //------------------------

    vector<string> v_label;
    vector<string>vect_label;
    vector<string> v_label_bool;  //  label  for  string  value
    v_label_bool.clear();
    vect_label.clear();
    v_label.clear();

    //  vect.push_back(0);

   
    // *************************************************************************************************
    //  rule<phrase_scanner_t>  for  phrase level parsing (e.g. with a  separator character (comma) )
    // rule <> for  space  separated  list (character level parsing) : NEEDS  TO READ  EXPLICITLY SPACES !!!
    // *************************************************************************************************

 //  NEW  *********************
    // rule<>list_of_numbers_space_sep = real_p[push_back_a(vect)] >> *( *(space_p) >> real_p[push_back_a(vect)]);
    //  try  with  ( 4.5  6.7  8.98)
    rule<>list_of_numbers_space_sep = ch_p('(')>> *(space_p) >>real_p[push_back_a(v)] >> *( *(space_p) >> real_p[push_back_a(v)])>> *(space_p)>> ch_p(')')  ;

    rule<>special_char =  (ch_p('_') | ch_p('-') |  ch_p('.') |  ch_p('/') );

    // rule<>label  = (+alnum_p)>>   * ( ch_p('_') >> *(+alnum_p) ) ;
 rule<>label  = (+alnum_p)>>   * ( (special_char ) >> *(+alnum_p) ) ;



 // rule<>label_file  = (+alnum_p)>>   * ( (( ch_p('_'))|( ch_p('.') ) | ( ch_p('/') ) )    >> *(+alnum_p) ) ;
 rule<>label_file  = (+alnum_p)>>   * (  (special_char )      >> *(+alnum_p) ) ;





   
 // rule<>assignement  =  (+alnum_p)[push_back_a(v_label)] >> *(space_p) >> ch_p('=')>> *(space_p) >> ( real_p[push_back_a(v)])  ; //ok !!
 rule<>assignement  =  (label)[push_back_a(v_label)] >> *(space_p) >> ch_p('=')>> *(space_p) >> ( real_p[push_back_a(v)])  ; // with _ !!




 //  label  instead  of +alnum_p !!! 
 // rule<>assignement_bool  =  (label)[push_back_a(v_label_bool)] >> *(space_p) >> ch_p('=')>> *(space_p) >> ((+alnum_p)[push_back_a(v_bool)])   ;

 rule<>assignement_bool  =  (label)[push_back_a(v_label_bool)] >> *(space_p) >> ch_p('=')>> *(space_p) >> ((label_file )[push_back_a(v_bool)])   ;



 rule<>assignement_vector =  (label)[push_back_a(vect_label)] >> *(space_p) >> ch_p('=')>> *(space_p) >> list_of_numbers_space_sep;

 // rule<>list_of_assignement = assignement >> *( *(space_p) >> assignement  );
 rule<>list_of_assignement = (assignement | assignement_bool)  >> *( *(space_p) >> (assignement| assignement_bool)       );


  //rule<>

 // rule<> r_command  = *(space_p) >> (list_of_assignement | assignement_vector)    >> *(anychar_p);

 rule<> r_command  = *(space_p) >> (list_of_assignement )    >> *(anychar_p);

 rule<> r_command_vector  = *(space_p) >> (assignement_vector)    >> *(anychar_p);


 // ******************
 //  END  NEW  
 // ********************


//  //  rule<phrase_scanner_t>composition = ch_p('(') >> real_p >> ch_p(')');
// //   rule<phrase_scanner_t>material = (+alnum_p)>> *(composition) >> *(alnum_p);

// //  rules  for  space  separated  list  (character level parsing)  !!!!!
//  rule<>composition = ch_p('(') >> real_p >> ch_p(')');
//   rule<>material = (+alnum_p)>> *(composition) >> *(alnum_p);



//   //ok !!  //  rule<phrase_scanner_t>list_of_numbers = real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)]);
//   //  rule<phrase_scanner_t> r = uint_p >> ch_p(',')>> material[assign_a(name)] >> *(ch_p(','))>>list_of_numbers >> *(ch_p('#'))>>*(alpha_p);

 
//   //  space  separated  list (instead  of  comma  separated)
//   // rule<>list_of_numbers_space_sep = real_p[push_back_a(v)] >> *( *(space_p) >> real_p[push_back_a(v)]);
//  rule<> r = *(space_p)   >> uint_p[assign_a(id)] >> *(space_p)   >> material[assign_a(mat_name)] >> *(space_p)>>list_of_numbers_space_sep  >> *(anychar_p);

//     //(space_p)>>*(ch_p('#'))>> *(alpha_p);




//  // rule<phrase_scanner_t>r = real_p >> *(ch_p(',') >> real_p);

//  //  OK  ***************
//  //  getline(in_stream, str);

//   //  while (!(str.empty()))
//  //  OK  ***************



    while ( getline(in_stream, str) )
    {

      //    if ( (!(parse(str.c_str(), *(ch_p('#')), space_p).full) )   /// not a comment !! comment_p("//")
      if ( (!(parse(str.c_str(), comment_p("#")   , space_p).full) ) 
	   && (! (parse(str.c_str(), (*alpha_p)[assign_a(name)]   , space_p).full) ) ) // not  End !



	{ // if !  comment_p("#")  

	  //  cout << name<< "  ???????????    " << endl;

	  if (  parse(str.c_str(),

		      //  Begin  grammar
		      // (
		       //     real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)])
		       //  (+alpha_p)[assign_a(name)] >> real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)])

		       //  (+alnum_p)[assign_a(name)] >> *(ch_p(','))>> real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)]) >> *(ch_p('#'))>>*(alpha_p) 
			//  (+alnum_p)[assign_a(name)] >> *(ch_p(','))>> list_of_numbers >> *(ch_p('#'))>>*(alpha_p)
		      //	r
		      r_command 
		       )

		      //  ,
		      //  End grammar

		      //  space_p).full )
		      .full )           //  not skipping spaces


	    { // if parse

	      //   cout << name<< "  material    " << endl;
	      // counter++;

	     
	      //   reg_values.push_back(v);
	      //   command_values.push_back(v);
	      //    cout << "*******" << v[0]<< v[1] << v[2]<<endl;


	      ////     cout << "labels   " << v_label[0]<< endl ;  //v_label[1] << v_label[2]<< endl;
	      ////	     cout << "labels_bool   " << v_label_bool[0]<< endl ;  //v_label[1] << v_label[2]<< endl;

	      ////	     cout << " ************* v_label.size(), v_label_bool.size()   " << v_label.size()<<"  "<< v_label_bool.size()<< endl;

	      for (int i =0; i< v.size();++i)
		{
		  prop_labels.push_back(v_label[i]);
		  prop_labels_map.insert(make_pair(v_label[i], v[i]) );
		}

	      for (int i =0; i< v_label_bool.size();++i)
		{
		  bool_prop_labels.push_back(v_label_bool[i]);
		  bool_prop_labels_map.insert(make_pair(v_label_bool[i], v_bool[i]) );
		}

	     


	      //   cout << "scalar  data   " << v[0] << "  " << v[1]<< "  "  << v[2]<< endl ; //<< "  "  << vect[3] <<  endl ;
	      //   cout  <<  "  v_bool ---------  "  << v_bool[0]<<  endl ;

		     //   cout << "vector  data   " << vect[0] ;

	      //  mat.push_back(mat_name);
	      //  reg_id.push_back(id);


	      //  cout << reg_values.size();
	      //   cout << reg_values[0][0] << " reg_values[0][0]     " << endl;

	   //    v.clear();
// 	      v_label.clear();
// 	      //   vect.clear();
// 	      vect_label.clear();
// 	      v_label_bool.clear();
// 	      v_bool.clear();


	      v_label_bool.clear();
	      vect_label.clear();
	      v_label.clear();
	      v_bool.clear();
 	      v.clear();



	     //  mat_name= "";
// 	      id = 0;

	    }

	  else if(  parse(str.c_str(),  r_command_vector   )   .full )           //  not skipping spaces
	    //  reads  list of  values (vector)

	    {


	      vector_prop_labels.push_back(vect_label[0]);
	      //  for (int i =0; i< v.size();++i)
	      //  	{
		 
	      vector_prop_labels_map.insert(make_pair(vect_label[0], v) );

		  // 	}

	      ////	     cout << " v_label[0]  " << v_label[0]<< endl ; 
	      ////   cout << " vect_label  " << vect_label[0]<< endl ;
	      ////    cout << "vector  data   " << v[0] << "  " << v[1]<< "  "  << v[2]<< endl ; //<< "  "  << vect[3] <<  endl ;

	     //  v.clear();
// 	      v_label.clear();
// 	      vect_label.clear();

	      
	      v_label_bool.clear();
	      vect_label.clear();
	      v_label.clear();
	      v_bool.clear();
 	      v.clear();





	    }


	}


     

      if (name == "End")
	{  
	  //cout << name  << endl ;
	  break;
	}

 //  OK  ***************
    //   else
// 	{
// 	  getline(in_stream, str);
// 	}
 //  OK  ***************



      

    }  //  end  while


}


// *********************************
// END  NEW: COMMAND  SECTION
// *********************************



// ******************************************
//  CASE WITH  SPACE SEPARATED  LIST !!
// ******************************************
//void InputParser::parse_1(ifstream& in_stream )
void InputParser::parse_regions(ifstream& in_stream )


  //void parse_1(ifstream& in_stream )

{

 

  string name, str, mat_name ;
  name = "Start";
  int id;

    vector<double> v;
    v.clear();

  
   
 

 //  rule<phrase_scanner_t>composition = ch_p('(') >> real_p >> ch_p(')');
//   rule<phrase_scanner_t>material = (+alnum_p)>> *(composition) >> *(alnum_p);

//  rules  for  space  separated  list  (character level parsing)  !!!!!
 rule<>composition = ch_p('(') >> real_p >> ch_p(')');
  rule<>material = (+alnum_p)>> *(composition) >> *(alnum_p);



  //ok !!  //  rule<phrase_scanner_t>list_of_numbers = real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)]);
  //  rule<phrase_scanner_t> r = uint_p >> ch_p(',')>> material[assign_a(name)] >> *(ch_p(','))>>list_of_numbers >> *(ch_p('#'))>>*(alpha_p);

 
  //  space  separated  list (instead  of  comma  separated)
 rule<>list_of_numbers_space_sep = real_p[push_back_a(v)] >> *( *(space_p) >> real_p[push_back_a(v)]);
 rule<> r = *(space_p)   >> uint_p[assign_a(id)] >> *(space_p)   >> material[assign_a(mat_name)] >> *(space_p)>>list_of_numbers_space_sep  >> *(anychar_p);

    //(space_p)>>*(ch_p('#'))>> *(alpha_p);




 // rule<phrase_scanner_t>r = real_p >> *(ch_p(',') >> real_p);

  getline(in_stream, str);

  while (!(str.empty()))
    {

      //    if ( (!(parse(str.c_str(), *(ch_p('#')), space_p).full) )   /// not a comment !! comment_p("//")
      if ( (!(parse(str.c_str(), comment_p("#")   , space_p).full) ) 
	   && (! (parse(str.c_str(), (*alpha_p)[assign_a(name)]   , space_p).full) ) ) // not  End !



	{

	  //  cout << name<< "  ???????????    " << endl;

	  if (  parse(str.c_str(),

		      //  Begin  grammar
		      // (
		       //     real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)])
		       //  (+alpha_p)[assign_a(name)] >> real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)])

		       //  (+alnum_p)[assign_a(name)] >> *(ch_p(','))>> real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)]) >> *(ch_p('#'))>>*(alpha_p) 
			//  (+alnum_p)[assign_a(name)] >> *(ch_p(','))>> list_of_numbers >> *(ch_p('#'))>>*(alpha_p)
			r

		       )
		      //  ,
		      //  End grammar

		      //  space_p).full )
		      .full )           //  not skipping spaces


	    {

	      //   cout << name<< "  material    " << endl;
	      // counter++;

	     
	      reg_values.push_back(v);
	      mat.push_back(mat_name);
	      reg_id.push_back(id);


	      //  cout << reg_values.size();
	      //   cout << reg_values[0][0] << " reg_values[0][0]     " << endl;
	      v.clear();
	      mat_name= "";
	      id = 0;

	    }


	}


     

      if (name == "End")
	{  
	  cout << name  << endl ;
	  break;
	}
      else
	{
	  getline(in_stream, str);
	}



      

    }  //  end  while




  // cout << reg_values[1].size() << endl;
   

 //  for (vector<double>::size_type i = 0; i < v.size(); ++i)
 //               cout << i << ": " << v[i] << endl;


//   for (int i =0; i< reg_values.size();++i)
// {

//   cout <<  "reg # " <<  reg_id[i]<< "   " <<  mat[i]<< "  material    " << endl;
//    for (int j = 0; j< reg_values[i].size();++j)
//   cout << i<< "  "<< j << ": " << reg_values[i][j] << endl;

// }



}




// ******************************************
//  CASE WITH  COMMA SEPARATED  LIST !!
// *********************************
//void InputParser::parse_2(ifstream& in_stream )
void InputParser::parse_comma(ifstream& in_stream )


{

 

  string name, str;
  name = "Start";

    vector<double> v;

  getline(in_stream, str);

  while (!(str.empty()))
    {

      //    if (!(parse(str.c_str(), *(ch_p('#')), space_p).full) )   /// not a comment !!

      if ( (!(parse(str.c_str(), *(ch_p('#')), space_p).full) )   /// not a comment !!
	   && (! (parse(str.c_str(), (*alpha_p)[assign_a(name)]   , space_p).full) ) )



	{

	  if (  parse(str.c_str(),

		      //  Begin  grammar
		      (
		       //     real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)])
		       //  (+alpha_p)[assign_a(name)] >> real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)])

		       (real_p[push_back_a(v)])  >>  *(space_p) >> real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)]) >> *(ch_p('#'))>>*(alpha_p) 


		       )
		      ,
		      //  End grammar

		      space_p).full )
	    {

	      cout << name<< "  B.C.    " << endl;

	


	    }


	}


     

      if (name == "End")
	{  
	  cout << name  << endl ;
	  break;
	}
      else
	{
	  getline(in_stream, str);
	}



      

    }  //  end  while




   

  for (vector<double>::size_type i = 0; i < v.size(); ++i)
                cout << i << ": " << v[i] << endl;


}





// parse  BC  section with  space  separation
void InputParser::parse_BC(ifstream& in_stream )

 

{

  int id;

  string name, str, BC_type;
  name = "Start";

    vector<double> v;

    rule<> BC_type_name = (str_p("Dir") | str_p("Neum") | str_p("Mixed"));

  //  space  separated  list (instead  of  comma  separated)
 rule<>list_of_numbers_space_sep = real_p[push_back_a(v)] >> *( *(space_p) >> real_p[push_back_a(v)]);
 // rule<> r = *(space_p)   >> uint_p[assign_a(id)] >> *(space_p)   >>  list_of_numbers_space_sep  >> *(space_p) >>(+alpha_p)>> *(comment_p("#"));
 rule<> r = *(space_p)   >> uint_p[assign_a(id)] >> *(space_p)   >>  list_of_numbers_space_sep  >> *(space_p) >> (+BC_type_name)[assign_a(BC_type)]  >> *(space_p) >>  *(comment_p("#"));





  getline(in_stream, str);

  while (!(str.empty()))
    {

      //    if (!(parse(str.c_str(), *(ch_p('#')), space_p).full) )   /// not a comment !!

      if ( (!(parse(str.c_str(), *(ch_p('#')), space_p).full) )   /// not a comment !!
	   && (! (parse(str.c_str(), (*alpha_p)[assign_a(name)]   , space_p).full) ) )



	{

	  if (  parse(str.c_str(),

		      //  Begin  grammar
		      //  (
		       //     real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)])
		       //  (+alpha_p)[assign_a(name)] >> real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)])

		      //   (real_p[push_back_a(v)])  >>  *(space_p) >> real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)]) >> *(ch_p('#'))>>*(alpha_p) 


		      //  )
		      // ,
		      r )
		      //  End grammar

		//   space_p)
		.full )
	    {

	      cout << name<< "  B.C.    "<< BC_type << endl;

	


	    }


	}


     

      if (name == "End")
	{  
	  cout << name  << endl ;
	  break;
	}
      else
	{
	  getline(in_stream, str);
	}



      

    }  //  end  while




   

  for (vector<double>::size_type i = 0; i < v.size(); ++i)
                cout << i << ": " << v[i] << endl;


}








void InputParser::read_data_section(char const* str,ifstream& in_stream )
// void InputParser::read_data_section(char const* str)


   // void read_data_section(char const* str,ifstream& in_stream )

{

 string name;

 if  (parse(str, if_p("$")[(+alpha_p)[assign_a(name)]].else_p[nothing_p]  , space_p).full)
   {
     //cout << name<< endl;
     if (name == "Regions")
       //	if (name == name_data )
       {
	 //  parse_1(in_stream);
	 parse_regions(in_stream);

       }

     else if  (name == "Boundary")
       {
	 //  parse_2(in_stream);
	 // parse_comma(in_stream);
	 parse_BC(in_stream);

       }

     else if  (name == "Options")
       {
	 //  parse_2(in_stream);
	 // parse_comma(in_stream);
	 parse_command(in_stream);

       }






   }



}





void InputParser::scan_input(string file_name)

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
  //cout << "\t\t read  from file in  section Data section \n\n";
  //
  //cout << "/////////////////////////////////////////////////////////\n\n";


  getline(in_stream, str);
  while (str.empty() )
    {
      getline(in_stream, str);
    }

    read_data_section(str.c_str(),in_stream);
    //  read_data_section(str.c_str());


  while (getline(in_stream, str))  //   
    {
           read_data_section(str.c_str(),in_stream);
	   //  read_data_section(str.c_str());


    }





}


void  InputParser::get_data ( vector< vector<double> >& glob_reg_values, vector< vector<double> >& glob_comm_values  , vector<int>& glob_reg_id,  vector<string>& glob_mat    )
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
	

	  //  cout << endl <<  "BC #"<< BC_reg_numbers[i] << "   " << endl; 
	  // 	for (int i =0; i< temp.size();++i)
	  // 		  {

	  // 		    if ( (i % 3) ==  0)     cout << endl; 

	  // 		    cout <<  temp[i] << "   " ;  //<< endl; 

	  // 		  }
	  // 		cout << endl; 
	  return value ;

	}

      else 
	{ 
	  //cout  <<  "*** Default ***  ";

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
	  //cout  <<  "*** Default ***  ";
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
      p = bool_prop_labels_map.find( label  );

      if  (p != bool_prop_labels_map.end() )

	{  
	  value  =  (p -> second) ;
	

	 
	  return value ;

	}

      else
	{
	  //cout  <<  "*** Default ***  ";
	  return Default;
	  //  cout  <<  "error"  ;
	}
      //   }

}



// // *****  bool  return
// bool   InputParser::read_input( string label , bool  Default)
// {
//   string  value;
//   map <string,string>  :: iterator  p;
//   bool true;
//   bool  false;




//   //  for (int i =0; i< prop_labels.size();++i)
//   //    {
//   p = bool_prop_labels_map.find( label  );

//   if  (p != bool_prop_labels_map.end() )

//     {  
//       value  =  (p -> second) ;
	 
//       //	if ( strncmp (str[n],"C3**",2) == 0)

//       if ( (strncmp ( value ,"true") == 0) )

// 	{ 
// 	  return true;

// 	}
//       else if  ( (strncmp ( value ,"false") == 0) )
// 	{ 
// 	  return false;

// 	}


	 
//       return value ;

//     }

//   else
//     {
//       cout  <<  "*** Default ***  ";
//       return Default;
//       //  cout  <<  "error"  ;
//     }
//   //   }

// }







void InputParser::read_input_vector( string label , vector<double>& return_vector)

{
 
  map <string, vector<double> >   :: iterator  p;


  //  for (int i =0; i< vector_prop_labels.size();++i)
  //   {
  p = vector_prop_labels_map.find( label  );

  if  (p != vector_prop_labels_map.end() )

    {  
      return_vector    =  (p -> second) ;
	


    }

  else
    cout  <<  "error --- read_input_vector  " << endl ;
  //   }

}




void InputParser::read_input_vector( string label , vector<int>& return_vector)

{
 
  map <string, vector<double> >   :: iterator  p;

  vector<double> temp;
  return_vector.clear();

  //  for (int i =0; i< vector_prop_labels.size();++i)
  //   {
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

	
    }

  else
    cout  <<  "error --- read_input_vector  "  ;
  //  }

}
