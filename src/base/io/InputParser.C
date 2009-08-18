// $Id$

#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/dynamic.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>
#include <cstring>

#include "InputParser.h"




///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace boost::spirit;




InputParser::InputParser(const std::string& input_file_name)

{

  start_symb = '{';
  end_symb  =  '}';

 
  filename = input_file_name;
  ifstream in_stream(filename.c_str());
  if (!in_stream.good())
    throw InitFailedException("Bad input file.");

  reset_all_maps();
  
}

InputParser::~InputParser(void)
{
}




// private  method: utility  to  find a  keyword in a section
bool InputParser::find_keyword_in_section(ifstream& in_stream, const std::string& keyword)

{

  bool found = false;
  std::string str,  label,dollar_symbol ;
  dollar_symbol = "$";
  char ch ;

  assert(in_stream.good());
  in_stream >>  label;

  //  it can be a comment  ! 
  //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
  while (skip_comments(in_stream,label) == true )
  {
    in_stream >> label; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
  } 


 
  cut_off_comment(label, in_stream); //  erase comment in  case label#my_comment....
  //  cerr <<  " cut_off_comment(label)  "  <<  label <<  endl;

  
  //       while (  ( found == false) && (  label != dollar_symbol ) && (!in_stream.eof()) ) //  
  while (  ( found == false) && (  std::strncmp ((label.c_str()),dollar_symbol.c_str(),1) != 0  )
           && (!in_stream.eof()) )   
 
  { //while
        
    //   cerr <<  " +++++++++keyword =  " << keyword << " ++++++++label  =  " << label <<  endl;

    if  (label == keyword  )
    { // if label  == keyword

      //cerr <<  " -----------keyword =  " << keyword << " ----------label  =  " << label <<  endl;
      //*********************************************
        //NEW :  check if  it is   really  beginning of  a block   
        // 

        //  *****************************************************

        in_stream.get(ch); // get next char
      do{

        //*******************************


          if (ch == '#' )  //  skip  comments 
        { 
          in_stream.ignore(256,'\n');  //  if  the   read keyword is # or  begins with #: ignore  all  the  line 
          //      and   read  the  next keyword !!!
        }
        in_stream.get(ch);


      } while ((ch == '\n') || (ch == '\r') || (ch == ' ')
               || (ch == '#') || (ch == '\t'));


      // then  check  if  first char (not  blank) is  "{"
      if (ch == '{')
      {
        //   cout  <<  "  ch { :   "<< ch <<  endl << endl;
        in_stream.putback(ch );

        found = true;
        break;

        //  OK

      }
      else in_stream.putback(ch );

    }  // end if 

    else   skip_block(in_stream);  //  if the label != keyword, skip the  whole  block 
                                   //  in {} (to avoid confusion with labels)



    in_stream >>  label;
    // cerr <<  " BEFORE cut_off_comment(label)  "  <<  label <<  endl;

    //  it can be a comment  ! 
    //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
    while (skip_comments(in_stream,label) == true )
    {
      in_stream >> label; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
    } 


    cut_off_comment(label, in_stream); //  erase comment in  case label#my_comment....
    // cerr <<  " AFTER cut_off_comment(label)  "  <<  label <<  endl;

    // -> next  check on  label 

  }   //  end  while


  //   if (found == false)

  //   {
  //     cerr << " Error: keyword "<< keyword <<  "  not  found !  " << endl ;
  //     exit(1);
   
  //   }


  return  found;

}



// public  method to  add  parameters  for a given model " model_name"

const  ModelOptions& InputParser::read_parameters(std::string section_name, const std::string& model_name)

{

  std::string  label ;
  std::ifstream in_stream (filename.c_str()) ;
  char ch;

  bool found_model, check_error ;
  found_model = false;

  reset_all_maps();

  assert(in_stream.good());

  section_name = "$"+section_name;

  find_keyword( in_stream,section_name  );
 
  //  // { 
 
  //   // while (  ( label  != end_symbol) && (!in_stream.eof()) ) //   
  //   //        {   

  //   //   if  ( !(  model_name == "") )
  //   //   {    

  // in_stream >>  label;   //   read   start_symb  !!
  //   while (skip_comments(in_stream,label) == true )
  //   {
  //     in_stream >> label; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
  //   } 

  // ****************************************************************
  // get comments and  blanks until char;  if  char = '{' then stop


  //       in_stream.get(ch); // get next char
  //       // cerr <<  " in_stream.get(ch)" <<  ch  <<  endl;
  //       while (ch == '\n')  //  if  starts with new line !!
  //       {
  //         cout  <<  "  ch = new line   "<<  endl << endl;
  //         in_stream.get(ch);
  //       }



  //       while(ch == ' ')        //  skip  blanks
  //       {
  //         in_stream.get(ch);

  //         //  OR  ch != ' ' ,  OR  ch = EOL
  //         while (ch == '\n')
  //         {
  //           cout  <<  "  ch = new line   "<<  endl << endl;
  //           in_stream.get(ch);
  //         }

  //         cout  <<  "  ch:   "<< ch <<  endl << endl;
  //       }

  //       while (ch == '#' )  //  skip  comments   // ***  if ?
  //       {
  //         in_stream.ignore(256,'\n');  //  if  the   read keyword is # or  begins with #: ignore  all  the  line 
  //         //   in_stream >> item;  //      and   read  the  next keyword !!!
  //         in_stream.get(ch);

  //         while(ch == ' ') //  skip  blanks
  //         {



  //           in_stream.get(ch);
  //           //  OR  ch != ' ' ,  OR  ch = EOL
  //           while (ch == '\n')
  //           {
  //             cout  <<  "  ch = new line   "<<  endl << endl;
  //             in_stream.get(ch);
  //           }

  //           cout  <<  "  ch:   "<< ch <<  endl << endl;
  //         }

  //       } 


  //       // then  check  if  first char (not  blank) is  "{"
  //       if (ch == '{')
  //       {
  //         cout  <<  "  ch { :   "<< ch <<  endl << endl;
  //        //  OK

  //       }

  //       else cerr <<  "ERROR  hereee"  << endl; 



  // ++++++++++++++  put   function skip_to_bracket !!!
  // bool InputParser::skip_to_bracket(ifstream& in_stream)

  check_error = skip_to_bracket(in_stream);  // go on reading until the  char  '{' is  read  

  if (check_error == true)

  {
    std::ostringstream stm;

    stm << "In input file: missing  block at the beginning of section " << section_name << endl;
    throw InitFailedException(stm.str());
  }



  // -----------------------------------------------------
  //  OK  ----------------
  //

  //       in_stream.get(ch); // get next char
  //   do{
      
  //     if  ( (ch == ' ') || (ch == '\n'))
  //     {
  //       in_stream.get(ch);
  //     }

  //     if (ch == '#' )  //  skip  comments 
  //     { 
  //       in_stream.ignore(256,'\n');  //  if  the   read keyword is # or  begins with #: ignore  all  the  line 
  //       //      and   read  the  next keyword !!!
  //       in_stream.get(ch);
             
  //     }


  //   } while ( (ch == '\n') || (ch == ' ') || (ch == '#') );


  //   // then  check  if  first char (not  blank) is  "{"
  //   if (ch == '{')
  //   {
  //     cout  <<  "  ch { :   "<< ch <<  endl << endl;
  //     //  OK

  //   }

  //   else cerr <<  "ERROR  hereee"  << endl; 




  // **********************************************


 


  //  find_keyword_in_section( in_stream,  model_name);  //   read   Model name  string
  found_model =  find_keyword_in_section( in_stream,  model_name);  //   read   Model name  string

    
  //  in_stream >>  label;  //    {  read  by   parse_options) !!
            
  //*************************************
     
      //    ModelOptions temp_options;
      //*************************************
         
          temp_options.clear();


  


  if (!found_model)
  {
    return  temp_options;  //  model keyword  not  found-> returns empty ModelOptions !!

  } 



 

  parse_options(in_stream,temp_options,  model_name, section_name );

 

  return  get_options();

 
 
}


//  overload  to  read model-independent  parameters

const  ModelOptions&  InputParser::read_parameters(std::string section_name)

{

  std::string  label ;
  std::ifstream in_stream(filename.c_str()) ;

  reset_all_maps();

  assert(in_stream.good());

  section_name = "$"+section_name;

  find_keyword( in_stream,section_name  );
 
 
 

  temp_options.clear();

 

  parse_options (in_stream,temp_options, section_name,section_name);

 

  return  get_options();
 
}







//  *******************


void InputParser::reset_all_maps()
{
 
  atomistic_regions_map.clear();

  model_BC_map.clear();
  physical_model_map.clear();
  //blocks_map.clear();


  // prop_labels_map.clear();
  string_prop_labels_map.clear();
  //  vector_prop_labels_map.clear();
  
  vector_string_prop_labels_map.clear();
            
 
}







// private  method: utility  to  find a  keyword 
void InputParser::find_keyword(ifstream& in_stream, const std::string& keyword)

{

  bool found = false;
  std::string str,  label;


  assert(in_stream.good());
  in_stream >>  label;

  while (skip_comments(in_stream,label ) == true )
  {
    in_stream >> label; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
  } 



  cut_off_CR(label, in_stream);

  cut_off_comment(label, in_stream); //  case  Region#commmm
  
  while (  ( found == false) && (!in_stream.eof()) ) //   
  { //while
        
    if  (label == keyword  )
    {
      found = true;


     
      break;
    }  

    in_stream >>  label;

    while (skip_comments(in_stream,label ) == true )
    {
      in_stream >> label; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
    } 

    cut_off_comment(label, in_stream);
              
  }   


  if (found == false)

  {





    ostringstream os;
    os << "keyword "<< keyword <<  "  not  found in input file!";
    throw InitFailedException(os.str()); 

  }

}

//**************************************************


  // private  method: utility  to  find an optional  keyword 
bool InputParser::find_optional_keyword(ifstream& in_stream, const std::string& keyword)

{

  bool found = false;
  std::string str,  label;


  assert(in_stream.good());
  in_stream >>  label;

  while (skip_comments(in_stream,label ) == true )
  {
    in_stream >> label; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
  } 

  cut_off_comment(label, in_stream); //  case  Region#commmm

  
  while (  ( found == false) && (!in_stream.eof()) ) //   
  { //while
        
    if  (label == keyword  )
    {
      found = true;
      break;
    }  

    in_stream >>  label;

    while (skip_comments(in_stream,label ) == true )
    {
      in_stream >> label; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
    } 

    cut_off_comment(label, in_stream);
              
  }   

  return found;


  //  if (found == false)

  //  {

  //    return found;



  //   cerr << " Error: keyword "<< keyword <<  "  not  found !  " << endl ;
  //   throw InitFailedException("ERROR "); 
  //   //  exit(1);

  //  }



}




//**************************************************




 
const multimap <const string, ModelStructure*>& InputParser::read_models()

{

  string section_name;
  section_name = "$Models";

  std::ifstream in_stream (filename.c_str()) ;

  assert(in_stream.good());

  find_keyword(in_stream,section_name );  //  looks  for "$Models$"; 

  parse_model(in_stream);
 
  in_stream.close();

  return model_structure_map;


}





//   *************************************************************
//  NEW  PARSE  OPTIONS  WHICH   READS  ALL  STRINGS:
//  ALL  THE  INPUT  ITEMS  ARE  READ  AS  STRINGS AND  
//  PUT   IN A  ModelOptions  object   
//   *************************************************************


// void InputParser::parse_options(ifstream& in_stream, ModelOptions& region_options )
void InputParser::parse_options(ifstream& in_stream, ModelOptions& region_options ,
                                const std::string& block, const std::string& section  )

{

  // *******************************************************************************************************
  //  ***** IDEA:  label = value  -> map <string,string> 
  // 
  //  
  // *******************************************************************************************************


 
  string name, str, mat_name ;
  name = "Start";
  int id;
  bool block_ended;

  block_ended = false;


  // ****************************
  //  assigned   values (local)
  // ****************************
  vector<double> v;              //  scalar  and   vector  values
  vector<int> v_int; 
  //  vector<double>  vect;
  vector<string>  v_string;
  vector<string>  v_list_string;   //   value of  list_string
  vector<string> block_termination;
  block_termination.clear();

  v_string.clear();
  //  vect.clear();
  v.clear();
  v_int.clear();
  // ******************


  //------------------------
  //  labels of  properties
  //------------------------


  vector<string> v_label;
  vector<string>vect_label;
  vector<string> v_label_string;  //  label  for  string  value
  vector<string>list_label; //  label  for  list_string  value

  v_label_string.clear();
  vect_label.clear();
  v_label.clear();
  list_label.clear();

  //  vect.push_back(0);

   
  // *************************************************************************************************
  //  rule<phrase_scanner_t>  for  phrase level parsing (e.g. with a  separator character (comma) )
  // rule <> for  space  separated  list (character level parsing) : NEEDS  TO READ  EXPLICITLY SPACES !!!
  // *************************************************************************************************




  //   rule<>list_of_numbers_space_sep = ch_p('(')>> *(space_p) >>*real_p[push_back_a(v)] >>
  //     *( *(space_p) >> real_p[push_back_a(v)])>> *(space_p)>> ch_p(')')  ;

 


  // *******  special  characters for  the  label
  rule<>special_char =  (ch_p('_') | ch_p('-') |  ch_p('.') |  ch_p('/')   |  ch_p('+') | ch_p(',') 
                         | ch_p('%')   | ch_p('@') | ch_p('[') | ch_p(']') );

  rule<>dot = ch_p('.');

  //    rule<>special_char =  (ch_p('_') | ch_p('-') |  ch_p('.') |  ch_p('/')   |  ch_p('+') | ch_p(',')    );
  //  rule<>special_char =  (ch_p('_') | ch_p('-') |  ch_p('.') |  ch_p('/')   |  ch_p('+')    );

  // ******* separator for  list  entry: 'property = ( x , y, z )'
  rule<>separator = ch_p(',');

  // rule<>label  = (+alnum_p)>>   * ( ch_p('_') >> *(+alnum_p) ) ;

  //  rule<>label  = (+alnum_p)>>   * ( (special_char ) >> *(+alnum_p) ) ;  //2006


  // ******** label = name of  property

  rule<>label  = *(special_char)>>  (+alnum_p)>>   * ( (special_char ) >> *(+alnum_p) ) ;

  //  rule<>list_string =  ch_p('(') >> *(space_p) >> (label)>>  *( *(space_p) >> (label) ) >> ch_p(')');
  //   rule<>list_string =  ch_p('(') >> *(space_p) >> *(separator) >> *(space_p) >> 
  //     (label) >> *(space_p) >>   *( *(space_p) >>  *(space_p) >> *(separator) 
  //                                   >> *(space_p) >>(label) ) >> *(space_p) >> ch_p(')');

  // ********* list  entry: 'property = ( x , y, z )'
  rule<>list_string =  ch_p('(') >> *(space_p) >> (label) >> *(space_p) >> 
    *(  *(space_p) >> *(separator) >> *(space_p) >>(label) ) >> *(space_p) >> ch_p(')');


  //  rule<>list_string =  ch_p('(') >> *(space_p) >> (label)>> >> *(space_p)>> 
  //  *( *(space_p) >> (label) >>*(space_p)  ) >> ch_p(')');

  //  rule<phrase_scanner_t>list_string =  ch_p('(') >> (label)>> *(',' >> label  )>>   ch_p(')');
  //  cannot  be  used  because other rules are  not phrase_scanner_t !!! (?) 


  // rule<>tag_value = if_p('(')[(list_string) ].else_p[label];
  // **********  property value  can  be  a  single string (label) or  a  list ( x , y, z )

  //  
  //  rule<>tag_value =  list_string | label;

  // for  the  case "searchpath =  ."  (path = local dir)
  rule<>extended_label = label | dot;

  rule<>tag_value =  list_string | extended_label;

  // ********************************************  NEW 30.11.06  *******************
  //tag_value

  // **********  general assignement 'property = value'
  rule<>assignement_string  =  (label)[push_back_a(v_label_string)] >> 
    *(space_p) >> ch_p('=')>> *(space_p) >> ((tag_value)[push_back_a(v_string)])   ;

  // ***********************************************************************




  // ******************************************************************* 
  //list of  only  strings !!
  // 

  // **********  in  general  on a  line there can be  several assignements 
  rule<>list_of_assignement =  (  assignement_string) >>   //  OK  23.11.06  !!!!!
    *( *(space_p) >> (assignement_string )) ;


  // ************************************************************************

  // rule<> r_command  = *(space_p) >> (list_of_assignement | assignement_vector)    >> *(anychar_p);


  // **********  highest  level  rule for  parsing
  rule<> r_command  = *(space_p) >> (list_of_assignement )    >> 
    *(space_p) >> !(comment_p("#") >> *(anychar_p)) >>*(space_p)  ;
  // *(anychar_p);  

  //  possibly COMMENTs after list_of_assignement !


  //  ************   highest  level  rule for  parsing when first    line begins with "{par =..."  
  // case {par = pippo ....  !
  //
  // if the first    line begins with "{par =..."  
  rule<> r_command_started  = *(space_p) >> (ch_p("{")) >> (*(space_p)) >> (list_of_assignement )    >> 
    *(space_p) >> !(comment_p("#") >> *(anychar_p)) >>*(space_p)  ;


  //  ************  highest  level  rule for  parsing when first    line is {par = pippo }
  // case {par = pippo }  ....  !
  //if the first    line is all:  "{par =...}"
  rule<> r_command_started_and_terminated  = *(space_p) >> (ch_p("{")) >> (*(space_p)) >> 
    (list_of_assignement )    >> *(space_p) >> (ch_p("}")) >> 
    *(anychar_p) >>*(space_p)  ;
  //   !(comment_p("#") >> *(anychar_p)) >>*(space_p)  ;




  //  ******** reads a  line  with  only start symbol "{    "
  rule<> r_start_symbol  = *(space_p) >> (ch_p("{")) >>   
    *(space_p)>> *(comment_p("#")) ;
  //>> !(comment_p("#") >> *(anychar_p)) >>*(space_p)  ;


  rule<> r_empty_line = *(space_p);

  rule<> r_empty_block = *(space_p) >> (ch_p("{")) >> (*(space_p)) >> (ch_p("}")) >>
    !(comment_p("#") >> *(anychar_p)) >>*(space_p);


  // ********* rule for  a   line   terminated with  }
  rule<> r_command_terminated  = *(space_p) >> (list_of_assignement )    >> 
    *(space_p) >> (ch_p("}")) >>  *(space_p) >> !(comment_p("#") >> *(anychar_p)) >>*(space_p)  ;


  //  ****************************************************************

  //  FIRST LINE !!!
  // 1) case "{         "
  //  2)  case  "{par = pippo ...."
  // 
  //

  int line_counter ;
  line_counter = 1;


  while ( getline(in_stream, str) )
  { //while getline
    if   (!(parse(str.c_str(), comment_p("#")   , space_p).full) ) 
    {// if !comment_p("#") 

      if (!(parse(str.c_str(), r_empty_line ).full) )
      { // if !empty line
   
        if  ((parse(str.c_str(), r_start_symbol  ).full) ) 
        { // if parse   {

          //  OK  
          break;
        }

        //   first line with "{par = "    !!!
        else if (  parse(str.c_str(), r_command_started ).full )  
        { // 


          if ( !(v_string.empty()) )
          {

            for (int i =0; i< v_label_string.size();++i)
            {
              string_prop_labels.push_back(v_label_string[i]);
              //  cout <<  v_label_string[i]<<  "    " <<  v_string[i] <<  endl ;

              string_prop_labels_map.insert(make_pair(v_label_string[i], v_string[i]) );

              //  put option pair in ModelOptions  object
              region_options.set_option(v_label_string[i],v_string[i]) ;

            }
          }
     

          v_label_string.clear();
          vect_label.clear();
          v_label.clear();
          v_string.clear();
          v.clear();
          v_int.clear();
          break;
    

        }

        //   first line is all! : termination condition  with  } on  the  line  !!!
        else if (  parse(str.c_str(), r_command_started_and_terminated).full )           //  not skipping spaces


        { // if parse


          if ( !(v_string.empty()) )
          {

            for (int i =0; i< v_label_string.size();++i)
            {
              string_prop_labels.push_back(v_label_string[i]);
              //  cout <<  v_label_string[i]<<  "    " <<  v_string[i] <<  endl ;

              string_prop_labels_map.insert(make_pair(v_label_string[i], v_string[i]) );

              //  put option pair in ModelOptions  object
              region_options.set_option(v_label_string[i],v_string[i]) ;

            }
          }
	     

          // cout << endl ;
          //     
          block_ended =  true;  //  flag     
          break;

    

        }

        //  empty block :  "{}"  or   "{  }" 
        else if  ((parse(str.c_str(),r_empty_block  ).full) ) 
        { // if parse   {

          //  OK  
          block_ended =  true;  //  flag   
          break;
        }

  

        else  

        {

          std::ostringstream stm;

          stm <<  "Syntax error in input file on line " << line_counter;
          if (block.size() > 0)
            stm << " of block " << block;
          stm << " in section " << section  << endl;

       
          throw InitFailedException(stm.str());

        }

       
      }

    }

    line_counter++ ; 

  }


  // *******************************************************************
  //  IF  THERE ARE OTHER LINES IN  THE  BLOCK (BLOCK NOT ENDED)
  //  *******************************************************************



  if (! (block_ended) )

  { // ! (block_ended)

    while ( getline(in_stream, str) )
    {  // while getline 2

   
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

	     

     


          if ( !(v_string.empty()) )
          {

            for (int i =0; i< v_label_string.size();++i)
            {
              string_prop_labels.push_back(v_label_string[i]);
              //    cout <<  "v_label_string[i] ****  " <<  v_label_string[i]<<  "    " <<  v_string[i] <<  endl ;

              string_prop_labels_map.insert(make_pair(v_label_string[i], v_string[i]) );

              //region_options.set_option(v_label_string[i],v_string[i])) ;
              region_options.set_option(v_label_string[i],v_string[i]) ;


                   

            }
          }
	     

       

	      

          v_label_string.clear();
          vect_label.clear();
          v_label.clear();
          v_string.clear();
          v.clear();
          v_int.clear();

  
	    

        }

     

        //  *******************************************************

        //   termination condition  !!!
        else if  (parse(str.c_str(), ( *(space_p)>>ch_p("}")  >> *(anychar_p) >>*(space_p) )  , space_p).full) 
        {  
          //  cout << endl ;
          //         cout <<  "Fine  !" << endl ;
          //   if (name == "a")
          break;
        }


         



        //   termination condition  with  } on  the  line  !!!
        else if (  parse(str.c_str(), r_command_terminated ).full )           //  not skipping spaces


        { // if parse


          if ( !(v_string.empty()) )
          {

            for (int i =0; i< v_label_string.size();++i)
            {
              string_prop_labels.push_back(v_label_string[i]);
              //  cout <<  v_label_string[i]<<  "    " <<  v_string[i] <<  endl ;

              string_prop_labels_map.insert(make_pair(v_label_string[i], v_string[i]) );

              //  put option pair in ModelOptions  object
              region_options.set_option(v_label_string[i],v_string[i]) ;

            }
          }
	     

          // cout << endl ;
          //         cout <<  "Fine device !!! !" << endl ;
          break;

    

        }





   



        //     skip   all  the  other  kinds  of  lines  !!!!
        else  if (parse(str.c_str(), if_p("{")[(+alpha_p) 
                                               ] , space_p ).full)
        {  
          //   cout << "SKIP" << endl ;
          //   if (name == "a")
          //   break;
        }


        else  

          //  throw InitFailedException("....................");


        {

          std::ostringstream stm;

          stm << "Syntax error in input file on line " << line_counter;
          if (block.size() > 0)
            stm << " of block " << block;
          stm << " (section " << section  << ")" << endl;
          stm <<  "Correct syntax is : 'label' = 'value' "
            "'label' = 'value' .......# 'comment' " <<
            endl << "A comment line  must be preceded by '#' "<< endl;
         
          throw InitFailedException(stm.str());

        }

      }


      v_label_string.clear();
      vect_label.clear();
      v_label.clear();
      v_string.clear();
      v.clear();
      v_int.clear();
     
 
      line_counter++ ; 

    }  //  end  while

  }

}




// public  method to  read  device regions

//const map <ID, RegionStructure>& InputParser::read_device(void)
void InputParser::read_device(void)

{

  std::string  label, keyword,region_name, atomistic_region_name, section_name   ;
  std::string cluster_name,device_block_name   ;

  std::ifstream in_stream (filename.c_str()) ;
  ID current_region_ID; //  to be  deleted
  bool  check_error;
  ID region_counter,atomistic_region_counter, cluster_counter ,interface_counter ;
  region_counter = 0;
  cluster_counter = 0;
  interface_counter = 0;

  atomistic_region_counter = 0;

  reset_all_maps();

  assert(in_stream.good());

  section_name =  "Device";

  section_name = "$"+section_name;

  find_keyword( in_stream,section_name  );
 
  // { 
 
  // while (  ( label  != end_symbol) && (!in_stream.eof()) ) //   
  //        {   

 
  check_error = skip_to_bracket(in_stream);  // go on reading until the  char  '{' is  read  
  if (check_error == true)
    throw InitFailedException("\'Device\' block missing in input file."); 


  //  in_stream >>  label; //  read   start_symb



  //  ********************************
  //  read keyword Region  (or Atomistic )
  in_stream >> keyword;

  //  skip_comments(in_stream,keyword );

  //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
  while (skip_comments(in_stream,keyword ) == true )
  {
    in_stream >> keyword; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
  } 

  cut_off_comment(keyword, in_stream); //  case  Region#commmm

  //  cout <<  " region_k = " << keyword <<  endl; 

 

  while  (keyword !=  end_symb)
  {


    if  (keyword == "Atomistic")
    {
      //      *********** read  atomistic ***********************

      //     read   region_name
      in_stream >>atomistic_region_name;
      atomistic_region_counter++;
      //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
      while (skip_comments(in_stream,atomistic_region_name) == true )
      {
        in_stream >> atomistic_region_name; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 

      cut_off_comment(atomistic_region_name, in_stream );  //  in  case  layer#commmm

      //cerr << " ************** ATOMISTIC ************ " <<  atomistic_region_name << endl;

      temp_options.clear();
      //   parse_options(in_stream,temp_options  ); //    read  the  block  between  { and  }


      parse_options(in_stream,temp_options, keyword, section_name);

      //....................

      // create new RegionStructure
      RegionStructure current_region_structure;

      string  atomistic_region_numb  , def1;

      // put  atomistic region_name  in RegionStructure
      current_region_structure.set_region_name(atomistic_region_name );

      current_region_structure.set_model_options(temp_options);
      //  temp_options contains  atomistic ID and all  other data (list of phis region ID) 

      //     atomistic_region_numb = temp_options.get_option( "atomistic_region_numb" ,def1);

      atomistic_regions_map.insert(make_pair (atomistic_region_counter, current_region_structure ));




    }


    else if (keyword == "Region")
    {  //  read  Physical Region

      //    if  (keyword != "Region")
      //     {
      //       throw InitFailedException("SYNTAX ERROR in input  file (device section): keyword Region is  missing! ");
    
      //     }

      //     read   region_name
      in_stream >>region_name;
      region_counter++;

      //  cout <<  " region_name = " << region_name <<  endl; 
    
      //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
      while (skip_comments(in_stream,region_name) == true )
      {
        in_stream >> region_name; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 

      cut_off_comment(region_name, in_stream );  //  in  case  layer#commmm
  
       
 
      //  parse_options(in_stream);   //    read  the  block  between  { and  }



      //  ********************  NEW options  **********************

   


      // create  new ModelOptions  
      // ModelOptions temp_region_options; private  member 
      //  

      temp_options.clear();
      //    parse_options(in_stream,temp_options  ); //    read  the  block  between  { and  }

      // region_name
      string block_name ;
      block_name = keyword + " " + region_name;

      parse_options(in_stream,temp_options,block_name, section_name   ); 
      // extract and delete name, ID and  mat   from  temp_region_options
      // ................
      string material, region_numb , def;

      def = "";

      //  temp for  back-compatibility
      //
      material  = temp_options.get_option( "mat" ,def);
      temp_options.delete_option ( "mat" );
      //      cerr << "****************** mat =  " << material << " ******************"<<  endl;


      if (material == "")
      {
        material  = temp_options.get_option( "material" ,def);
        temp_options.delete_option ( "material" );
        //        cerr << "****************** material !! =  " << material << " ******************"<<  endl;

      }

      region_numb = temp_options.get_option( "reg_numb" ,def);
      temp_options.delete_option ( "reg_numb" );
      //      cerr << "****************** reg_numb =  " << region_numb << " ******************"<<  endl;



      if ( region_numb == "")
      {
        region_numb = temp_options.get_option( "mesh_regions" ,def);
        temp_options.delete_option ( "mesh_regions" );
        //      cerr << "****************** mesh_regions !!! =  " << region_numb << " ******************"<<  endl;

      }




    


      // create new RegionStructure
      RegionStructure current_region_structure;


      //put  region_name  current_region_ID, material  and temp_region_options in RegionStructure

      current_region_structure.set_region_name(region_name);
      current_region_structure.set_region_ID(region_numb);

      current_region_structure.set_material_name(material);
      current_region_structure.set_model_options(temp_options);

      //   current_region_structure.set_model_options(temp_region_options);
      //    current_region_ID   = atoi(region_numb.c_str());   OBSOLETE ????

      // ****************************************************
      // instead of  current_region_ID,  put  incremental ID in map device_map: region_counter
      // device_map.insert(make_pair (region_counter, current_region_structure )); 

      device_map.insert(make_pair (region_counter, current_region_structure ));

      //    device_map.insert(make_pair (current_region_ID, current_region_structure )); 

      // ***************************************


      //  put region_name in map  map_region <string,string >
      // region_map.insert(make_pair ("region_name",region_name  ));  
      //   string_prop_labels_map.insert(make_pair ("region_name",region_name  )); //  Obsolete ???


      //  get   ID of   the  current region :
      //    ID = map_region[reg_mumb]
      //  current_region_ID   = atoi(string_prop_labels_map["reg_numb"].c_str());
      // cout <<  current_region_ID << endl;


      //  put  in  map  <ID, map_region   >
      //   the map map_region <string,string >
      //  device_map.insert(make_pair (ID,map_region  ));   

      //    device_map.insert(make_pair (current_region_ID, string_prop_labels_map ));  


      //cout <<  "(device_map[current_region_ID])[mat]  =  " << (device_map[current_region_ID])["mat"] << endl  ;

     
      //     string_prop_labels_map.clear(); //  Obsolete ???


    } //  end  read Region




    else if (keyword == "Cluster")
    {  //  read  Cluster

      //  Cluster section contains a list of  mesh_regions :
      //  e.g., mesh_regions = ( 2, 5,6) ;  these  regions  CAN have  different materials associated !!!
      // 
      //  clusters will be  put  in  a  separated map  cluster_map 
      //  to  avoid  confusion with  Regions of the device

      string  region_numb , def;

      //     read cluster_name
      in_stream >>cluster_name;
      //   region_counter++;
      cluster_counter++;

      //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
      while (skip_comments(in_stream,cluster_name) == true )
      {
        in_stream >> cluster_name; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 

      cut_off_comment(cluster_name, in_stream );  //  in  case  layer#commmm

      temp_options.clear();
      //   parse_options(in_stream,temp_options  ); //    read  the  block  between  { and  }
      parse_options(in_stream,temp_options,keyword, section_name    ); 

      // get mesh_region(s)
      region_numb = temp_options.get_option( "mesh_regions" ,def);

      // create new RegionStructure  with  empty fields (material, region_ID, material_name)
      RegionStructure current_region_structure;

      // put mesh region(s) associated to cluster
      current_region_structure.set_region_ID(region_numb);

      // put options of cluster 
      current_region_structure.set_model_options(temp_options);

      // put name of the cluster
      current_region_structure.set_region_name(cluster_name);

      // ***NO !     // insert in  map of the device as made with the simple Region !!!
      //      device_map.insert(make_pair (region_counter, current_region_structure ));

      // insert in  map cluster_map  
      cluster_map.insert(make_pair (cluster_counter, current_region_structure ));

    }// end read cluster


    // *******************  Interface *************************


    else if (keyword == "Interface")
    {  //  read  Interface  (generic for boundary  conditions, interface  models,  ecc.)
      //  interfaces  will be  put  in  a  separated map  interface_map 

      //     read interface_name
      in_stream >>device_block_name;
      //   region_counter++;
      //  cluster_counter++;
      interface_counter++;

      //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
      while (skip_comments(in_stream,device_block_name) == true )
      {
        in_stream >> device_block_name ; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 

      cut_off_comment(device_block_name , in_stream );  //  in  case  layer#commmm

      temp_options.clear();
      //    parse_options(in_stream,temp_options  ); //    read  the  block  between  { and  }
      parse_options(in_stream,temp_options,keyword, section_name   ); 

 
      // create new RegionStructure  with  empty fields (material, region_ID, material_name)
      RegionStructure current_region_structure;

      // put options of interface
      current_region_structure.set_model_options(temp_options);

      // put name of the interface 
      current_region_structure.set_region_name(device_block_name);

      // insert in  map cluster_map  
      interface_map.insert(make_pair (interface_counter, current_region_structure ));



    }// end read Interface


    else if (keyword == "Options")
    {

      // read Options block  for  General Options common  to  ALL  the  device  regions
      // 
      temp_options.clear();
      parse_options(in_stream,temp_options, keyword ,  section_name   );

      set_device_options(temp_options);


    }// end read General Options



    else 
    {
      throw InitFailedException("In input file (\'Device\' section): "
                                "keyword \'Region\' or \'Atomistic\'  is  missing! ");
    
    }


    //    next   Region  
    in_stream >> keyword;

    //   skip_comments(in_stream,keyword );

    while (skip_comments(in_stream,keyword ) == true )
    {
      in_stream >> keyword; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
    } 

    cut_off_comment(keyword, in_stream); //  case  Region#commmm


  }  //end   while 
     
  // map <ID, RegionStructure>& InputParser::get_device_map(void)
  //    return device_map;  //  to  be  read  from  get_device_map  !!!!


}      //  end   Device   section 





map <ID, RegionStructure>& InputParser::get_device_map(void) 

{

  return device_map;

} 

map <ID, RegionStructure>& InputParser::get_cluster_map(void) 

{

  return cluster_map ;

} 


map <ID, RegionStructure>& InputParser::get_interface_map(void) 

{

  return interface_map ;

} 





map <ID, RegionStructure>& InputParser::get_atomistic_map(void)
{

  return atomistic_regions_map;

} 



// private  method: utility  to  skip  comments (everything on a line, after "#" ) 
bool InputParser::skip_comments(ifstream& in_stream, const std::string& item)

{

  bool  skip  = false;

  if (item == "#" | (std::strncmp ((item.c_str()), "#",1) == 0) )
  {
    in_stream.ignore(256,'\n');  //  if  the   read keyword is # or  begins with #: ignore  all  the  line 
    //   in_stream >> item;  //      and   read  the  next keyword !!! 

    skip  = true;

  } 

  return  skip;


} 



//  ***************************************************

//    parse  model  section :  for  each model a   list  of  phys regions and a list of BC definitions
//   are  read from input  file
void
InputParser::parse_model(ifstream& in_stream)
{// method  for   parsing  of   model section

  ID BC_region_counter ;
  bool  check_error;
  string BC_regions_keyword_string , keyword_BC_Region_string, model_keyword_string,
    phys_model_label, options_label    ;

  //  numb_regions_keyword_string = "numb_regions";
  //  phys_regions_keyword_string = "phys_regions";
  BC_regions_keyword_string = "BC_Regions";
  keyword_BC_Region_string = "BC_Region";
  model_keyword_string  =  "model";
  phys_model_label = "physical_model";
  //  simulation_name_label = "name";
  //  physical_regions_label = "phys_regions";
  options_label   = "options";

  set<string>  model_section_keywords;
  // set<string>::iterator   s_it;

  model_section_keywords.insert("BC_Regions");
  model_section_keywords.insert("model");
  model_section_keywords.insert("physical_model");
  model_section_keywords.insert("options");
  model_section_keywords.insert("}");




  //  vector <ID> list_physical_regions  ;

  //  vector <string> list_physical_regions  ;
  string list_physical_regions  ;

  //  vector <string> model_list;  ->  private  member

  string model_name, label , keyword,physical_model_name,simulation_name ;
  string  item, model_keyword, numb_regions_keyword, equal,regions_keyword ,BC_regions_keyword;
  //    start_symbol,end_symbol;
  ID numb_regions ;
  ID reg_numb;
 
  ID current_BC_region_ID ;

  string dummy;
  char c;

  string  start_symbol,end_symbol, BC_region_name ;

  string str;

 
  check_error = skip_to_bracket(in_stream);  // go on reading until the  char  '{' is  read  
  if (check_error == true)
    throw InitFailedException("In input file: \'Models\' section  missing."); 
  //   cerr <<  "ERROR "  << endl; 


  
  //   find_keyword(in_stream,start_symb ); //  looks  for  "{"

  in_stream >> model_keyword;
  //  cout << model_keyword<< endl;



  while (skip_comments(in_stream,model_keyword ) == true )
  {
    in_stream >> model_keyword; // if  the  whole  line has  ben  skipped:
    // read  the  next keyword !!! 
  } 

  if  ( model_keyword != model_keyword_string)
  {
    throw InitFailedException("In input file: model keyword missing "
                              "in \'Models\' section");  

   
  }


    
  //while (  (model_keyword != "End") && (!in_stream.eof()) ) 
  while (  (model_keyword != end_symb ) && (!in_stream.eof()) ) 
  { //  while loop  models

    // Model n
    if  (model_keyword != model_keyword_string)
    {

      throw InitFailedException("In input file: model keyword missing "
                                "in \'Models\' section");  
   
    }




    in_stream >>  model_name ;
  



    while (skip_comments(in_stream, model_name ) == true )
    {
      in_stream >>  model_name ; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
    } 
    //**********************************
        //cut_off_comment(string& label )//   erase  possible comments attached to model_name

        cut_off_comment(model_name, in_stream);


    reset_all_maps();  //  clear  the  IP  maps !!!!
    model_BC_map.clear();

    //   new   ModelStructure  object  !!!!
    current_model_point = new ModelStructure(model_name);
    //  and ModelStructure.set_model_name !!!!!!!!!!

   

    check_error = skip_to_bracket(in_stream);  // go on reading until the  char  '{' is  read  
    if (check_error == true)
    {

      std::ostringstream stm;
      stm << "In input file: block missing for model " << model_name <<  endl;
      throw InitFailedException(stm.str());


    }

    //  ***********************************************************************

    //  NEW :  
    //  read  optional label  

    // **************************************************************************************
    //  Read  the  next item: it  can  be :
    //  1) a  closing  bracket ->  go  to  next  model !
    //  2) one  of  the  keywords  in  the  set  for  model section (model_section_keywords)
    //  ->  read related section
    //  3)  an  unknown   label  ->  throw  exception !
    //
    // *************************************************************************************




    in_stream >>  label ;
    while (skip_comments(in_stream,  label ) == true )
    {
      in_stream >> label  ; // if  the  whole  line has
      //  been  skipped: read  the next keyword !!! 
    } 

    cut_off_comment(label, in_stream); //  case  options#commmm vv

    // bool InputParser::check_label(set<string>& section_keywords, const string& label )

    if (!check_label(model_section_keywords, label))

    {

      std::ostringstream stm;
      stm << "In input file: unknown keyword \'"
          << label << "\' in models section, model "
          << model_name <<  endl;
      throw InitFailedException(stm.str());

    }

    // label is one  of  the  keywords  in  the  set  for  model section!!
  
    //  read  optional label  "options"
    if  ( label  == options_label )
    {
               
      temp_options.clear();
      //    parse_options(in_stream,temp_options);
      parse_options(in_stream,temp_options, label, model_name   );

      current_model_point->set_model_options(temp_options );

      // after "options"  block,   read  the  next  label

      in_stream >>  label ;
      while (skip_comments(in_stream,  label ) == true )
      {
        in_stream >> label  ; // if  the  whole  line has
        //  been  skipped: read  the next keyword !!! 
      } 

      if (!check_label(model_section_keywords, label))
      {

        std::ostringstream stm;
        stm << "In input file: unknown keyword \'"
            << label << "\' in models section, model "
            << model_name <<  endl;
        throw InitFailedException(stm.str());

      }




    }


    //  NEW ***********************
    //  READ  OPTIONAL PHYSICAL_MODEL AND  OPTIONAL  BC_REGIONS

    //



  



    //  READ  OPTIONAL PHYSICAL_MODEL SECTION(S)

    while (  (label == phys_model_label  ) && (!in_stream.eof()) ) 
    
    {  //  while loop  physical models

      // phisical Model n
               

      in_stream >>  physical_model_name ;

      while (skip_comments(in_stream, physical_model_name ) == true )
      {
        in_stream >>  physical_model_name ; // if  the  whole  line has  ben  skipped: 
        //read  the  next keyword !!! 
      } 

      cut_off_comment(physical_model_name, in_stream); //  case  recombination#commmm

      //     reset_all_maps();  //  clear  the  IP  maps !!!!




           

      temp_options.clear();

      string_prop_labels_map.clear();  //  obsolete !!!!
   
      //    parse_options(in_stream);   //    read  the  block  between  { and  }
      //    parse_options(in_stream,temp_options  );

      parse_options(in_stream,temp_options,physical_model_name, model_name  );

      //parse_options(in_stream);   //    read  the  block  between  { and  }
      // 
      // 
      //   current_model_point->set_phys_model_map( string_prop_labels_map);

      physical_model_map.insert(make_pair (physical_model_name,temp_options));
      //    physical_model_map[physical_model_name]=temp_options;



      //  current_model_point->set_phys_model_options(  temp_options );// old

      current_model_point->set_physical_model_map(physical_model_map);


      // void set_phys_model_options(  ModelOptions& physical_model_options  );



      string_prop_labels_map.clear();   //  obsolete !!!!  
 

      in_stream >>  label; 
   
      while (skip_comments(in_stream, label ) == true )
      {
        in_stream >> label ; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 
  
      cut_off_comment(label, in_stream); //  case  Recomb#commmm


      if (!check_label(model_section_keywords, label))
      {

        std::ostringstream stm;
        stm << "In input file: unknown keyword \'"
            << label << "\' in models section, model "
            << model_name <<  endl;
        throw InitFailedException(stm.str());

      }




    }// end while loop  physical models


    //     //   ----------------------------------------------------------------



    //----------------------------------------
    //  NOW READ  OPTIONAL BC_REGIONS SECTION(S)
    // ---------------------------------------

    if  ( label  == BC_regions_keyword_string)
    { //  IF label  == BC_regions_keyword_string************

      check_error = skip_to_bracket(in_stream);  // go on reading until the  char  '{' is  read  
      if (check_error == true)
      {
        std::ostringstream stm;
        stm << "In input file: BC regions block missing for model "
            << model_name <<  endl;
        throw InitFailedException(stm.str());   
        
      }

      BC_region_counter = 0 ;//  reset  counter for  boundary regions for this model
      //  one single   BC region can refer to more than one BC_reg_ID !!!!
 
      //  read keyword BC_Region
      in_stream >> keyword;
      //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
      while (skip_comments(in_stream,keyword ) == true )
      {
        in_stream >> keyword; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 

      cut_off_comment(keyword, in_stream); //  case  BC_Region#commmm


      while  (keyword !=  end_symb)   //    loop  BC_region   block
      {//    loop  BC_region   block

        if  (keyword != keyword_BC_Region_string)
        {
          std::ostringstream stm;
          stm << "In input file: missing keyword BC_Region in model "
              << model_name <<  endl;
          throw InitFailedException(stm.str());   
        

          //    cerr << " SYNTAX ERROR in input file: missing keyword BC_Region in model  " << model_name <<  endl;
          //    throw InitFailedException("SYNTAX ERROR ");

        }

        //     read   BC_region_name
        in_stream >>BC_region_name;
        BC_region_counter++; // ***********************


        //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
        while (skip_comments(in_stream,BC_region_name) == true )
        {
          in_stream >> BC_region_name; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
        } 

        cut_off_comment(BC_region_name, in_stream); //  case  cathode#commmm

        //*************************************          
            // ModelOptions temp_BC_options;  private  member 
            //*************************************

                temp_options.clear();
     
        //    parse_options(in_stream);   //    read  the  block  between  { and  }
        //   parse_options(in_stream,temp_options);
        parse_options(in_stream,temp_options, keyword , model_name );

        string  BC_region_numb , def;
        BC_region_numb = temp_options.get_option( "BC_reg_numb" ,def);

        // create new RegionStructure
        RegionStructure current_BC_region_structure;

        current_BC_region_structure.set_region_name(BC_region_name);
        current_BC_region_structure.set_region_ID(BC_region_numb);

        current_BC_region_structure.set_material_name("");
        current_BC_region_structure.set_model_options(temp_options);

        current_BC_region_ID   = atoi(BC_region_numb.c_str());

        // ****************************************************
        // instead of  current_BC_region_ID,  put  incremental ID in map model_BC_map: BC_region_counter
        // model_BC_map.insert(make_pair (BC_region_counter,current_BC_region_structure  )); 

        model_BC_map.insert(make_pair (BC_region_counter,current_BC_region_structure  )); 

        //        model_BC_map.insert(make_pair (current_BC_region_ID,current_BC_region_structure  ));  

        // *******************************  END  2.2.07/1.3.07 *****************

        //  put region_name in map  map_region <string,string >
   
        string_prop_labels_map.insert(make_pair ("BC_region_name",BC_region_name  )); 

        //  get   ID of   the  current BC_region :
        //    ID = map_region[BC_reg_mumb]
        //      current_BC_region_ID   = atoi(string_prop_labels_map["BC_reg_numb"].c_str());
        //  cout <<  current_BC_region_ID << endl;

        //  put  in  map  <ID, map_region   >
        //   the map map_region <string,string >
        //  device_map.insert(make_pair (ID,map_region  ));   

        //      model_BC_map.insert(make_pair (current_BC_region_ID, string_prop_labels_map ));  
        //  cout <<  "([ model_BC_map.current_BC_region_ID])[value]  =  " 
        //  << (model_BC_map[current_BC_region_ID])["value"] << endl  ;

        //    next   BC_Region    
        string_prop_labels_map.clear();
        in_stream >> keyword;

        //   skip_comments(in_stream,keyword );

        while (skip_comments(in_stream,keyword ) == true )
        {
          in_stream >> keyword; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
        } 

        cut_off_comment(keyword, in_stream); //  case  BC_Region#commmm



      }  //end   while  loop  BC_region   block

      //   END   READ  BC  REGIONS  



      //  ******************************
      //  PUT   model_BC_map IN   OBJECT  ModelStructure

      //  ModelStructure->set_model_BC_map(model_BC_map)  **********************************
      //  ModelStructure.set_model_BC_map( map  <unsigned int,  map <string,string> >& 
      // id_BC_regions_map   )   !!!!!!!!!
      current_model_point->set_model_BC_map( model_BC_map);



   
      //   in_stream >>  end_symbol;

      //keyword  !!!!
      end_symbol = keyword;   //   closing  bracket  of  BC_Regions  section !!

      while (skip_comments(in_stream, end_symbol) == true )
      {
        in_stream >> end_symbol; // if  the  whole  line has  been  skipped: read  the  next keyword !!! 
      } 


      

      //    if  (strncmp ((end_symbol.c_str()),"}",1) != 0)
      if  (std::strncmp ((end_symbol.c_str()),end_symb.c_str(),1) != 0)

      {

        std::ostringstream stm;
        stm << "In input file: bad BC regions block in model " << model_name <<  endl;
        throw InitFailedException(stm.str());  

      
      }


      //  read  next  item :

      in_stream >>  label; 
      // cout << model_keyword<< endl;

      while (skip_comments(in_stream, label) == true )
      {
        in_stream >>  label ; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 

      cut_off_comment(label, in_stream); //  case  model#commmm

      if (!check_label(model_section_keywords, label)) 
      {
        std::ostringstream stm;
        stm << "In input file: unknown keyword \'" << label
            << "\' in BC regions block of model "
            << model_name <<  endl;
        throw InitFailedException(stm.str());  

        //   cerr << " SYNTAX ERROR in input  file: unknown keyword in $Models section, in   " <<  endl;
        //   throw InitFailedException(" BC_Regions block ");
      }

    } //  end  of  if  BC_Regions_keyword

    // *****************************************************************************


    // *****************************************************
    // if label   == }  then   NEXT  MODEL  !!!
    // ******************************************************
  
    if  (std::strncmp ((label.c_str()),end_symb.c_str(),1) == 0)  // read  a  closing  bracket !!

    {

      //END  OF  MODEL
      //  end  of  the  current model  !!
      //  put ModelStructure in  map  <model_name, *ModelStructure>
      model_structure_map.insert(make_pair (model_name, current_model_point )); 

      in_stream >>  model_keyword;  //  read  next model model_keyword OR closing  bracket

      // read  next item
      while (skip_comments(in_stream, model_keyword ) == true )
      {
        in_stream >>  model_keyword ; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 

      cut_off_comment(model_keyword, in_stream); //  case  model#commmm

      if (!check_label(model_section_keywords, model_keyword))
      {
        std::ostringstream stm;
        stm << "In input file: unknown keyword \'" << model_keyword <<
          "\' in \'Models\' section (close to model " << 
          model_name << ")" <<  endl;
        throw InitFailedException(stm.str());  

      }

    }




    //  NEXT  MODEL; if  model_keyword =  closing  bracket } ->  END  !!!
	      
  }  // end while   while (  (model_keyword != end_symb ) && (!in_stream.eof()) ) ***********  end   Model
      
  
  // ***********  end   $Model section
 

} //  end  method   parse_model




  
//   ************************************************************
//  Model  section  with  "boundary model"   scheme  only  blocks and  subblocks
//

void
InputParser::NEW_parse_model(ifstream& in_stream)
{// method  for   parsing  of   model section

  bool  check_error;
  string model_name, model_keyword, label ,  model_keyword_string ,options_label, model_block_name   ;

  set<string>  model_section_keywords;
  // set<string>::iterator   s_it;

  //  model_section_keywords.insert("BC_Regions");
  model_section_keywords.insert("model");
  model_section_keywords.insert("physical_model");
  model_section_keywords.insert("options");
  model_section_keywords.insert("}");


  options_label   = "options";
  model_keyword_string  =  "model";

  check_error = skip_to_bracket(in_stream);  // go on reading until the  char  '{' is  read  
  if (check_error == true)
    throw InitFailedException("In input file: \'Models\' section  missing."); 
  

  in_stream >> model_keyword;

  while (skip_comments(in_stream,model_keyword ) == true )
  {
    in_stream >> model_keyword; // if  the  whole  line has  ben  skipped:
    // read  the  next keyword !!! 
  } 

  if  ( model_keyword != model_keyword_string)
  {
    throw InitFailedException("In input file: model keyword missing in models section");  

    
  }


  while (  (model_keyword != end_symb ) && (!in_stream.eof()) ) 
  { //  while loop  models

    // Model n
    if  (model_keyword != model_keyword_string)
    {

      throw InitFailedException("In input file: model keyword missing in models section");  
   
    }


    in_stream >>  model_name ;

    while (skip_comments(in_stream, model_name ) == true )
    {
      in_stream >>  model_name ; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
    } 
   
    //   erase  possible comments attached to model_name
    cut_off_comment(model_name, in_stream);

    reset_all_maps();  //  clear  the  IP  maps !!!!


    //   new   ModelStructure  object  !!!!
    current_model_point = new ModelStructure(model_name);

    check_error = skip_to_bracket(in_stream);  // go on reading until the  char  '{' is  read  
    if (check_error == true)
      throw InitFailedException("In input file: model block missing.");


    // ***********************************************************************

    //  NEW :  
    //  read  optional block-label  

    //  **************************************************************************************
    //  Read  the  next item: it  can  be :
    //  1) a  closing  bracket ->  go  to  next  model !
    //  2) one  of  the  keywords  in  the  set  for  model section (model_section_keywords)
    //  ->  parse the  block
    //  3)  an  unknown   label  ->  throw  exception !
    //
    //  *************************************************************************************






    in_stream >>  label ;
    while (skip_comments(in_stream,  label ) == true )
    {
      in_stream >> label  ; // if  the  whole  line has
      //                       been  skipped: read  the next keyword !!! 
    } 

    cut_off_comment(label, in_stream); //  case  options#commmm vv


    if (!check_label(model_section_keywords, label))
    {
      ostringstream os;
      os << "In input file: unknown keyword \'" << label
         << "\' in models section!";
      throw InitFailedException(os.str());
    }

    // label is one  of  the  keywords  in  the  set  for  model section!!



    // ********  loop  on  all  the  blocks of  this  model ************

    while (!(std::strncmp ((label.c_str()),end_symb.c_str(),1) == 0))  // dont read  a  closing  bracket !!
    {  //  loop while there  are  blocks in  model

      //  read  optional label  "options"
      if  ( label  == "options") // options_label )
      {

        temp_options.clear();
        //   parse_options(in_stream,temp_options);
        parse_options(in_stream,temp_options,label, model_name  );

        current_model_point->set_model_options(temp_options );

      }

      else if   ( label  == "physical_model") //
      {

        model_block_name = "";

        in_stream >>   model_block_name ;

        while (skip_comments(in_stream,  model_block_name ) == true )
        {
          in_stream >>   model_block_name ; // if  the  whole  line has  ben  skipped: 
          //read  the  next keyword !!! 
        } 

        cut_off_comment( model_block_name , in_stream); //  case  recombination#commmm

        temp_options.clear();

        //     parse_options(in_stream,temp_options  );
        parse_options(in_stream,temp_options,label,model_name);

        physical_model_map.insert(make_pair ( model_block_name,temp_options));


        //  solo  alla  fine  ???
        current_model_point->set_physical_model_map(physical_model_map);


      }


      else if   ( label  == "boundary_model") //
      {
        model_block_name = "";

        in_stream >>  model_block_name ;

        while (skip_comments(in_stream,  model_block_name) == true )
        {
          in_stream >>   model_block_name ; // if  the  whole  line has  ben  skipped: 
          //read  the  next keyword !!! 
        } 

        cut_off_comment( model_block_name, in_stream); //  case  recombination#commmm

        temp_options.clear();

        //  parse_options(in_stream,temp_options  );
        parse_options(in_stream,temp_options,label,label  );


        boundary_model_map.insert(make_pair ( model_block_name ,temp_options));

        current_model_point->set_boundary_model_map( boundary_model_map);

      }

      else if   ( label  == "simulation") //
      {

        // ??????????
      }



      // after block,   read  the  next  label

      in_stream >>  label ;
      while (skip_comments(in_stream,  label ) == true )
      {
        in_stream >> label  ; // if  the  whole  line has
        //  been  skipped: read  the next keyword !!! 
      } 

      if (!check_label(model_section_keywords, label))
      {
        ostringstream os;
        os << "In input file: unknown keyword \'" << label
           << "\' in models section!";
        throw InitFailedException(os.str());
      }



    }

    // *****************************************************
    // if label   == }  then   NEXT  MODEL  !!!
    // ******************************************************

    if (std::strncmp ((label.c_str()),end_symb.c_str(),1) == 0)  // read  a  closing  bracket !!

    {

      //END  OF  MODEL
      //  end  of  the  current model  !!
      //  put ModelStructure in  map  <model_name, *ModelStructure>
      model_structure_map.insert(make_pair (model_name, current_model_point )); 

      in_stream >>  model_keyword;  //  read  next model model_keyword OR closing  bracket

      // read  next item
      while (skip_comments(in_stream, model_keyword ) == true )
      {
        in_stream >>  model_keyword ; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 

      cut_off_comment(model_keyword, in_stream); //  case  model#commmm

      if (!check_label(model_section_keywords, model_keyword))
      {
        ostringstream os;
        os << "In input file: unknown keyword \'" << label
           << "\' in models section!";
        throw InitFailedException(os.str());
      }

    }



    
    //  NEXT  MODEL; if  model_keyword =  closing  bracket } ->  END  !!!

  }  // end  while  models










} //  END   NEW_parse_model


//  **********************************************************








multimap <const string, ModelStructure*>& InputParser:: get_model_structure_map(void)
{
  
  
  return model_structure_map;

} //  end  method





const  ModelOptions& InputParser::get_options(void)

{

  return temp_options;

} 




// // **********************

// //  method to parse  list  of  IDs 
// //  NEW:  returns just a string "(3 4 5 )" 

// //  obsolete !!!!!!!
// void InputParser::parse_list_phys_ID(ifstream& in_stream, string& list_regions   )

// {



//   string str ;


//   vector<string>  v_string;



//   v_string.clear();


 
//   rule<>special_char =  (ch_p('_') | ch_p('-') |  ch_p('.') |  ch_p('/')   |  ch_p('+') | ch_p(',')    );

 
//   rule<>label  = *(special_char)>>  (+alnum_p)>>   * ( (special_char ) >> *(+alnum_p) ) ;


//   // ***********************************************************************
//   //      1-12-06  !!!

//   rule<>list_string =  ch_p('(') >> *(space_p) >> (label)>>  *( *(space_p) >> (label) ) >> ch_p(')');
//   rule<>tag_value =  list_string | label;
//   //  rule<>assignement_string  =  (label) >>  *(space_p) >> ch_p('=')>> *(space_p) >> 
//   //    ((tag_value)[push_back_a(v_string)])   ;
//   rule<>assignement_string  =  *(space_p) >> ch_p('=')>> *(space_p) >> 
//     ((tag_value)[push_back_a(v_string)])   ;

//   rule<> r_command  = *(space_p) >> (assignement_string )    >> 
//     *(space_p) >> !(comment_p("#") >> *(anychar_p)) >>*(space_p)  ;

//   //      14-2-07  !!!



 

//   //  rule<>list_of_strings_space_sep  = ch_p('(')>> *(space_p) >> (label)[push_back_a(v_string)]>>
//   //     *( *(space_p) >> (label)[push_back_a(v_string)] )>> *(space_p)>> ch_p(')')  ;

  
//   //   rule<>assignement = (label) >> 
//   //     *(space_p) >> ch_p('=')>> *(space_p) >> list_of_strings_space_sep;

//   //   rule<> r_command  = *(space_p) >> (assignement )    >> 
//   //     *(space_p) >> !(comment_p("#") >> *(anychar_p)) >>*(space_p)  ;



//   // ************************************************************************

 
 


//   while ( getline(in_stream, str) )
//   {

   
//     if  (!(parse(str.c_str(), comment_p("#")   , space_p).full) ) 

//     { // if !  comment_p("#")  

	 


//       if(  parse(str.c_str(),  r_command   )   .full ) //  not skipping spaces
//         // reads  list of  values (vector)

//       {

//         if ( !(v_string.empty()) )
//         {
//           //  vector_prop_labels.push_back(vect_label[0]);
		  
//           //           vector_prop_labels_map.insert(make_pair(vect_label[0], v) );
//           list_regions = v_string[0];
//           break;   //  return to  input file

//         }

	     

	      
       

//         v_string.clear();

        

//       }

     

//       //     skip   all  the  other  kids  of  lines  !!!!
      
//       else if (parse(str.c_str(), if_p("{")[(+alpha_p) 
//                                             ] , space_p ).full)
//       {  
//         //   cout << "SKIP" << endl ;
//         //   if (name == "a")
//         //   break;
//       }


//       else  

//       {


//         cerr <<  "  SYNTAX ERROR in input  file   " <<  endl;
//         cerr << " Correct syntax is : 'label' = 'value' 'label' = 'value' .......# 'comment' " 
//              << endl;
//         cerr << " A comment line  must be preceded by '#'. "<< endl;
//         cerr <<"'Value' is a string of alphanumerics char or a list of strings between par., e.g. (p1 p2)"
//              << endl;
//         exit(1); 
//       }


//     }


//     v_string.clear();

   

  
//   }  //  end  while


// }



bool InputParser::check_label(set<string>& section_keywords, const string& label )
{

  set<string>::iterator   s_it;
  bool  found;
  found = false;

  
  //for ( s_it = section_keywords.begin( ); s_it != section_keywords.end( ); s_it++ )
  //     cout << " " << *s_it;
  //   cout << "." << endl;

  s_it = section_keywords.find(label);

  if (s_it != section_keywords.end())
  {

    //   cout <<  " Found  " <<  endl;

    found = true;
 

  } 
  else 
  {
    //   cout <<  " NOT  Found  " <<  endl;
    found = false;
  } 


  return found;

} //  end  method




//void InputParser::cut_off_comment(string& label )
void InputParser::cut_off_comment(string& label,  ifstream& in_stream)
{

  string::size_type loc = label.find_first_of( "#", 0 );
  if( loc != string::npos )
  {
    //    cout << "Found # at " << loc << endl;
    label.erase(loc);

    in_stream.ignore(256,'\n');  //  if  the   read keyword is # or  begins with #: ignore  all  the  line

  }
  //  else
  //    cout << "Didn't find # " << endl;        
  // lab1.erase(loc);

} //  end  method




//  to handle CR DOS files line termination  
void InputParser::cut_off_CR(string& label,  ifstream& in_stream)
{

  string::size_type loc = label.find_first_of( "CR", 0 );
  if( loc != string::npos )
  {
    //    cout << "Found # at " << loc << endl;
    label.erase(loc);

    in_stream.ignore(256,'\n');  //  if  the   read keyword is # or  begins with #: ignore  all  the  line

  }
  //  else
  //    cout << "Didn't find # " << endl;        
  // lab1.erase(loc);

} //  end  method



// ****************************************



bool InputParser::skip_to_bracket(ifstream& in_stream)

{

  bool check_error = false;  
  char ch ;

  //*************************


      in_stream.get(ch); // get next char
  do{
      

    //********************************
      if (ch == '#' )  //  skip  comments 
    { 
      in_stream.ignore(256,'\n');  //  if  the   read keyword is # or  begins with #: ignore  all  the  line 
      //      and   read  the  next keyword !!!
    }
    in_stream.get(ch);


  } while ( (ch == '\n') || (ch == '\r')
            || (ch == ' ') || (ch == '\t') || (ch == '#') );


  // then  check  if  first char (not  blank) is  "{"
  if (ch == '{')
  {
    //  cout  <<  "  ch { :   "<< ch <<  endl << endl;
    //  OK
   
    return  check_error ;

  }

  else  
  {
    in_stream.unget();

    check_error =true;
    return  check_error ;

  }

} //  end  method



//  *************************************************
//  read  $Scale section 


void InputParser::read_scale(void)
{// read_scale



  std::string  label, keyword,region_name, atomistic_region_name, section_name   ;
  std::ifstream in_stream (filename.c_str()) ;
  ID atomistic_region_counter ;
  bool  check_error;
  bool found = false;
  atomistic_region_counter = 0;


  assert(in_stream.good());

  section_name =  "Scale";

  section_name = "$"+section_name;

  //  find_keyword( in_stream,section_name  );
  found = find_optional_keyword( in_stream,section_name  );

  // ************************
  if  (found == false) return; // $Scale section not  found;  atomistic_regions_map will be  void
  // ************************

  check_error = skip_to_bracket(in_stream);  // go on reading until the  char  '{' is  read  
  if (check_error == true)
    throw InitFailedException("In input file: \'Scale\'  block  missing."); 

  //  read keyword Atomistic 
  in_stream >> keyword;



  //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
  while (skip_comments(in_stream,keyword ) == true )
  {
    in_stream >> keyword; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
  } 

  cut_off_comment(keyword, in_stream); //  case  Atomistic#commmm


  while  (keyword !=  end_symb)
  { // while

    //cerr << " keyword = "<< keyword << endl ;

    if  (keyword == "Atomistic")
    {
      //      *********** read  atomistic ***********************

      //     read   region_name
      in_stream >>atomistic_region_name;
      atomistic_region_counter++;
      //  if  the   read keyword is # or  begins with #: ignore  all  the  line !!
      while (skip_comments(in_stream,atomistic_region_name) == true )
      {
        in_stream >> atomistic_region_name; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
      } 

      cut_off_comment(atomistic_region_name, in_stream );  //  in  case  layer#commmm

      //    cerr << " ************** ATOMISTIC ************ " <<  atomistic_region_name << endl;

      temp_options.clear();
      //    parse_options(in_stream,temp_options  ); //    read  the  block  between  { and  }
      parse_options(in_stream,temp_options, keyword,  section_name);

      //....................

      // create new RegionStructure
      RegionStructure current_region_structure;

      string  atomistic_region_numb  , def1;

      // put  atomistic region_name  in RegionStructure
      current_region_structure.set_region_name(atomistic_region_name );

      current_region_structure.set_model_options(temp_options);
      //  temp_options contains  atomistic ID and all  other data (list of phis region ID) 

      //     atomistic_region_numb = temp_options.get_option( "atomistic_region_numb" ,def1);

      atomistic_regions_map.insert(make_pair (atomistic_region_counter, current_region_structure ));




    }

    else if (keyword == "Lumped")
    {
   

      throw InitFailedException("In input file (Scale section): keyword "
                                "\'Lumped\' not implemented! ");

      // TO  BE  IMPLEMENTED 

    }

    else 
    {
      throw InitFailedException("In input file (Scale section): keyword "
                                "\'Atomistic\' is missing! ");
    
    }


    //    next   block  
    in_stream >> keyword;

    //   skip_comments(in_stream,keyword );

    while (skip_comments(in_stream,keyword ) == true )
    {
      in_stream >> keyword; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 
    } 

    cut_off_comment(keyword, in_stream); //  case  Atomistic#commmm




  }//end   while 



} //  end   Scale   section 



// ********************************************************************




// parse a number n of  subblocks of kind 
//  "block_name"
//  { .........
//  }
//  and  put the contents in a map<block_name,ModelOption>
void
InputParser::parse_n_subblocks(ifstream& in_stream, ModelOptions& block_options)
{

  string  label, block_name ;

  // read  the  next  label

  in_stream >>  label ;

  while (skip_comments(in_stream, label ) == true )
  {
    in_stream >> label  ; // if  the  whole  line has
    //                       been  skipped: read  the next keyword !!! 
  } 

  cut_off_comment( label, in_stream); //  case  Recomb#commmm
 
  //  READ the block-keyword

  while (  ( label !=  end_symb    ) && (!in_stream.eof()) ) 
    
  {  //  while loop  blocks

    // block n
               
    block_name = label;

    cut_off_comment(block_name, in_stream); //

    temp_options.clear();

    //    ****** read  the  subblock  between  { and  } *************
    //   parse_options(in_stream,temp_options  );
    parse_options(in_stream,temp_options, label,label );

    // ****  crea  ModelOptions figlio
    block_options.add_submodel(block_name, temp_options);


    //  blocks_map.insert(make_pair (block_name,temp_options));


    // read  the  next  label
    in_stream >>  label ;
    while (skip_comments(in_stream,  label ) == true )
    {
      in_stream >> label  ; // if  the  whole  line has
      //  been  skipped: read  the next keyword !!! 
    } 
    cut_off_comment(label, in_stream); //  case  Recomb#commmm

  }


}


//   ---------------------------------------------------



void InputParser::read_subblocks(string section_name,
                                 map<string, ModelOptions>& options_map)

{

  string  keyword ;
  ifstream in_stream (filename.c_str()) ;

  ModelOptions   block_options;  //  local  object ModelOptions

  block_options.clear();

  assert(in_stream.good());

  section_name = "$"+section_name;
 
 
  //find section name
  find_keyword( in_stream,section_name  );

  // -----------------------------


  // // ModelOptions   temp_options

  //   if  (section_name == "$Simulation")
  //   {

  //     temp_options.clear();
  //     parse_options(in_stream,temp_options, keyword, section_name);  
  //     //  add  to  map <string,ModelOptions>  keyword, temp_options
  //     options_map.insert(make_pair (keyword,temp_options  ));

  //     return;

  //   }


  // -----------------------------

  // go to the  first '{'
  skip_to_bracket(in_stream);

  in_stream >> keyword;


  while (skip_comments(in_stream,keyword ) == true )
  {
    in_stream >> keyword; // if  the  whole  line has  ben  skipped: read  the  next keyword !!! 


  } 

  cut_off_comment(keyword, in_stream); //  case  Region#commmm

 

  while (  (keyword  !=  end_symb    ) && (!in_stream.eof()) ) 
    
  {  //  while loop  blocks

    // block n

    if  ( (keyword == "Sweep")  || (keyword == "Selfconsistent" )  )
      // reserved label ->  parse_n_subblocks
    {

      skip_to_bracket(in_stream);

      block_options.clear();
      //read all the  subblocks of  type "block_type" in the section
      //  (included header)
      parse_n_subblocks(in_stream,block_options);
      // 

   
      //  add  to  map <string,ModelOptions>  keyword, temp_options
      options_map.insert(make_pair (keyword,block_options  ));

    }

    else 

      // simulation/model  name  -> parse_options
    {

     

      block_options.clear();
      parse_options(in_stream,block_options, keyword, section_name);  
      //  add  to  map <string,ModelOptions>  keyword, temp_options
      options_map.insert(make_pair (keyword,block_options));
 
    }

    // read  the  next keyword  
    in_stream >> keyword ;
    while (skip_comments(in_stream,keyword  ) == true )
    {
      in_stream >> keyword   ; // if  the  whole  line has
      //  been  skipped: read  the next keyword !!! 
    } 
    cut_off_comment(keyword, in_stream); //  case  Recomb#commmm

    //  next  block
  } 

}



void InputParser::get_solver_options_map(
                                         map<string, ModelOptions>& options_map)
{
  read_subblocks("Solver", options_map);
}


void InputParser::get_physics_options_map(
                                          map<string, ModelOptions>& options_map )
{
  read_subblocks("Physics", options_map);
}


void InputParser::get_simulation_options( ModelOptions& options )
{


  string  keyword, section_name ;
  ifstream in_stream (filename.c_str()) ;

 
  assert(in_stream.good());

  section_name = "$Simulation";
 

  //find section name
  find_keyword( in_stream,section_name  );

  keyword = section_name;
  parse_options(in_stream, options, keyword, section_name); 

  

}



void InputParser::skip_block(ifstream& in_stream)

{

  bool check_error = false;  
  char ch ;

  //*************************


      in_stream.get(ch); // get next char
  do{
         
    in_stream.get(ch);
  
  } while (ch != '}');


} //  end  method



const ModelOptions&
InputParser:: get_device_options(void)

{

  return _device_options;


} //  end  method

void InputParser:: set_device_options(ModelOptions& device_opt)

{

  _device_options = device_opt;


} //  end  method
