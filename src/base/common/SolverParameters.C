#include <iostream>  
#include <sstream>
#include <vector>
#include <string>

#include "SolverParameters.h"
using namespace std;

SolverParameters::SolverParameters()
{

}


SolverParameters::~SolverParameters()
{

}



void SolverParameters::set_parameters( map <string,double>& num_map, map <string,string>&  string_map, 
				  map <string, vector<double> >&  vector_map   )
{

  // makes maps for  parameters read from  section in  input  file



  numerical_param_map =  num_map;
  string_param_map = string_map;
  vector_param_map = vector_map;

}



// template <typename T> 
// T  SolverParameters::get_parameters_tmp(string&  label) 

// {

//   // template .... get_param_value(string)
//   //  used in  solver  module to  get relevant  parameters 
 
//   T  value,Default ;
 
//   map <string,T>  :: iterator  p;

//   Default = 0;

  
//   p = numerical_param_map.find( label  );

//   if  (p != numerical_param_map.end() )

//     {  
//       value  =  (p -> second) ;
	

     
//         return value ;

//     }

//   else 
//     { 
//       cout  <<  "*** Default ***  ";

//          return Default;
//       //  cout  <<  "error"  ;
//     }
 

// }


// needs  instantiation for  each type  actually  implemented, otherwise put 
//  all  the  template 
//template <typename T> 
//T  SolverParameters::get_parameters_tmp(string&  label) 
//  in .h  file  (without  implement it  in  .C) 
//  so  that  template  can be  included with  include SolverParameters.h


// template
// double  SolverParameters::get_parameters_tmp<double>(string&  label);






 double   SolverParameters::get_parameters(string&  label, double  Default) 

{

  // template .... get_param_value(string)
  //  used in  solver  module to  get relevant  parameters 


 
   double value;
 
  map <string,double>  :: iterator  p;

  Default = 0;

  
  p = numerical_param_map.find( label  );

  if  (p != numerical_param_map.end() )

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
 

}





//int  InputParser::read_input( string label , int  Default)

int   SolverParameters::get_parameters(string&  label, int  Default ) 


{
  int value ;
  map <string,double>  :: iterator  p;


  //  for (int i =0; i< prop_labels.size();++i)
  //    {
  p = numerical_param_map.find( label  );

  if  (p != numerical_param_map.end() )

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




//string   InputParser::read_input( string label , string  Default)
string  SolverParameters::get_parameters(string&  label, string  Default ) 

{
  string  value;
  map <string,string>  :: iterator  p;


  //  for (int i =0; i< prop_labels.size();++i)
  //    {
  p = string_param_map.find( label  );

  if  (p != string_param_map.end() )

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



// overload  for   vector 

vector<double>  SolverParameters::get_parameters(string&  label,vector<double>& def_vector)

{

  map <string, vector<double> >   :: iterator  p;

  vector<double>  return_vector;

  p = vector_param_map.find( label  );

  if  (p != vector_param_map.end() )

    {  
      return_vector    =  (p -> second) ;
      return  return_vector;

    }

  else
    cout  <<  "Warning --- read_input_vector double: vector data empty or  missing  " << endl ;


}


vector<unsigned int>  SolverParameters::get_parameters(string&  label,vector<unsigned int>& def_vector)
{

  map <string, vector<double> >   :: iterator  p;

  vector<unsigned int>  return_vector;
  vector<double> temp;
  return_vector.clear();


  p = vector_param_map.find( label  );

  if  (p != vector_param_map.end() )

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

