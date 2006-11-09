#ifndef __SOLVERPARAMETERS_H__
#define __SOLVERPARAMETERS_H__

#include <string>
#include <iostream>  
#include <sstream>
#include <vector>
#include <string>
#include <map>
using namespace std;

class SolverParameters{

 public:

  SolverParameters();
  ~SolverParameters();

  void set_parameters( map <string,double>& num_map, map <string,string>&  string_map, 
				    map <string, vector<double> >&  vector_map   );


  // *************************************************************

  template <typename T> 
    T  get_parameters_tmp(string&  label) 

{

  // template .... get_param_value(string)
  //  used in  solver  module to  get relevant  parameters 
 
  T  value,Default ;
 
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






  // ****************************************************************

  
  int   get_parameters(string&  label, int  Default)  ;

  double  get_parameters(string&  label, double  Default) ;

  string  get_parameters(string&  label, string  Default ) ;

  vector<double>  get_parameters(string&  label,vector<double>& def_vector);
  vector<unsigned int>  get_parameters(string&  label,vector<unsigned int>& def_vector);




 /*  unsigned int RegionDefinition::get_region_number() const; */
/*   string  RegionDefinition::get_region_name() const; */
/*   string  RegionDefinition::get_material_name()const ; */
/*   double  RegionDefinition::get_doping_concentration() const ; */
/*   string  RegionDefinition::get_doping_type() const ; */




 private:

  map <string,double>  numerical_param_map;
  map <string,string>  string_param_map;

  map <string, vector<double> >  vector_param_map;




};



// *******************************************************





//**************************************************************





/* inline unsigned int BcRegionDefinition::get_BC_region_number() const */

/* { */

/*   return BC_region_number; */

/* } */


//**************************************************




#endif // #define define __SOLVERPARAMETERS_H__ 
