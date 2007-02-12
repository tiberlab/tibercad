#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <map>
#include "ModelStructure.h"

using namespace std;

ModelStructure::ModelStructure(const string& model_name)
{
  set_model_name(model_name);	
}

ModelStructure::~ModelStructure()
{

}


// map  <unsigned int,  map <const string,string>  >&  ModelStructure::get_model_BC_map()

// {
 
//   return model_BC_map;

// }

map <ID, RegionStructure>& ModelStructure::get_model_BC_map()
{
 
  return model_BC_map;

}


// map  <const string,string>&  ModelStructure::get_phys_model_map()

// {
 
//   return phys_model_map;

// }
 
// const ModelOptions& ModelStructure::get_phys_model_options( )
// {
 
//   return physical_model_options  ;
 
// }

multimap<const string,ModelOptions>& ModelStructure::get_physical_model_map()
{
 
  return  physical_model_map;

}



vector<string> ModelStructure::get_physical_regions()
{
 
  return physical_regions;
}
 

void ModelStructure::set_model_name(const string& model)
{
 
  model_name=model ;
}

string  ModelStructure::get_model_name()
{
 
  return model_name;


}


//void ModelStructure::set_model_BC_map( map  <unsigned int,  map <const string,string> >& id_BC_regions_map   )
void ModelStructure::set_model_BC_map( map <ID, RegionStructure>& mod_BC_map)
{
 
  model_BC_map =mod_BC_map  ;

  //  cout <<  "(model_BC_map[1])[value] =  "  <<   (model_BC_map[1])["value"] <<  endl;
}


// void ModelStructure::set_phys_model_map( map <const string,string>& physical_model_map )
// {
 
 
//   phys_model_map = physical_model_map;
 
// }


// void ModelStructure::set_phys_model_options(  ModelOptions& phys_model_options  )
// {
 
//   physical_model_options = phys_model_options ;
 
// }

void ModelStructure::set_physical_model_map( multimap <const string,ModelOptions>& phys_model_map)
{
 
  physical_model_map   = phys_model_map ;
 
}





void ModelStructure::set_physical_regions( vector<string>& list_phys_regions)
{
 

  physical_regions = list_phys_regions;

}
