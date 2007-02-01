#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <map>
#include "ModelStructure.h"


ModelStructure::ModelStructure(string& model_name)
{
  set_model_name(model_name);	
}

ModelStructure::~ModelStructure()
{

}


map  <unsigned int,  map <const string,string>  >&  ModelStructure::get_model_BC_map()

{
 
  return model_BC_map;

}

map  <const string,string>&  ModelStructure::get_phys_model_map()

{
 
  return phys_model_map;

}
 

vector<string> ModelStructure::get_physical_regions()
{
 
  return physical_regions;
}
 

void ModelStructure::set_model_name(string& model)
{
 
  model_name=model ;
}

string  ModelStructure::get_model_name()
{
 
  return model_name;


}


void ModelStructure::set_model_BC_map( map  <unsigned int,  map <const string,string> >& id_BC_regions_map   )
{
 
  model_BC_map = id_BC_regions_map ;

  //  cout <<  "(model_BC_map[1])[value] =  "  <<   (model_BC_map[1])["value"] <<  endl;
}


void ModelStructure::set_phys_model_map( map <const string,string>& physical_model_map )
{
 
 
  phys_model_map = physical_model_map;
 
}



void ModelStructure::set_physical_regions( vector<string>& list_phys_regions)
{
 

  physical_regions = list_phys_regions;

}
