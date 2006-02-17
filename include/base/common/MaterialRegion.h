#ifndef __MATERIALREGION_H__
#define __MATERIALREGION_H__

#include <iostream>  
#include <sstream>
#include <vector>
#include <string>
#include "DeviceRegion.h"
#include "Material.h"
#include "RegionDefinition.h"

//using namespace std;


class MaterialRegion:public DeviceRegion    
{

 public:

  MaterialRegion(const RegionDefinition&  mat_reg, Material* mat_pointer);
  // MaterialRegion(const RegionDefinition&  mat_reg);


  ~MaterialRegion();

  unsigned int get_region_number();


 private:

  RegionDefinition  material_region_definition;

  Material* matpoint ;


};


//inline  MaterialRegion::MaterialRegion(const RegionDefinition&  mat_reg)
inline  MaterialRegion::MaterialRegion(const RegionDefinition&  mat_reg, Material* mat_pointer)
{
  material_region_definition =  mat_reg;
  //  string mat_name = material_region_definition.get_material_name();

  // matpoint = new Material( mat_name); 
  
  matpoint = mat_pointer;
  

}



inline   MaterialRegion::~MaterialRegion()
{

  delete matpoint;

}


#endif //  

