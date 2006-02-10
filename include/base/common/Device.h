
#ifndef __DEVICE_H__
#define __DEVICE_H__

#include <iostream>  
#include <sstream>
#include <vector>
#include <string>
#include "RegionDefinition.h"
#include "BcRegionDefinition.h"
#include "MaterialRegion.h"
#include "DeviceRegion.h"


using namespace std;

class Device{



 public:

  Device();

 
  ~Device();




  MaterialRegion* get_material_region(unsigned int  ID);
  //  DopingRegion* get_doping_region(unsigned int  ID);
  //  ModelRegion* get_model_region(unsigned int  ID);




   void set_material_regions( const  vector<RegionDefinition>& dev_reg   );

 /*  void Device::set_device_structure(vector<string>& region_name_v, vector<unsigned int>&  region_number_v,  */
/* 				    vector<string>&  material_name_v, vector<double>& doping_concentration_v,  */
/* 				    vector<string>&  doping_type_v   ); */

  // vector<RegionDefinition> device_regions;

  /* string  Device::get_device_data(unsigned int reg_query, string& field); */

  /*  double  Device::get_device_data2(unsigned int reg_query, string& field); */

  const  RegionDefinition& get_device_data(unsigned int region_query);

  void set_device_boundary_cond(vector<string>&  BC_region_name_v,
					vector<unsigned int>& BC_region_number_v,
					vector<string>&  BC_type_v, 
					vector<double>& BC_value_v   );


  const  BcRegionDefinition&  get_device_boundary_cond(unsigned int region_query);



 private:


  /* struct RegionDefinition { */     //     ->   class   RegionDefinition

  /*     string reg_name ; */
  /*     unsigned int  reg_numb; */
  /*     string  mat_name; */
  /*     double  dop_conc  ; */
  /*     string dop_type; */

  /*   }; */

  //   //  region_definition region_def;

  // vector<RegionDefinition> device_regions;
  vector<RegionDefinition> material_regions;

  vector<BcRegionDefinition> BC_device_regions;

 

  //***************************************

      //DeviceRegion*  dev_region_point;

/* vector<MaterialRegion> mat_region_vec; */
/* vector<DopingRegion> dop_region_vec; */
/* vector<ModelRegion> model_region_vec; */



 MaterialRegion* mat_reg_point;  
 map < unsigned int,MaterialRegion* > ID_mat_reg_map;

/* DopingRegion* dop_reg_point; */
/*  map < unsigned int,DopingRegion* > ID_dop_reg_map; */

/* ModelRegion* model_reg_point; */
/*  map < unsigned int,ModelRegion* > ID_model_reg_map; */


 void set_map_ID_material_region();
 void set_map_ID_doping_region();
 void set_map_ID_model_region();




/*  map < unsigned int,DopingRegion > ID_doping_reg_map; */
/*  map < unsigned int,ModelRegion > ID_model_reg_map; */


// ******************************************


};


inline Device::~Device()
{

//  for   all  MaterialRegions in  map :
//  delete  p


}





#endif //  #define __DEVICE_H__  
