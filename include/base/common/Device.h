
#ifndef __DEVICE_H__
#define __DEVICE_H__

#include <iostream>  
#include <sstream>
#include <vector>
#include <string>
//#include "RegionDefinition.h"
#include "BcRegionDefinition.h"
#include "MaterialRegion.h"
//#include "DeviceRegion.h"
#include "Material.h"
#include "Alloy.h"

using namespace std;

//! Higher-level definition of the  structure to  be  simulated.
/*!
 *  Structure of  the  device  regions (input data) is  read  with \c set_material_regions.
 *   With \c set_materials the  necessary materials (maybe  alloys) are  created.
 * \c set_map_ID_material_region instatiates an object \c MaterialRegion for  each 
 * RegionDefinition present in 
 * the  device and makes a map to  associate each \c MaterialRegion  to  the  region ID. 
 * Also, each \c MaterialRegion has  a  pointer to  the  corresponding Material
 * 
 * 
 * 
 */
 
class Device{



 public:

  Device();

 
  ~Device();

  //! Reads  a list of  material regions 
  /*!
   * Used to  pass  description of physical regions  (from input)  to  \c Device  object (through \c RegionDefinition objects).
   * 
   */
  void set_material_regions( const  vector<RegionDefinition>& dev_reg   );
 

  //! Gets  \c  MaterialRegion corresponding  to  region  ID.
  /*!
   * Public  method  to get \c  MaterialRegion corresponding  to  region  ID.
   * 
   */
  MaterialRegion* get_material_region(unsigned int  ID);
  //  DopingRegion* get_doping_region(unsigned int  ID); //  to be  implemented
  //  ModelRegion* get_model_region(unsigned int  ID);  //  to be  implemented




  

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

  //! Defines list of materials 
  /*!
   * Defines all \c Material objects to be  associated with  present  \c Device.
   * If  material is  an  alloy, instantiates an \c Alloy  object with the  relative pointers to 
   * alloy material components .
   */
  void set_materials();
	
	
	
  //! Defines a map of  \c MaterialRegion 
  /*!
   *   Makes all \c MaterialRegion objects  and associates them with a map with  region IDs .
   * 
   */
  void set_map_ID_material_region();
  void set_map_ID_doping_region();
  void set_map_ID_model_region();
 


  /* struct RegionDefinition { */     //     ->   class   RegionDefinition

  /*     string reg_name ; */
  /*     unsigned int  reg_numb; */
  /*     string  mat_name; */
  /*     double  dop_conc  ; */
  /*     string dop_type; */

  /*   }; */

  //   //  region_definition region_def;

  // vector<RegionDefinition> device_regions;
  
  
  /*!
   *   local  structures with  RegionDefinition (input data)
   */
  vector<RegionDefinition> material_regions;

  vector<BcRegionDefinition> BC_device_regions;

 

  //***************************************

      //DeviceRegion*  dev_region_point;

      /* vector<MaterialRegion> mat_region_vec; */
    /* vector<DopingRegion> dop_region_vec; */
    /* vector<ModelRegion> model_region_vec; */


    /*!
     *   Pointer to \c MaterialRegion  object.
     */
    MaterialRegion* mat_reg_point;  
 
    /*!
     *   Map  associating  \c MaterialRegion  objects  with  region IDs (from input).
     */
    map < unsigned int,MaterialRegion* > ID_mat_reg_map;

    // TO BE IMLEMENTED
    /* DopingRegion* dop_reg_point; */
    /*  map < unsigned int,DopingRegion* > ID_dop_reg_map; */

    /* ModelRegion* model_reg_point; */
    /*  map < unsigned int,ModelRegion* > ID_model_reg_map; */


    /*!
     *   Map  associating  \c Material  objects  with  material names (from input).
     */
    map < string ,Material* > name_mat_map;
    //	Material* matpoint, matpoint_alloy ;

    /*!
     *   Pointers to \c Material  and  \c Alloy objects.
     */
    Material* matpoint;
    Alloy* matpoint_alloy ;


    /*!
     *   Info on the material components of the  alloy (to be read from database) .
     */
    vector<string>  alloy_components;

};


//inline Device::~Device()
//{
//
////  for   all  MaterialRegions in  map :
////  delete  p
//
//
//}


//inline Device::Device()
//{
//	cout <<  "ciao" << endl;
//}


#endif //  #define __DEVICE_H__  
