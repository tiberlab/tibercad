
#ifndef __DEVICE_H__
#define __DEVICE_H__

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
//#include "RegionDefinition.h"
#include "BcRegionDefinition.h"
#include "MaterialRegion.h"
//#include "DeviceRegion.h"
#include "Material.h"
#include "Alloy.h"
#include "AlloyModel.h"

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

class Device
{

 public:
  
  //! Constructor. 
  /*!
   * 
   * 
   */
  Device ();

  //! Void  Destructor. 
  /*!
   * 
   * 
   */
  ~Device ();

  
  
  //! Initialize  \c  Device  object. 
  /*!
   * Public  method  to initialize  a  \c  Device  with a vector of  physical regions  
   * \c RegionDefinition  (read from input file).
   * 
   */
  void init_device(const  vector<RegionDefinition>& dev_reg);
  
      
  //! Gets  \c  MaterialRegion corresponding  to  region  ID.
  /*!
   * Public  method  to get \c  MaterialRegion corresponding  to  region  ID.
   * 
   */
  MaterialRegion *get_material_region (unsigned int ID);
  //  DopingRegion* get_doping_region(unsigned int  ID); //  to be  implemented
  //  ModelRegion* get_model_region(unsigned int  ID);  //  to be  implemented



  /*  void Device::set_device_structure(vector<string>& region_name_v, vector<unsigned int>&  region_number_v,  */
  /*                                vector<string>&  material_name_v, vector<double>& doping_concentration_v,  */
  /*                                vector<string>&  doping_type_v   ); */

  // vector<RegionDefinition> device_regions;

  /* string  Device::get_device_data(unsigned int reg_query, string& field); */

  /*  double  Device::get_device_data2(unsigned int reg_query, string& field); */

  const RegionDefinition & get_device_data (unsigned int region_query);

  void set_device_boundary_cond (vector < string > &BC_region_name_v,
                                 vector < unsigned int >&BC_region_number_v,
                                 vector < string > &BC_type_v,
                                 vector < double >&BC_value_v);


  const BcRegionDefinition & get_device_boundary_cond (unsigned int
                                                       region_query);
        
        
  //! Sets  map reg_alloy_model_map.
  /*!
   * Public  method  to write in \c Device a map which associates
   *  region ID with \c AlloyModel objects.
   * 
   */
  void set_map_alloy_model(const map <unsigned int, AlloyModel*>&  map_alloy_model);
  


  //! Returns set of all materials of the  device.
  /*!
   * Public  method  to read the  \set of all the  
   * \c  Material and/or  \c Alloy objects  created in \c  Device.
   * 
   */
  set<Material *>&  get_set_all_materials();

 private:
  
  // ------------------
  // PRIVATE METHODS
  //-------------------
  
  //! Reads  a list of  material regions 
  /*!
   * Used to  pass  description of physical regions  (from input)  to  \c Device  object (through \c RegionDefinition objects).
   * 
   */
  void set_material_regions (const vector < RegionDefinition > &dev_reg);
  
  

  //! Defines list of materials 
  /*!
   * Defines all \c Material objects to be  associated with  present  \c Device.
   * If  material is  an  alloy, instantiates an \c Alloy  object with the  relative pointers to 
   * alloy material components .
   */
  void set_materials ();



  //! Defines a map of  \c MaterialRegion 
  /*!
   *   Makes all \c MaterialRegion objects  and associates them with a map with  region IDs .
   * 
   */
  void set_map_ID_material_region ();
    
  // TO BE implemented
  void set_map_ID_doping_region ();
  void set_map_ID_model_region ();
    
  void makes_set_of_materials(); 
    
  //-----------------------------
  //  PRIVATE  Data MEMBERS
  //-----------------------------
  //
  //  MAPS
  //
  /*!
   *   Map  associating  \c Material  objects  with  material names (from input).
   */
  map < string, Material* >name_mat_map;
    
  /*!
   *   Map  associating  \c Alloy  objects  with region ID(from input).
   */
  map < unsigned int, Alloy* >reg_id_alloy_map;
    
  /*!
   *   Map  associating  \c MaterialRegion  objects  with  region IDs (from input).
   */
  map < unsigned int, MaterialRegion * >ID_mat_reg_map;

  // TO BE IMPLEMENTED
  /* DopingRegion* dop_reg_point; */
  /*  map < unsigned int,DopingRegion* > ID_dop_reg_map; */

  /* ModelRegion* model_reg_point; */
  /*  map < unsigned int,ModelRegion* > ID_model_reg_map; */
    
    
  /*!
   *   Map  associating  region IDs  with \c AlloyModel  objects \c.
   */
  map <unsigned int, AlloyModel*>  reg_alloy_model_map;
    
    
    
  /*!
   *   Set containing all the \c Material objects  present in name_mat_map and reg_id_alloy_map
   * that is all the  simple Materials plus the  Alloy material plus the alloy components.
   */
  set<Material *>  set_all_materials;
    
  //        
  // local vectors
  //
  /*!
   *   Local  structures with  RegionDefinition (input data)
   */
  vector < RegionDefinition > material_regions;

  /*!
   *   Local  structures with  BcRegionDefinition (input data)
   */
  vector < BcRegionDefinition > BC_device_regions;
    
  /*!
   *   Names of  the material components of the  alloy (to be read from database) .
   */
  vector < string > alloy_components;



  //
  // Pointers
  //
  /*!
   *   Pointer to \c MaterialRegion  object.
   */
  MaterialRegion * mat_reg_point;

  /*!
   *   Pointers to \c Material  object.
   */
  Material *matpoint;
    
  /*!
   *   Pointers to \c Alloy object.
   */
  Alloy *matpoint_alloy;
    
    
  /*!
   *   Pointers to \c AlloyModel object.
   */
  AlloyModel* alloy_model_point;

  
    
    
   
    
  
  
  /* struct RegionDefinition { *///     ->   class   RegionDefinition

  /*     string reg_name ; */
  /*     unsigned int  reg_numb; */
  /*     string  mat_name; */
  /*     double  dop_conc  ; */
  /*     string dop_type; */

  /*   }; */

  //   //  region_definition region_def;

  // vector<RegionDefinition> device_regions;





  //***************************************

      //DeviceRegion*  dev_region_point;

      /* vector<MaterialRegion> mat_region_vec; */
    /* vector<DopingRegion> dop_region_vec; */
    /* vector<ModelRegion> model_region_vec; */

    

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
//      cout <<  "ciao" << endl;
//}


#endif //  #define __DEVICE_H__
