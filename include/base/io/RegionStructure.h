#ifndef _REGIONSTRUCTURE_H_
#define _REGIONSTRUCTURE_H_


#include <iostream>
#include <fstream>

#include <map>
#include <vector>
#include <string>
#include "TypeDefs.h"
#include "ModelOptions.h"



//!  Class  containing  input-defined REGION STRUCTURE.
/*!
 * Contains  name,  region ID, material, Options. 
 */

class RegionStructure{


 public:

  //!  Constructor. 
  /*!
   * Creates a new region (physical or  BC) structure .
   */
  RegionStructure();

  //!  Destructor 
  /*!
   *
   */
  ~RegionStructure();

  //!  Returns the  region name. 
  /*!
   * 
   */
  std::string  get_region_name();

  //!  Returns the  region ID. 
  /*!
   * 
   */
  std::string  get_region_ID();

  //!  Returns the name of the material associated to the region (if  physical) . 
  /*!
   * 
   */
  std::string  get_material_name();

  //!  Returns the options associated to the region. 
  /*!
   * Returns a reference to object ModelOptions.
   */
  const ModelOptions&  get_options();

 

  
  //!  Set the name of the  region to "region". 
  /*!
   * 
   */
  void set_region_name(const std::string& region);


  //!  Set the region ID to "region_numb". 
  /*!
   * 
   */
  void set_region_ID(const std::string& region_numb);


  //!  Set the material name  to "material". 
  /*!
   * 
   */
  void set_material_name(const std::string& material);


  //!  Set  the options for the present  region. 
  /*!
   * 
   */
  void set_model_options(ModelOptions& options);

  //----------------------------------------------------------------------------------


 private:

  //! Compulsory items of  region  description.
  /*!
   * 
   */
  std::string region_name, region_ID, material_name;

  //!  Contains a ModelOptions description of  region options
  /*!
   * 
   */
  ModelOptions    region_options;

};



inline const
ModelOptions&  RegionStructure::get_options()
{
   
  return region_options;

}



inline
std::string  RegionStructure::get_region_name()
{
   
  return region_name;

}


inline
std::string  RegionStructure::get_region_ID()
{
   
  return region_ID;

}


inline
std::string  RegionStructure::get_material_name()
{
   
  return material_name;

}


inline
void
RegionStructure::set_model_options(ModelOptions& options)
{
   
  region_options = options;

}


inline
void
RegionStructure::set_region_name(const std::string& region)
{
   
  region_name = region;

}


inline
void
RegionStructure::set_region_ID(const std::string& region_numb)
{
   
  region_ID = region_numb;  

}


inline
void
RegionStructure::set_material_name(const std::string& material)
{
   
  material_name =  material;  

}



#endif // endif define  _REGIONSTRUCTURE_H_
