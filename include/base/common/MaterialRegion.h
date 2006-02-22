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

//! Contains all input-derived properties   of  a material region
/*!
 * A set of  basic  properties and a  pointer to the \c Material 
  * present in the  region is  associated   to  any 
  * physical region defined in the  device  structure.
 */
class MaterialRegion:public DeviceRegion    
{

  public:

    //! Constructor
  /*!
   * In the  constructor  a \c RegionDefinition and 
    * a pointer to the  region material  are  passed to  the  class.  
 */
    MaterialRegion(const RegionDefinition&  mat_reg, Material* mat_pointer);
  // MaterialRegion(const RegionDefinition&  mat_reg);

  //! Void  Destructor
   /*!
   */
    ~MaterialRegion();

    //! Returns a pointer to the material associated with this \c MaterialRegion.
    /*!
   */
    Material*  get_material() const;
  
    
    //! Returns id  number of  this \c MaterialRegion.
    /*!
     */
    unsigned int get_region_number() const ;
    //! Returns  name of  this \c MaterialRegion.
    /*!
     */
    string get_region_name() const ;
    //! Returns name of the material  of  this \c MaterialRegion.
    /*!
     */
    string   get_material_name() const;
    //! Returns name of the crystal structure  of  this \c MaterialRegion.
    /*!
     */
    string   get_crystal_name()  const;
    //! Returns the doping concentration in    this \c MaterialRegion.
    /*!
     */
    double  get_doping_concentration()   const;
    //! Returns the type ( p or n ) of  doping in    this \c MaterialRegion.
    /*!
     */
    string  get_doping_type()  const;
  

  private:

    //! Description of the material region through a \c RegionDefinition object.
    /*!
     */
    RegionDefinition  material_region_definition;

    //! Pointer to the  \c Material  associated to this \c MaterialRegion. 
    /*!
     */
    Material* matpoint ;


};

//
// inline  members functions
//


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


inline
    Material* MaterialRegion::get_material() const
{
  
  return   matpoint;
  
}



inline 
    unsigned int MaterialRegion::get_region_number() const
{


  return material_region_definition.get_region_number();

}  

inline 
    string  MaterialRegion::get_region_name() const
{


  return material_region_definition.get_region_name();

}  


inline 
    string MaterialRegion::get_material_name() const
{


  return material_region_definition.get_material_name();

}  


inline
    string  MaterialRegion::get_crystal_name() const
{
  
  return material_region_definition.get_crystal_name() ;
  
}


inline
    double  MaterialRegion::get_doping_concentration()   const
{
  
  return material_region_definition.get_doping_concentration() ;
  
}

inline
    string  MaterialRegion::get_doping_type()  const
{
  
  return material_region_definition.get_doping_type();
  
}
  
  


    
#endif //  

