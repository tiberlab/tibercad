#ifndef __REGIONDEFINITION_h__
#define __REGIONDEFINITION_h__

#include <string>
#include <iostream>  
#include <sstream>
#include <vector>
#include <string>
using namespace std;


//! Contains all  data of a physical region read from input  device  description. 
/*!
 * A set of  basic  properties of  a  single  physical region is  read 
 * from input  device  description in a  \c  RegionDefinition object.
 */
class RegionDefinition{

 public:

   //! Void constructor. 
  /*!
 */
  RegionDefinition();

  //! Void destructor. 
  /*!
   */
  ~RegionDefinition();



//! Set the name of the region. 
/*!
 *
 */
  void set_region_name(string& field_value);
  //! Set the number id of the region. 
/*!
   *
 */
  void set_region_number(unsigned int field_value);
  //! Set the name of the material of the region. 
/*!
   *
 */
  void set_material_name(string& field_value);
  //! Set the name of the crystal structure of the material in the region. 
/*!
   *
 */
  void set_crystal_name(string& field_value);
  //! Set the doping concentration  of the region. 
/*!
   *
 */
  void set_doping_concentration(double field_value);
  //! Set the doping type of the region. 
/*!
   *
 */
  void set_doping_type(string& field_value);

//! Returns  the number id of the region. 
/*!
   *
 */
  unsigned int get_region_number() const;
  //! Returns  the name  of the region. 
/*!
   *
 */
  string  get_region_name() const;
  
  //! Returns  the name of the material  of the region. 
/*!
   *
 */
  string  get_material_name()const ;
  //! Returns  the crystal structure of the material of the region. 
/*!
   *
 */
  string  get_crystal_name() const;
  //! Returns  the  doping concentration in the region. 
/*!
   *
 */
  double  get_doping_concentration() const ;
  //! Returns  the doping type of the region. 
/*!
   *
 */
  string  get_doping_type() const ;


 private:

//! Private member data. 
/*!
  *
 */

  string region_name ;
  unsigned int  region_number;
  string  material_name;
  string  crystal_name;
  double  doping_concentration  ;
  string doping_type;

};

// ------------------------------------------------------------
// RegionDefinition  inline members
// ---------------------------------------------------------------

inline unsigned int RegionDefinition::get_region_number() const

{

  return region_number;

}

inline string  RegionDefinition::get_region_name() const

{

  return region_name;

}

inline string  RegionDefinition::get_material_name() const

{

  return material_name;

}

inline string  RegionDefinition::get_crystal_name() const

{

  return crystal_name;

}    
    

inline double  RegionDefinition::get_doping_concentration() const

{

  return doping_concentration;

}


inline string  RegionDefinition::get_doping_type() const

{

  return doping_type;

}

// *******************************************************



inline void RegionDefinition::set_region_name(string& field_value)

{

  region_name = field_value;

}

inline void RegionDefinition::set_region_number(unsigned int field_value)

{

  region_number = field_value;

}

inline void RegionDefinition::set_material_name(string& field_value)

{

  material_name = field_value;

}

inline void RegionDefinition::set_crystal_name(string& field_value)

{

  crystal_name = field_value;

}

inline void RegionDefinition::set_doping_concentration(double field_value)

{

  doping_concentration = field_value;

}
inline void RegionDefinition::set_doping_type(string& field_value)

{

  doping_type = field_value;

}








#endif // #define define __REGIONDEFINITION_h__ 
