#ifndef __REGIONDEFINITION_h__
#define __REGIONDEFINITION_h__

#include <string>
#include <iostream>  
#include <sstream>
#include <vector>
#include <string>
using namespace std;

class RegionDefinition{

 public:

  RegionDefinition::RegionDefinition();

  RegionDefinition::~RegionDefinition();




  void RegionDefinition::set_region_name(string& field_value);
  void RegionDefinition::set_region_number(unsigned int field_value);
  void RegionDefinition::set_material_name(string& field_value);
  void RegionDefinition::set_doping_concentration(double field_value);
  void RegionDefinition::set_doping_type(string& field_value);

  void RegionDefinition::set_BC_region_name(string& field_value);
  void RegionDefinition::set_BC_region_number(unsigned int field_value);
  void RegionDefinition::set_BC_type(string& field_value);
  void RegionDefinition::set_BC_value(double field_value);




  unsigned int RegionDefinition::get_region_number() const;
  string  RegionDefinition::get_region_name() const;
  string  RegionDefinition::get_material_name()const ;
  double  RegionDefinition::get_doping_concentration() const ;
  string  RegionDefinition::get_doping_type() const ;

  unsigned int RegionDefinition::get_BC_region_number() const;
  string  RegionDefinition::get_BC_region_name() const;
  string  RegionDefinition::get_BC_type() const;
  double  RegionDefinition::get_BC_value() const;










 private:


  /* struct RegionDefinition { */     //     ->   class   RegionDefinition

  string region_name ;
  unsigned int  region_number;
  string  material_name;
  double  doping_concentration  ;
  string doping_type;



  //  BC regions

 string BC_region_name ;
 unsigned int BC_region_number;
 string  BC_type;
 double BC_value;





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

inline void RegionDefinition::set_doping_concentration(double field_value)

{

  doping_concentration = field_value;

}
inline void RegionDefinition::set_doping_type(string& field_value)

{

  doping_type = field_value;

}


//**************************************************************

// BC regions inline  members



inline unsigned int RegionDefinition::get_BC_region_number() const

{

  return BC_region_number;

}

inline string  RegionDefinition::get_BC_region_name() const

{

  return BC_region_name;

}

inline string  RegionDefinition::get_BC_type() const

{

  return BC_type;

}

inline double  RegionDefinition::get_BC_value() const

{

  return  BC_value;

}






// ------------------------------------------------



inline void RegionDefinition::set_BC_region_name(string& field_value)

{

  BC_region_name = field_value;

}

inline void RegionDefinition::set_BC_region_number(unsigned int field_value)

{

  BC_region_number = field_value;

}

inline void RegionDefinition::set_BC_type(string& field_value)

{

  BC_type = field_value;

}

inline void RegionDefinition::set_BC_value(double field_value)

{

  BC_value   = field_value;

}


//**************************************************




#endif // #define define __REGIONDEFINITION_h__ 
