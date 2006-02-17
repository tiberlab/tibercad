#ifndef __BCREGIONDEFINITION_H__
#define __BCREGIONDEFINITION_H__

#include <string>
#include <iostream>  
#include <sstream>
#include <vector>
#include <string>
using namespace std;

class BcRegionDefinition{

 public:

  BcRegionDefinition();

  ~BcRegionDefinition();




 /*  void BcRegionDefinition::set_region_name(string& field_value); */
/*   void BcRegionDefinition::set_region_number(unsigned int field_value); */
/*   void BcRegionDefinition::set_material_name(string& field_value); */
/*   void BcRegionDefinition::set_doping_concentration(double field_value); */
/*   void BcRegionDefinition::set_doping_type(string& field_value); */

  void set_BC_region_name(string& field_value);
  void set_BC_region_number(unsigned int field_value);
  void set_BC_type(string& field_value);
  void set_BC_value(double field_value);




 /*  unsigned int RegionDefinition::get_region_number() const; */
/*   string  RegionDefinition::get_region_name() const; */
/*   string  RegionDefinition::get_material_name()const ; */
/*   double  RegionDefinition::get_doping_concentration() const ; */
/*   string  RegionDefinition::get_doping_type() const ; */

  unsigned int get_BC_region_number() const;
  string  get_BC_region_name() const;
  string  get_BC_type() const;
  double  get_BC_value() const;



 private:


/*   /* struct RegionDefinition { */     //     ->   class   RegionDefinition */

/*   string region_name ; */
/*   unsigned int  region_number; */
/*   string  material_name; */
/*   double  doping_concentration  ; */
/*   string doping_type; */



  //  BC regions

 string BC_region_name ;
 unsigned int BC_region_number;
 string  BC_type;
 double BC_value;





};



// *******************************************************





//**************************************************************

// BC regions inline  members



inline unsigned int BcRegionDefinition::get_BC_region_number() const

{

  return BC_region_number;

}

inline string  BcRegionDefinition::get_BC_region_name() const

{

  return BC_region_name;

}

inline string  BcRegionDefinition::get_BC_type() const

{

  return BC_type;

}

inline double  BcRegionDefinition::get_BC_value() const

{

  return  BC_value;

}






// ------------------------------------------------



inline void BcRegionDefinition::set_BC_region_name(string& field_value)

{

  BC_region_name = field_value;

}

inline void BcRegionDefinition::set_BC_region_number(unsigned int field_value)

{

  BC_region_number = field_value;

}

inline void BcRegionDefinition::set_BC_type(string& field_value)

{

  BC_type = field_value;

}

inline void BcRegionDefinition::set_BC_value(double field_value)

{

  BC_value   = field_value;

}


//**************************************************




#endif // #define define __BCREGIONDEFINITION_H__ 
