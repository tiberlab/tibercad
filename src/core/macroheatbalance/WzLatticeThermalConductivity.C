#include "WzLatticeThermalConductivity.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"  





//--------------------------------------------------------//
void  WzLatticeThermalConductivity::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  ModelOptions & options = get_options ();
  _kappa_model = options.get_option("kappa_model", "constant");
  

  if (_kappa_model != "constant")
  {
   _kappa_a_x = data("therm_lat_cond_a_x", 0.0);
   _kappa_b_x = data("therm_lat_cond_b_x", 0.0);
   _kappa_c_x = data("therm_lat_cond_c_x", 0.0);
  
   _kappa_a_z = data("therm_lat_cond_a_z", 0.0);
   _kappa_b_z = data("therm_lat_cond_b_z", 0.0);
   _kappa_c_z = data("therm_lat_cond_c_z", 0.0);
  }
  else
  { 
   _kappa_x = data("therm_lat_cond_x", 0.0);
   _kappa_z = data("therm_lat_cond_z", 0.0);
  }
  

}

//---------------------------------------------------------//

void WzLatticeThermalConductivity::do_init(void)
{
 
  ModelOptions & options = get_options ();

  double k_x;
  double k_z;
 
if (_kappa_model != "constant")

{ 
 
 _kappa_a_x = options.get_option("therm_lat_cond_x", _kappa_a_x );
 _kappa_b_x = options.get_option("therm_lat_cond_x", _kappa_b_x );
 _kappa_c_x = options.get_option("therm_lat_cond_x", _kappa_c_x );
 
 _kappa_a_z = options.get_option("therm_lat_cond_x", _kappa_a_z );
 _kappa_b_z = options.get_option("therm_lat_cond_x", _kappa_b_z );
 _kappa_c_z = options.get_option("therm_lat_cond_x", _kappa_c_z );

  k_x = 1.0 / (_kappa_a_x + _kappa_b_x * _temperature +  _kappa_c_x * _temperature * _temperature );

  k_z = 1.0 / (_kappa_a_z + _kappa_b_z * _temperature +  _kappa_c_z * _temperature * _temperature );

 }
 else
 {

  k_x = options.get_option("therm_lat_cond_x",_kappa_x);
  k_z = options.get_option("therm_lat_cond_z",_kappa_z);

}


   _conductivity(1,1) = k_x;
   _conductivity(2,2) = k_x;
   _conductivity(3,3) = k_z;

  Material* mat = get_material();

  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

  rotate_to_calculation_system(cr.RotMatrix);

}


void WzLatticeThermalConductivity::update_tensor(void) 
{

 

  double k_x;
  double k_z;
 
if (_kappa_model != "constant")

{ 
 
  k_x = 1.0 / (_kappa_a_x + _kappa_b_x * _temperature +  _kappa_c_x * _temperature * _temperature );

  k_z = 1.0 / (_kappa_a_z + _kappa_b_z * _temperature +  _kappa_c_z * _temperature * _temperature );

  _conductivity(1,1) = k_x;
  _conductivity(2,2) = k_x;
  _conductivity(3,3) = k_z;


  Material* mat = get_material();

  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

  rotate_to_calculation_system(cr.RotMatrix);

 }



  
  

}
