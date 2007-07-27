#include "WzOptDielectricConstant.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h" 



//--------------------------------------------------------//
void  WzOptDielectricConstant::read_database(void)
{

 const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  ModelOptions & options = get_options ();
  _eps_model = options.get_option("eps_model", "constant");

 if (_eps_model != "constant")
  {
   _eps_a_x = data("eps_a_x", 0.0);
   _eps_b_x = data("eps_b_x", 0.0);
   _eps_c_x = data("eps_c_x", 0.0);
  
   _eps_a_z = data("eps_a_z", 0.0);
   _eps_b_z = data("eps_b_z", 0.0);
   _eps_c_z = data("eps_c_z", 0.0);
  }
  else
  { 
   _eps_x = data("eps_x", 0.0);
   _eps_z = data("eps_z", 0.0);
  }
  //  ...........................................

}

//---------------------------------------------------------//

void  WzOptDielectricConstant::do_init(void)
{
 ModelOptions & options = get_options ();

  double eps_x;
  double eps_z;
 
if (_eps_model != "constant")

  { 
 
//  _eps_a_x = options.get_option("therm_lat_cond_x", _eps_a_x );
//  _eps_b_x = options.get_option("therm_lat_cond_x", _eps_b_x );
//  _eps_c_x = options.get_option("therm_lat_cond_x", _eps_c_x );
 
//  _eps_a_z = options.get_option("therm_lat_cond_x", _eps_a_z );
//  _eps_b_z = options.get_option("therm_lat_cond_x", _eps_b_z );
//  _eps_c_z = options.get_option("therm_lat_cond_x", _eps_c_z );

// //  k_x = 1.0 / (_eps_a_x + _kappa_b_x * _temperature +  _kappa_c_x * _temperature * _temperature );

// //  k_z = 1.0 / (_kappa_a_z + _kappa_b_z * _temperature +  _kappa_c_z * _temperature * _temperature );

//  }
//  else
//  {

//   k_x = options.get_option("therm_lat_cond_x",_kappa_x);
//   k_z = options.get_option("therm_lat_cond_z",_kappa_z);

  }






}


void  WzOptDielectricConstant::update_tensor(void) 
{
}


