#include "ZbLatticeThermalConductivity.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"  
//--------------------------------------------------------//
void  ZbLatticeThermalConductivity::read_database(void)
{
 
  const ModelOptions& options = get_options();
  _kappa_model = options.get_option("kappa_model","constant");
    

  const Material* mat = get_material();

  GetPot data((mat->get_database()).get_data_file());

 

  if (_kappa_model != "constant")
  {

   _kappa_a = data("therm_lat_cond_a", 0.0);
   _kappa_b = data("therm_lat_cond_b", 0.0);
   _kappa_c = data("therm_lat_cond_c", 0.0);

  }
  else
  { 
   _kappa = data("therm_lat_cond", 0.0);
  
  }
  
 
}

//---------------------------------------------------------//



void ZbLatticeThermalConductivity::do_init(void)
{

   const ModelOptions& options = get_options();
   double k;

  if (_kappa_model != "constant")

  {
       

   _kappa_a = options.get_option("therm_lat_cond_a", _kappa_a);
   _kappa_b = options.get_option("therm_lat_cond_b", _kappa_b);
   _kappa_c = options.get_option("therm_lat_cond_c", _kappa_c);

   k = 1.0 / (_kappa_a + _kappa_b * _temperature +  _kappa_c * _temperature * _temperature );

  }
  else
  {

   k = options.get_option("therm_lat_cond",_kappa);

  }
   
   _conductivity(1,1) = k;
   _conductivity(2,2) = k;
   _conductivity(3,3) = k;
  
 
}


void ZbLatticeThermalConductivity::update_tensor(void) 
{


  
  double k;

  if (_kappa_model != "constant")
  {
       
   k = 1.0 / (_kappa_a + _kappa_b * _temperature +  _kappa_c * _temperature * _temperature );
  
   _conductivity(1,1) = k;
   _conductivity(2,2) = k;
   _conductivity(3,3) = k;

  }
 
   
   


}
