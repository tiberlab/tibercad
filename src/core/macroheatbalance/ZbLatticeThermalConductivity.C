#include "ZbLatticeThermalConductivity.h"
#include "Database.h"
#include "RotatedCrystal.h"  
//--------------------------------------------------------//
void  ZbLatticeThermalConductivity::read_database(void)
{
 

  Database& db = get_database();
  db.set_section("thermal_conductivity/constant");

  _kappa = db.get("therm_lat_cond", 0.0);

//---------------------------------------------------------//
}


void ZbLatticeThermalConductivity::do_init(void)
{

   const ModelOptions& options = get_options();

   double k;

   k = options.get_option("therm_lat_cond",_kappa);

   _conductivity(1,1) = k;
   _conductivity(2,2) = k;
   _conductivity(3,3) = k;
  
 
}

