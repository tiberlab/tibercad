#include "ZbLatticeThermalConductivity.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"  
//--------------------------------------------------------//
void  ZbLatticeThermalConductivity::read_database(void)
{
 
    
  const Material* mat = get_material();

  GetPot data((mat->get_database()).get_data_file());

  _kappa = data("therm_lat_cond", 0.0);
  
 
}

//---------------------------------------------------------//



void ZbLatticeThermalConductivity::do_init(void)
{

   const ModelOptions& options = get_options();

   double k;

   k = options.get_option("therm_lat_cond",_kappa);

   _conductivity(1,1) = k;
   _conductivity(2,2) = k;
   _conductivity(3,3) = k;
  
 
}

