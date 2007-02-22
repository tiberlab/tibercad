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

  double k = data("therm_lat_cond", 0.0);

  _conductivity(1,1) = k;
  _conductivity(2,2) = k;
  _conductivity(3,3) = k;
  
}

//---------------------------------------------------------//

void ZbLatticeThermalConductivity::do_init(void)
{
  ModelOptions & options = get_options ();
  double k = options.get_option("therm_lat_cond", _conductivity(1,1) );

  _conductivity(1,1) = k;
  _conductivity(2,2) = k;
  _conductivity(3,3) = k;
  
  Material* mat = get_material();

  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

  rotate_to_calculation_system(cr.RotMatrix);
}
