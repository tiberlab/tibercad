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

  double K11 = data("cond_term_lat_x", 0.0);
  double K33 = data("cond_term_lat_z", 0.0);

  _conductivity(1,1) = K11;
  _conductivity(2,2) = K11;
  _conductivity(3,3) = K33;
  
}

//---------------------------------------------------------//

void WzLatticeThermalConductivity::do_init(void)
{
  ModelOptions & options = get_options ();
  double K11 = options.get_option("cond_term_lat_x", _conductivity(1,1) );
  double K33 = options.get_option("cond_term_lat_z", _conductivity(3,3) ); 

  _conductivity(1,1) = K11;
  _conductivity(2,2) = K11;
  _conductivity(3,3) = K33;
  
  Material* mat = get_material();

  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

  rotate_to_calculation_system(cr.RotMatrix);
}
