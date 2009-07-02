#include "WzLatticeThermalConductivity.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"





//--------------------------------------------------------//
void  WzLatticeThermalConductivity::read_database(void)
{
  //  const Material* mat = get_material();
  //  GetPot data((mat->get_database()).get_data_file());


  Database& db = get_database();
  db.set_section("thermal_conductivity/constant");


 _kappa_x = db.get("therm_lat_cond_x", 0.0);
 _kappa_z = db.get("therm_lat_cond_z", 0.0);







}

//---------------------------------------------------------//

void WzLatticeThermalConductivity::do_init(void)
{

  ModelOptions & options = get_options ();

  double k_x;
  double k_z;

  k_x = options.get_option("therm_lat_cond_x",_kappa_x);
  k_z = options.get_option("therm_lat_cond_z",_kappa_z);

  _conductivity(1,1) = k_x;
  _conductivity(2,2) = k_x;
  _conductivity(3,3) = k_z;

  Material* mat = get_material();

  const RotatedCrystal&   cr = mat->get_rotated_crystal ();

  rotate_to_calculation_system(cr.RotMatrix);


}



void WzLatticeThermalConductivity::do_init_alloy (const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa)
{
  const WzLatticeThermalConductivity* modA = dynamic_cast<const WzLatticeThermalConductivity*>(comp_A);

  const WzLatticeThermalConductivity* modB = dynamic_cast<const WzLatticeThermalConductivity*>(comp_B);

  _kappa_x = alloy(modA->_kappa_x, modB->_kappa_x, xa);
  _kappa_z = alloy(modA->_kappa_z, modB->_kappa_z, xa);

}
