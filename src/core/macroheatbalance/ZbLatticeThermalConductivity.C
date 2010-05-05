#include "ZbLatticeThermalConductivity.h"
#include "Database.h"
#include "RotatedCrystal.h"
//--------------------------------------------------------//
void  ZbLatticeThermalConductivity::read_database(void)
{


  Database& db = get_database();
  db.set_section("thermal_conductivity/constant");

  _kappa = db.get("therm_lat_cond_x", 0.0);

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


void ZbLatticeThermalConductivity::do_init_alloy (const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa)
{
  const ZbLatticeThermalConductivity* modA = dynamic_cast<const ZbLatticeThermalConductivity*>(comp_A);

  const ZbLatticeThermalConductivity* modB = dynamic_cast<const ZbLatticeThermalConductivity*>(comp_B);

  _kappa = alloy(modA->_kappa, modB->_kappa, xa);

}

