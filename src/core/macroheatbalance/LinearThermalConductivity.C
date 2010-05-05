#include "LinearThermalConductivity.h"
#include "Database.h"
#include "RotatedCrystal.h"

using namespace std;
//--------------------------------------------------------//
void  LinearThermalConductivity::read_database(void)
{


  //Database& db = get_database();
  //db.set_section("thermal_conductivity/constant");

  //_kappa = db.get("therm_lat_cond", 0.0);

//---------------------------------------------------------//
}


void LinearThermalConductivity::do_init(void)
{

   const ModelOptions& options = get_options();


   kx0 = 0.0;
   kz0 = 0.0;
   mx  = 0.0;
   mz  = 0.0;
   z0  = 0.0;

 
   get_parameter("kx0",kx0);
   get_parameter("mx",mx);
   get_parameter("kz0",kz0);
   get_parameter("mz",mz);
   get_parameter("z0",z0);
   
 

}

void
LinearThermalConductivity::calculate(void)
{

  Point p = _elem->centroid();
  
  double x = p(0);
  double y = p(1);
  double z = p(2);


  double kx = kx0  + mx * (z-z0);
  double kz = kz0  + mz * (z-z0);

 _conductivity(1,1) = kx;
 _conductivity(2,2) = kx;
 _conductivity(3,3) = kz;


}

void LinearThermalConductivity::do_init_alloy (const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa)
{
  const LinearThermalConductivity* modA = dynamic_cast<const LinearThermalConductivity*>(comp_A);

  const LinearThermalConductivity* modB = dynamic_cast<const LinearThermalConductivity*>(comp_B);

  _kappa = alloy(modA->_kappa, modB->_kappa, xa);

}

