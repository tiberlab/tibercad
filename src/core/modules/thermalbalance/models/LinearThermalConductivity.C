// $Id$


#include "Material.h"
#include "LinearThermalConductivity.h"
#include "RotatedCrystal.h"


TIBER_MODULE(LinearThermalConductivity, thermal_conductivity, linear)

using namespace std;


LinearThermalConductivity::LinearThermalConductivity(const ModelOptions& options):ThermalConductivityModel(options)
{
  _kx = 0.0;
  _kz = 0.0;
}

void 
LinearThermalConductivity::read_database(void)
{

  Database& db = get_database();
  db.set_section("thermal_conductivity/constant");

  _kx = db.get("therm_lat_cond_x",0.0, true);

  if (get_material()->get_structure() == "wz")
    _kz = db.get("therm_lat_cond_z",0.0, true);
  else
    _kz = _kx;

}


void
LinearThermalConductivity::do_init(void)
{

  get_parameter("therm_lat_cond_x",_kx);
  
  if (get_material()->get_structure() == "wz")
    get_parameter("therm_lat_cond_z",_kz);
  else
    _kz = _kx;


  RealTensor kappa(0);
  kappa(0,0) = _kx;
  kappa(1,1) = _kx;
  kappa(2,2) = _kz;

  set_thermal_conductivity(kappa);

 //  // cout<<kappa<<endl;
//   if (get_material()->get_structure() == "wz")
//   {
//     const RotatedCrystal&   cr = get_material()->get_rotated_crystal();
//     rotate_to_calculation_system(cr.RotMatrix);
//   }
}


void
LinearThermalConductivity::do_init_alloy(const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa)
{
   const LinearThermalConductivity* modA = dynamic_cast<const LinearThermalConductivity*>(comp_A);
   const LinearThermalConductivity* modB = dynamic_cast<const LinearThermalConductivity*>(comp_B);

   _kx = alloy(modA->_kx, modB->_kx, xa);
   _kz = alloy(modA->_kz, modB->_kz, xa);

   //Read Alloy Database does not work because here however take the linear approximantion. TO FIX IT. 
   get_parameter("therm_lat_cond_x",_kx);
   get_parameter("therm_lat_cond_z",_kz);

}


//--------------------------------------------------------//
//void  
//ConstantThermalConductivity::read_database_alloy(void)
//{

//    Database& db = get_database();
//    //Sound Velocity
//    db.set_section("thermal_conductivity/constant");
//    _kx = db.get("therm_lat_cond_x",_kx, true);

//    if (get_material()->get_structure() == "wz")
//      _kz = db.get("therm_lat_cond_z",_kz, true);
//    else
//      _kz = _kx;

//    RealTensor kappa(0);
//    kappa(0,0) = _kx;
//    kappa(1,1) = _kx;
//    kappa(2,2) = _kz;
//    set_thermal_conductivity(kappa);
   
//    // cout<<kappa<<endl;
//    if (get_material()->get_structure() == "wz")
//    {
//      const RotatedCrystal&   cr = get_material()->get_rotated_crystal();
//      rotate_to_calculation_system(cr.RotMatrix);
//    }

// } 
