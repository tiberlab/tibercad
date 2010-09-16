// $Id$


#include "Material.h"
#include "Database.h"
#include "ConstantThermalConductivity.h"
#include "RotatedCrystal.h"


TIBER_MODULE(ConstantThermalConductivity, thermal_conductivity, constant)

using namespace std;


ConstantThermalConductivity::ConstantThermalConductivity(const ModelOptions& options):ThermalConductivityModel(options)
{
 _kappa.resize(3,0.0);
}

void
ConstantThermalConductivity::read_database(void)
{

  Database& db = get_database();
  db.set_section("thermal_conductivity/constant");
  db.get("ThermCond", _kappa);

}

void
ConstantThermalConductivity::do_init(void)
{
  get_parameter("ThermCond",_kappa,true);



  RealTensor kappa_tens(0);
  kappa_tens(0,0) = _kappa[0];
  kappa_tens(1,1) = _kappa[1];
  kappa_tens(2,2) = _kappa[2];

  set_thermal_conductivity(kappa_tens);


  if (get_material()->get_structure() == "wz")
  {
    const RotatedCrystal&   cr = get_material()->get_rotated_crystal();
    rotate_to_calculation_system(cr.RotMatrix);
  }

}


//void
// ConstantThermalConductivity::do_init_alloy(const PhysicalModelInterface *comp_A,
//                                                 const PhysicalModelInterface *comp_B, double xa)
// {
//    const ConstantThermalConductivity* modA = dynamic_cast<const ConstantThermalConductivity*>(comp_A);
//    const ConstantThermalConductivity* modB = dynamic_cast<const ConstantThermalConductivity*>(comp_B);

//    _kappa[0] = alloy(modA->_kappa[0], modB->_kappa[0], xa);
//    _kappa[1] = alloy(modA->_kappa[1], modB->_kappa[1], xa);
//    _kappa[2] = alloy(modA->_kappa[2], modB->_kappa[2], xa);

//    // _kz = alloy(modA->_kz, modB->_kz, xa);

//    //Read Alloy Database does not work because here however take the linear approximantion. TO FIX IT.
//    //get_parameter("therm_lat_cond_x",_kx);
//    //get_parameter("therm_lat_cond_z",_kz);

// }


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
