// $Id$

#include "WzOptDielectricConstant.h"
#include "Database.h"



//--------------------------------------------------------//
void  WzOptDielectricConstant::read_database(void)
{

  const Database& db = get_database();
  db.set_section("permittivity");

  ModelOptions & options = get_options ();

  _eps_a = db.get("optical_epsilon_x", 1.0);
  _eps_c = db.get("optical_epsilon_z", 1.0);
  
  _dielectric_constant_real(1,1) = _eps_a;
  _dielectric_constant_real(2,2) = _eps_a;

  _dielectric_constant_real(3,3) = _eps_c;

}

//---------------------------------------------------------//

void  WzOptDielectricConstant::do_init(void)
{



    _eps_a = get_option("optical_epsilon_x", _eps_a);
    _eps_c = get_option("optical_epsilon_z", _eps_c);


  _dielectric_constant_real(1,1) = _eps_a;
  _dielectric_constant_real(2,2) = _eps_a;

  _dielectric_constant_real(3,3) = _eps_c;



}



