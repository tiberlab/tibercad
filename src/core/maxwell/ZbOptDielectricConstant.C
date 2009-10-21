// $Id$

#include "ZbOptDielectricConstant.h"
#include "Database.h"



//--------------------------------------------------------//
void  ZbOptDielectricConstant::read_database(void)
{

  Database& db = get_database();
  db.set_section("permittivity");


  _eps = db.get("optical_epsilon_x", 1.0);


  _dielectric_constant_real(1,1) = _eps;
  _dielectric_constant_real(2,2) = _eps;
  _dielectric_constant_real(3,3) = _eps;

}

//---------------------------------------------------------//

void  ZbOptDielectricConstant::do_init(void)
{



get_parameter("optical_epsilon_x", _eps);


}



