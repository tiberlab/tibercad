// $Id$

#include "ZbOptDielectricConstant.h"
#include "Database.h"



//--------------------------------------------------------//
void  ZbOptDielectricConstant::read_database(void)
{

  Database& db = get_database();
  db.set_section("permittivity");
 

  if (_eps_model == "constant")
  {
    _eps = db.get("optical_epsilon", 1.0);
  
  }
 
  _dielectric_constant_real(1,1) = _eps;
  _dielectric_constant_real(2,2) = _eps;
  _dielectric_constant_real(3,3) = _eps; 

}

//---------------------------------------------------------//

void  ZbOptDielectricConstant::do_init(void)
{
  ModelOptions & options = get_options ();

  _eps_model = options.get_option("model", "constant");
 
  if (_eps_model == "constant")

  { 
    _eps = get_parameter("optical_epsilon", _eps);
    

  }
  else
  {
    InitFailedException("ZbOptDielectricConstant::do_init ()  incorrect model " + _eps_model);

  }


}



