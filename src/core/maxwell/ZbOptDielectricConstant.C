#include "ZbOptDielectricConstant.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h" 



//--------------------------------------------------------//
void  ZbOptDielectricConstant::read_database(void)
{

  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

 

  if (_eps_model != "constant")
  {
    _eps = data("optical_epsilon", 1.0);
  
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
 
  if (_eps_model != "constant")

  { 
    _eps = options.get_option("optical_epsilon", 1.0);
    

  }
  else
  {
    InitFailedException("ZbOptDielectricConstant::do_init ()  incorrect model " + _eps_model);

  }


  
  _dielectric_constant_real(1,1) = _eps;

  _dielectric_constant_real(2,2) = _eps;

  _dielectric_constant_real(3,3) = _eps; 

  




}



