#include "WzOptDielectricConstant.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h" 



//--------------------------------------------------------//
void  WzOptDielectricConstant::read_database(void)
{

  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  ModelOptions & options = get_options ();

  _eps_model = options.get_option("model", "constant");

  if (_eps_model == "constant")
  {
    _eps_a = data("optical_epsilon_x", 1.0);
    _eps_c = data("optical_epsilon_z", 1.0);
  }
 
  _dielectric_constant_real(1,1) = _eps_a;
  _dielectric_constant_real(2,2) = _eps_a;

  _dielectric_constant_real(3,3) = _eps_c; 

}

//---------------------------------------------------------//

void  WzOptDielectricConstant::do_init(void)
{
 
  

  

  ModelOptions & options = get_options ();

  _eps_model = options.get_option("model", "constant");
 
  if (_eps_model == "constant")

  { 
    _eps_a = get_parameter("optical_epsilon_x", _eps_a);
    _eps_c = get_parameter("optical_epsilon_z", _eps_c);

  }
  else
  {
    InitFailedException("WzOptDielectricConstant::do_init ()  incorrect model " + _eps_model);

  }

 
  const Material* mat = get_material();

  const ModelOptions & options_mat = mat->get_options ();

  
  if (_eps_model == "constant")

  { 
    _eps_a =  options_mat.get_option("optical_epsilon_x", _eps_a);
    _eps_c =  options_mat.get_option("optical_epsilon_z", _eps_c);

  }

  


  

  _dielectric_constant_real(1,1) = _eps_a;
  _dielectric_constant_real(2,2) = _eps_a;

  _dielectric_constant_real(3,3) = _eps_c; 

  




}



