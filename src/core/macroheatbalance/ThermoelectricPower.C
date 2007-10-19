#include "ThermoelectricPower.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "Constants.h"



//-------------------------------------------------------------------------//

void ThermoelectricPower::copy_from(const PhysicalModelInterface *rhs)
{
  const ThermoelectricPower* mod = dynamic_cast<const  ThermoelectricPower*> (rhs);

  _eTEpower = mod->_eTEpower;
  _hTEpower = mod->_hTEpower;
}

//-------------------------------------------------------------------------//


void  ThermoelectricPower::calculate_VCA (const PhysicalModelInterface *comp_A, 
                                                const PhysicalModelInterface *comp_B, double xa) 
{ 
  const  ThermoelectricPower* modA = dynamic_cast<const  ThermoelectricPower*>(comp_A);

  const ThermoelectricPower* modB = dynamic_cast<const  ThermoelectricPower*>(comp_B);


   alloy(_eTEpower,modA->_eTEpower, modB->_eTEpower, xa);  

   alloy(_hTEpower,modA->_hTEpower, modB->_hTEpower, xa); 
  
}

//-------------------------------------------------------------------------//

//--------------------------------------------------------//
void  ThermoelectricPower::read_database(void)
{

  const ModelOptions& options = get_options();
  _TEmodel = options.get_option("model","constant");



  if  (_TEmodel == "constant")
  {
    const Material* mat = get_material();
    GetPot data((mat->get_database()).get_data_file());
    
    _eTEpower = data("eTEpower", 0.0);
    _hTEpower = data("hTEpower", 0.0); 

  }	
  
}

//---------------------------------------------------------//

void ThermoelectricPower::do_init(void)
{

     
  if (_TEmodel=="constant")
  {
    
    ModelOptions& options = get_options ();
    
    _eTEpower = options.get_option("eTEpower", _eTEpower );
        
    _hTEpower = options.get_option("hTEpower", _hTEpower );

    _eTEpower = -std::abs(_eTEpower);

      
    _hTEpower = -std::abs(_hTEpower); 
   
    
  }
  else if  (_TEmodel.compare("diffusivity_model") == 0 ) 
  {
    
    _eTEpower = - Constants::k_B * (5.0 / 2.0 + _e_mobility_term +  (_eQfermi +  _Ec) / (Constants::k_B * _Tloc)   );
    
    _hTEpower =  Constants::k_B * (5.0 / 2.0 + _h_mobility_term -  (_hQfermi + _Ev) / (Constants::k_B * _Tloc)  );

    _eTEpower = -std::abs(_eTEpower);
      
    _hTEpower = -std::abs(_hTEpower);

    
  }
  else
  {


  throw InitFailedException("Unrecongnition the thermoelectric power model:  " + _TEmodel);
  

  } 

  
  
}

void ThermoelectricPower::re_init(void)
{
 
  this->do_init();
 

  
}







//-------------------------------------------------------------------------//
