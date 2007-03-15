#include "ThermoelectricPower.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"



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
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

   _eTEpower = data("eTEpower", 0.0);
   _hTEpower = data("hTEpower", 0.0); 

  
}

//---------------------------------------------------------//

void ThermoelectricPower::do_init(void)
{
  ModelOptions & options = get_options ();
  _eTEpower = options.get_option("eTEpower", _eTEpower );
  _hTEpower = options.get_option("hTEpower", _hTEpower );

  _eTEpower *= 0.001;
  _hTEpower *= 0.001;
  
}







//-------------------------------------------------------------------------//
