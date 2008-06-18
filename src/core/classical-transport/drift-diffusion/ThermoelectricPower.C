// $Id$

#include "ThermoelectricPower.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "Constants.h"


ThermoelectricPower::ThermoelectricPower(void)
  : _eQfermi(0.0),
    _hQfermi(0.0),
    _Ec(0.0),
    _Ev(0.0),
    _ElPot(0.0),
    _eFermiGrad(0.0),
    _hFermiGrad(0.0),
    _ElectricField(0.0),
    _TEmodel(CONSTANT),
    _eTEpower(0.0),
    _hTEpower(0.0)
{
}



ThermoelectricPower*
ThermoelectricPower::create(void)
{
  return new ThermoelectricPower();
}



void
ThermoelectricPower::destroy(PhysicalModelInterface* mod)
{
  delete dynamic_cast<ThermoelectricPower*>(mod);
}



ThermoelectricPower*
ThermoelectricPower::create_model(const std::string& model,
    const ModelOptions& options)
{
  return dynamic_cast<ThermoelectricPower*>(
      PhysicalModelInterface::create((create_t) create, destroy, options));
}




void
ThermoelectricPower::copy_from(const PhysicalModelInterface *rhs)
{
  const ThermoelectricPower* mod =
    dynamic_cast<const ThermoelectricPower*>(rhs);

  _eTEpower = mod->_eTEpower;
  _hTEpower = mod->_hTEpower;
}



void
ThermoelectricPower::calculate_VCA (const PhysicalModelInterface *comp_A, 
    const PhysicalModelInterface *comp_B, double xa) 
{ 
  const  ThermoelectricPower* modA =
    dynamic_cast<const ThermoelectricPower*>(comp_A);

  const ThermoelectricPower* modB =
    dynamic_cast<const ThermoelectricPower*>(comp_B);


  alloy(_eTEpower,modA->_eTEpower, modB->_eTEpower, xa);  

  alloy(_hTEpower,modA->_hTEpower, modB->_hTEpower, xa); 

}




void
ThermoelectricPower::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  _eTEpower = data("eTEpower", 0.0);
  _hTEpower = data("hTEpower", 0.0); 

}



void
ThermoelectricPower::do_init(void)
{

  std::string TEmodel = get_parameter("model", "diffusivity_model");

  if (TEmodel == "constant")
  {
 
    _TEmodel = CONSTANT;

    _eTEpower = get_parameter("eTEpower", _eTEpower);
    _hTEpower = get_parameter("hTEpower", _hTEpower);
  }
  else if (TEmodel == "diffusivity_model")
    _TEmodel = DIFFUSIVITY;
  else 
    throw InitFailedException("Unknown thermoelectric power model: " + _TEmodel);

}




void
ThermoelectricPower::calculate(void)
{
  if (_TEmodel == DIFFUSIVITY)
  {
    _eTEpower = -Constants::k_B * (5.0 / 2.0 + (_eQfermi +  _Ec - _ElPot) / _Tloc);

    _hTEpower =  Constants::k_B * (5.0 / 2.0 - (_hQfermi + _Ev - _ElPot) / _Tloc);
  }

}


void
ThermoelectricPower::calculate_derivatives(void)
{

  if (_TEmodel == DIFFUSIVITY)
  {

    _eTEpowerGrad = - Constants::k_B / _Tloc * (_eFermiGrad +  _ElectricField);

    _hTEpowerGrad = - Constants::k_B / _Tloc * (_hFermiGrad +  _ElectricField);
   
   
  }

}





