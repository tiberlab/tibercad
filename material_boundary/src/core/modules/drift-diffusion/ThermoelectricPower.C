// $Id$

#include "ThermoelectricPower.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "Constants.h"


//TIBER_MODULE(ThermoelectricPower, default)

namespace
{
  TiberModelObject* thpow_create(void)
  {
    return new ThermoelectricPower();
  }

  void thpow_destroy(TiberModelObject* mod)
  {
    delete mod;
  }
}


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
ThermoelectricPower::create_model(const std::string& model,
    const ModelOptions& options)
{
  return static_cast<ThermoelectricPower*>(
      PhysicalModelInterface::create(thpow_create, thpow_destroy, options));
}





void
ThermoelectricPower::do_init_alloy (const PhysicalModelInterface *comp_A,
    const PhysicalModelInterface *comp_B, double xa)
{
  const  ThermoelectricPower* modA =
    dynamic_cast<const ThermoelectricPower*>(comp_A);

  const ThermoelectricPower* modB =
    dynamic_cast<const ThermoelectricPower*>(comp_B);

  _TEmodel = modA->_TEmodel;

  alloy(_eTEpower,modA->_eTEpower, modB->_eTEpower, xa);

  alloy(_hTEpower,modA->_hTEpower, modB->_hTEpower, xa);

}




void
ThermoelectricPower::read_database(void)
{


  Database& db = get_database();
  db.set_section("thermoelectric_power/constant");

  _eTEpower = db.get("eTEpower", 0.0);
  _hTEpower = db.get("hTEpower", 0.0);

}



void
ThermoelectricPower::do_init(void)
{

  const std::string& TEmodel = get_option("model", "diffusivity_model");

  if (TEmodel == "constant")
  {

    _TEmodel = CONSTANT;

    get_parameter("eTEpower", _eTEpower);
    get_parameter("hTEpower", _hTEpower);
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





