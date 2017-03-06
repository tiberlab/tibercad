// $Id: ThermoelectricPower.C 3599 2013-05-07 08:38:43Z maufder $

#include "ThermoelectricPower.h"
#include "Material.h"
#include "Database.h"
#include "Constants.h"


//TIBER _ MODULE(ThermoelectricPower, default)

namespace
{
  TiberModelObject* thpow_create(const ModelOptions& options, const void*)
  {
    return new ThermoelectricPower(options);
  }

  void thpow_destroy(TiberModelObject* mod)
  {
    delete mod;
  }
}


ThermoelectricPower::ThermoelectricPower(const ModelOptions& options)
  : DriftDiffusionModelInterface(options),
    _eQfermi(0.0),
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
ThermoelectricPower::create_model(const std::string&,
    const Material* mat, const ModelOptions& options)
{
  return static_cast<ThermoelectricPower*>(
      PhysicalModelInterface::create(thpow_create, thpow_destroy, mat, options));
}







void
ThermoelectricPower::read_database(void)
{


  const Database& db = get_database();
  db.set_section("thermoelectric_power/constant");

  _eTEpower = db.get("eTEpower", 0.0);
  _hTEpower = db.get("hTEpower", 0.0);

}



void
ThermoelectricPower::do_init(void)
{

  std::string TEmodel = get_option("type", "diffusivity_model");
  TEmodel = get_option("model", TEmodel);

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





