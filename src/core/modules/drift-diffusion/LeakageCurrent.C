// $Id$


#include "LeakageCurrent.h"
#include "Variable.h"

#include "TunnelingCurrent.h"




void
LeakageCurrent::do_init(void)
{
  ElectricalContact::do_init();

  get_parameter("A", _A);
  get_parameter("c", _c);
  _outer_voltage = get_option("outer_voltage", _outer_voltage);

  std::string tc_name = get_option("tunneling_simulation", "");
  _tc = dynamic_cast<TunnelingCurrent*>(
      SimulationInterface::find_simulation(tc_name));
}



void
LeakageCurrent::get_normal_derivative(DriftDiffusionDefs::DDVariable variable,
        double& a, double& c)
{
  a = 0.0;
  c = 0.0;

  double Vg = Variable::get_variable_value<double>(_outer_voltage);
  //double Ef = get_material().get_electron_electro_chemical_potential();
  double Ec = get_material().get_conduction_band_edge() -
    get_material().get_electric_potential();
  double Vdiff = std::abs(-Vg - Ec + 0.81);

  double I = 0.0;

  if (_tc == NULL)
  {
    I = _A * Vdiff * Vdiff * std::exp(Vdiff / _c);
    if (Vdiff < 1.0)
      I = 0.0;
  }
  else
    I = _tc->get_current(Vdiff);


  if (-Vg < Ec)
    I = -I;

  if (variable == DriftDiffusionDefs::FERMIE)
    c = - I / Constants::e;
}


