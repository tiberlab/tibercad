// $Id$


#include "LeakageCurrent.h"

#include "TunnelingCurrent.h"




void
LeakageCurrent::do_init(void)
{
  ElectricalContact::do_init();

  _A = get_options().get_option("A", _A);
  _c = get_options().get_option("c", _c);
  _outer_voltage = get_options().get_option("outer_voltage", _outer_voltage);

  std::string tc_name = get_options().get_option("tunneling_simulation", "");
  _tc = dynamic_cast<TunnelingCurrent*>(
      SimulationInterface::find_simulation(tc_name));
}



void
LeakageCurrent::get_normal_derivative(DriftDiffusionDefs::Variable variable,
        double& a, double& c)
{
  a = 0.0;
  c = 0.0;

  double Vg = Variable::get_variable_value(_outer_voltage);
  //double Ef = get_material().get_electron_electro_chemical_potential();
  double Ec = get_material().get_conduction_band_edge() -
    get_material().get_electric_potential();
  double Vdiff = std::abs(-Vg - Ec + 0.81);

  double I = 0.0;

  if (_tc == NULL)
    I = _A * Vdiff * Vdiff * std::exp(Vdiff / _c);
  else
    I = _tc->get_current(Vdiff);

  if (-Vg < Ec)
    I = -I;

  if (variable == DriftDiffusionDefs::FERMIE)
    c = - I / Constants::e;
}


