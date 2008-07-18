// $Id$

#include "MaterialInterface.h"
#include "Constants.h"
#include "DriftDiffusionProperties.h"



MaterialInterface::MaterialInterface(void)
  : _Ns(0.0),
    _Es(-1.0),
    _g_factor(2.0)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::NEUMANN);
  is_real_contact(false);
}




void
MaterialInterface::do_init(void)
{
  ElectricalContact::do_init();

  _Ns = get_options().get_option("Ns", _Ns);
  _Es = get_options().get_option("Es", _Es);
  _g_factor = get_options().get_option("g", _g_factor);
}



void
MaterialInterface::get_normal_derivative(DriftDiffusionDefs::Variable variable,
        double& a, double& c)
{
  a = 0.0;
  c = 0.0;

  if (variable == DriftDiffusionDefs::POTENTIAL)
  {
    if (_Es > 0.0)
    {
      double kT = get_material().get_lattice_temperature();
      double Ec = get_material().get_conduction_band_edge();
      double V = get_material().get_electric_potential();
      double Efn = get_material().get_electron_electro_chemical_potential();

      double arg = -(V - Efn - Ec + _Es) / kT;
      double denom = 1.0 + std::exp(arg) / _g_factor;
      c = -_Ns / denom;
    }
    else
      c = _Ns;

    if (is_internal_boundary())
      c /= 2.0;
  }
}


void
MaterialInterface::get_derivatives_of_normal_derivative(
        DriftDiffusionDefs::Variable variable,
        std::vector<double>& da, std::vector<double>& dc)
{
  da = std::vector<double>(3, 0.0);
  dc = std::vector<double>(3, 0.0);

  if (variable == DriftDiffusionDefs::POTENTIAL)
  {
    if (_Es > 0.0)
    {
      double kT = get_material().get_lattice_temperature();
      double Ec = get_material().get_conduction_band_edge();
      double V = get_material().get_electric_potential();
      double Efn = get_material().get_electron_electro_chemical_potential();

      double arg = -(V - Efn - Ec + _Es) / kT;
      double tmp = std::exp(arg) / _g_factor;
      double denom = 1.0 + tmp;
      denom *= denom;

      dc[0] = -_Ns * tmp / (denom * kT);
      if (is_internal_boundary())
        dc[0] /= 2.0;
    }
    else
      dc[0] = 0.0;

    dc[1] = -dc[0];
  }
}



