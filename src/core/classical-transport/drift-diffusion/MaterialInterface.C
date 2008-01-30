// $Id$

#include "MaterialInterface.h"
#include "Constants.h"
#include "DriftDiffusionProperties.h"



MaterialInterface::MaterialInterface(void)
  : _Ns(0.0),
    _Es(-1.0),
    _g_factor(0.5)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::NEUMANN);
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
      double denom = 1.0 + _g_factor * std::exp(arg);
      c = -_Ns / denom;
    }
    else
      c = _Ns;
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
    double kT = get_material().get_lattice_temperature();
    double Ec = get_material().get_conduction_band_edge();
    double V = get_material().get_electric_potential();
    double Efn = get_material().get_electron_electro_chemical_potential();

    double arg = -(V - Efn - Ec + _Es) / kT;
    double tmp = _g_factor * std::exp(arg);
    double denom = 1.0 + tmp;
    denom *= denom;
    
    dc[0] = -_Ns * tmp / (denom * kT);
    dc[1] = -dc[0];
  }
}



