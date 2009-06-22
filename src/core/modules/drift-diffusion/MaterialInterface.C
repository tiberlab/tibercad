// $Id$

#include "MaterialInterface.h"
#include "Constants.h"
#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"



MaterialInterface::MaterialInterface(void)
  : _Ns(0.0),
    _Es(-1.0),
    _g_factor(2.0),
    _srec(NULL)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::NEUMANN);
  is_real_contact(false);
}


MaterialInterface::~MaterialInterface(void)
{
  destroy(_srec);
}




void
MaterialInterface::do_init(void)
{
  ElectricalContact::do_init();

  _Ns = get_options().get_option("Ns", _Ns);
  _Es = get_options().get_option("Es", _Es);
  _g_factor = get_options().get_option("g", _g_factor);

  if (get_options().get_option("surface_rec", false))
  {
    double vrec_e = get_options().get_option("rec_vel_e", 1e3);
    double vrec_h = get_options().get_option("rec_vel_h", 1e3);

    ModelOptions srh_opts;
    srh_opts.set_option("tau_n", 1.0 / vrec_e);
    srh_opts.set_option("tau_p", 1.0 / vrec_h);
    _srec = RecombinationModelInterface::create("srh", srh_opts);
    if (_srec != NULL)
    {
      _srec->set_driftdiffusionproperties(&get_reference_material());
      _srec->set_material(get_reference_material().get_material());
      _srec->set_simulator_id(get_reference_material().get_simulator_id());
      _srec->init();
    }
  }
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

  }


  if (_srec != NULL)
  {
    double Re, Rh;
    _srec->set_driftdiffusionproperties(&get_material());
    _srec->get_net_recombination_rates(Re, Rh);

    if (variable == DriftDiffusionDefs::FERMIE)
      c = Re;

    if (variable == DriftDiffusionDefs::FERMIH)
      c = Rh;
  }

  if (is_internal_boundary())
    c /= 2.0;
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
    }
    else
      dc[0] = 0.0;

    dc[1] = -dc[0];
  }


  if (_srec != NULL)
  {
    std::vector<double> dRe(3), dRh(3);
    _srec->set_driftdiffusionproperties(&get_material());
    _srec->get_net_recombination_rate_derivatives(dRe, dRh);

    if (variable == DriftDiffusionDefs::FERMIE)
    {
      dc[1] = -dRe[0] * get_material().get_electron_density_derivative();
      dc[2] = -dRe[1] * get_material().get_hole_density_derivative();
      dc[0] = -(dRe[0] + dRe[1]);
    }

    if (variable == DriftDiffusionDefs::FERMIH)
    {
      dc[1] = -dRh[0] * get_material().get_electron_density_derivative();
      dc[2] = -dRh[1] * get_material().get_hole_density_derivative();
      dc[0] = -(dRh[0] + dRh[1]);
    }
  }

  if (is_internal_boundary())
  {
    dc[0] /= 2.0;
    dc[1] /= 2.0;
    dc[2] /= 2.0;
  }

}



