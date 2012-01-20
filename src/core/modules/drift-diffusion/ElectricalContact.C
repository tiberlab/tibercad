// $Id$

#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"
#include "ParticleDensity.h"
#include "Initializer.h"


ElectricalContact::ElectricalContact(const ModelOptions& options)
  : DDInterfaceModel(options),
    _voltage(0.0),
    _surfres(0.0),
    _contact_fermilevel(0.0),
    _vrec_n(0),
    _vrec_p(0),
    _fixed_vrec_n(false),
    _fixed_vrec_p(false)
{
  // this is a real contact
  has_current(true);
}


void
ElectricalContact::do_init(void)
{

  if (get_option("zero_field", false))
    set_type(0, NEUMANN);

  DDInterfaceModel::do_init();

  get_parameter("rec_velocity_n", _vrec_n);
  get_parameter("rec_velocity_p", _vrec_p);

  if (get_option("zero_grad_fermi_e", false) || (_vrec_n > 0))
  {
    set_type(1, NEUMANN);
    _fixed_vrec_n = true;
  }
  if (get_option("zero_grad_fermi_h", false) || (_vrec_p > 0))
  {
    set_type(2, NEUMANN);
    _fixed_vrec_p = true;
  }


  get_parameter("contact_resistance", _surfres);

  get_parameter("voltage", _voltage);
}


void
ElectricalContact::do_compute(void)
{
  if (get_type(0) != NEUMANN)
    coeff_g(0) = _contact_fermilevel + get_inner_voltage();

  if (get_type(1) != NEUMANN)
    coeff_g(1) = get_inner_voltage();
  else
  {
    const DriftDiffusionProperties::PointData& pd =
        get_dd_properties()->get_point_data();

    // we have to take the correct equilibrium density!
    ParticleDensity& el = get_dd_properties()->get_electrons();
    BandProperties& cb = get_dd_properties()->get_conduction_band();
    cb.set_temperature(pd.electron_vt);
    el.set_classical_parameters(cb.get_effective_DOS(),
        get_dd_properties()->get_conduction_band_edge() - _contact_fermilevel, 0,
        pd.electron_vt);
    double n0 = el.get_particle_density();

    double Rn = _vrec_n * (pd.electron_density - n0);
    double dRn = _vrec_n * pd.electron_density_derivative;

    coeff_g(1) += Rn;
    jacobian(1, 0) += dRn;
    jacobian(1, 1) -= dRn;
  }

  if (get_type(2) != NEUMANN)
    coeff_g(2) = get_inner_voltage();
  else
  {
    const DriftDiffusionProperties::PointData& pd =
        get_dd_properties()->get_point_data();

    // we have to take the correct equilibrium density!
    ParticleDensity& hl = get_dd_properties()->get_holes();
    BandProperties& vb = get_dd_properties()->get_valence_band();
    vb.set_temperature(pd.hole_vt);
    hl.set_classical_parameters(vb.get_effective_DOS(),
        _contact_fermilevel - get_dd_properties()->get_valence_band_edge(), 0,
        pd.hole_vt);
    double p0 = hl.get_particle_density();

    double Rp = _vrec_p * (pd.hole_density - p0);
    double dRp = _vrec_p * pd.hole_density_derivative;

    coeff_g(2) += Rp;
    jacobian(2, 0) += dRp;
    jacobian(2, 2) -= dRp;
  }

}

/*
double
ElectricalContact::get_contact_voltage_drop(void) const
{
  double j = Constants::e * (get_normal_hole_flux() -
      get_normal_electron_flux());

  // a negative current means inflowing current
  return -_surfres * j;
}
*/



