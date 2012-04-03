// $Id: SchottkyContact.C 1912 2010-04-16 10:13:35Z maufder $

#include "SchottkyContact.h"
#include "DriftDiffusionProperties.h"
#include "ModelErrorException.h"

#include "TiberModule.h"



inline
SchottkyContact::SchottkyContact(const ModelOptions& options)
  : ElectricalContact(options),
    _band('c'),
    _fixed_barrier(true)
{
  set_type(0, DIRICHLET);
  set_type(1, DIRICHLET);
  set_type(2, DIRICHLET);
}


void
SchottkyContact::do_init(void)
{

  ElectricalContact::do_init();

  double barrier = get_option("barrier_height", 0.8);
  barrier = get_option("barrier", barrier);

  std::string band = get_options().get_option("band", "c");
  if (band == "c")
    _metal_Ef = get_dd_properties()->get_conduction_band_edge() - barrier;
  else if (band == "v")
    _metal_Ef = get_dd_properties()->get_valence_band_edge() + barrier;
  else
    throw ModelErrorException("SchottkyContact: undefined band specified");

  _band = band[0];

  // we override when work function is given explicitly
  if (has_option("work_function"))
  {
    // if the workfunction is set then we assume that this one should
    // be fixed and not the resulting barrier
    _fixed_barrier = false;
    _metal_Ef = -get_option("work_function", _metal_Ef);
  }

  if (has_option("metal_fermilevel"))
  {
    // sanity check
    if (has_option("work_function"))
      throw ModelErrorException("You may not specify both work_function "
          "and metal_fermilevel for Schottky contact.");

    // if the workfunction is set then we assume that this one should
    // be fixed and not the resulting barrier
    _fixed_barrier = false;
    _metal_Ef = get_option("metal_fermilevel", _metal_Ef);
  }

  _fixed_barrier = get_option("fixed_barrier", _fixed_barrier);

  if (_fixed_barrier)
  {
    // this is a dirty trick, but assures the barrier in strained case
    if (_band == 'c')
      _metal_Ef -= get_dd_properties()->get_conduction_band_edge();
    else
      _metal_Ef -= get_dd_properties()->get_valence_band_edge();
  }

  _thermionic_emission = get_option("thermionic_emission", true);
  if (_thermionic_emission)
  {
    double temp = Constants::k_B * SimulationOptions::T;
    double vth_n = get_dd_properties()->get_conduction_band().get_thermal_velocity(temp);
    double vth_p = get_dd_properties()->get_valence_band().get_thermal_velocity(temp);

    // the effective emission velocity is fac * v_thermal

    const double fac = 0.23032943; // std::sqrt(1.0 / (6.0 * M_PI));
    set_recombination_velocities(fac * vth_n, fac * vth_p);
  }

  _tunneling = get_option("tunneling", true);

}


void
SchottkyContact::do_compute(void)
{
  double val = _metal_Ef;
  if (_fixed_barrier)
  {
    if (_band == 'c')
      val += get_dd_properties()->get_conduction_band_edge();
    else
      val += get_dd_properties()->get_valence_band_edge();
  }
  set_contact_fermilevel(val);

  if (_thermionic_emission)
  {
    const DriftDiffusionProperties* dd = get_dd_properties();
    const DriftDiffusionProperties::PointData& pd = dd->get_point_data();

    double vth_n = dd->get_conduction_band().get_thermal_velocity(pd.electron_vt);
    double vth_p = dd->get_valence_band().get_thermal_velocity(pd.hole_vt);

    const double fac = 0.23032943; // std::sqrt(1.0 / (6.0 * M_PI));
    set_recombination_velocities(fac * vth_n, fac * vth_p);
  }

  if (_tunneling)
  {

  }

  ElectricalContact::do_compute();
}



