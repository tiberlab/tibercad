// $Id: SchottkyContact.C 1912 2010-04-16 10:13:35Z maufder $

#include "SchottkyContact.h"
#include "DriftDiffusionProperties.h"
#include "ModelErrorException.h"



TIBER_MODULE(SchottkyContact, ddbnd, schottky)


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
    _workfunction = get_dd_properties()->get_conduction_band_edge() - barrier;
  else if (band == "v")
    _workfunction = get_dd_properties()->get_valence_band_edge() + barrier;
  else
    throw ModelErrorException("SchottkyContact: undefined band specified");

  _band = band[0];

  // we override when work function is given explicitly
  if (has_option("work_function"))
  {
    // if the workfunction is set then we assume that this one should
    // be fixed and not the resulting barrier
    _fixed_barrier = false;
    _workfunction = get_option("work_function", _workfunction);
  }

  _fixed_barrier = get_option("fixed_barrier", _fixed_barrier);

  if (_fixed_barrier)
  {
    // this is a dirty trick, but assures the barrier in strained case
    if (_band == 'c')
      _workfunction -= get_dd_properties()->get_conduction_band_edge();
    else
      _workfunction -= get_dd_properties()->get_valence_band_edge();
  }

  _thermionic_emission = get_option("thermionic_emission", false);
  if (_thermionic_emission)
  {
    double m = 0.1 * Constants::me;
    double vth = SimulationOptions::T * Constants::e / (2 * M_PI * m);
    set_recombination_velocities(vth, -1);
  }

}


void
SchottkyContact::do_compute(void)
{
  double val = _workfunction;
  if (_fixed_barrier)
  {
    if (_band == 'c')
      val += get_dd_properties()->get_conduction_band_edge();
    else
      val += get_dd_properties()->get_valence_band_edge();
  }
  set_barrier(val);

  ElectricalContact::do_compute();
}



