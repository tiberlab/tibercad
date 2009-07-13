// $Id$

#include "SchottkyContact.h"
#include "DriftDiffusionProperties.h"
#include "ModelErrorException.h"




void
SchottkyContact::do_init(void)
{
  ElectricalContact::do_init();

  double barrier = get_option("barrier_height", 0.8);
  barrier = get_option("barrier", barrier);


  std::string band = get_options().get_option("band", "c");
  if (band == "c")
    _workfunction = get_reference_material().get_conduction_band_edge() - barrier;
  else if (band == "v")
    _workfunction = get_reference_material().get_valence_band_edge() + barrier;
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
      _workfunction -= get_reference_material().get_conduction_band_edge();
    else if (_band == 'v')
      _workfunction -= get_reference_material().get_valence_band_edge();
  }

}



double
SchottkyContact::get_boundary_value(DriftDiffusionDefs::DDVariable variable)
{
  double val = 0.0;
  switch (variable)
  {
    case DriftDiffusionDefs::POTENTIAL:
      val = _workfunction;
      if (_fixed_barrier)
      {
        if (_band == 'c')
          val += get_reference_material().get_conduction_band_edge();
        else if (_band == 'v')
          val += get_reference_material().get_valence_band_edge();
      }
    case DriftDiffusionDefs::FERMIE:
      break;
    case DriftDiffusionDefs::FERMIH:
      break;
  }

  return val;
}


