// $Id$

#include "SchottkyContact.h"
#include "DriftDiffusionProperties.h"
#include "ModelErrorException.h"




void
SchottkyContact::do_init(void)
{
  ElectricalContact::do_init();

  double barrier = get_options().get_option("barrier_height", 0.8);
  barrier = get_options().get_option("barrier", barrier);

  std::string band = get_options().get_option("band", "c");
  if (band == "c")
    _workfunction = get_reference_material().get_conduction_band_edge() - barrier;
  else if (band == "v")
    _workfunction = get_reference_material().get_valence_band_edge() + barrier;
  else
    throw ModelErrorException("SchottkyContact: undefined band specified");

  _band = band[0];

  // we override when work function is given explicitly
  _workfunction = get_options().get_option("work_function", _workfunction);

  // TODO this is dirty, but assures the barrier in strained case
  if (_band == 'c')
    _workfunction -= get_reference_material().get_conduction_band_edge();
  else if (_band == 'v')
    _workfunction -= get_reference_material().get_valence_band_edge();

}



double
SchottkyContact::get_boundary_value(DriftDiffusionDefs::Variable variable)
{
  double val = 0.0;
  switch (variable)
  {
    case DriftDiffusionDefs::POTENTIAL:
      val = _workfunction;
      if (_band == 'c')
        val += get_reference_material().get_conduction_band_edge();
      else if (_band == 'v')
        val += get_reference_material().get_valence_band_edge();
    case DriftDiffusionDefs::FERMIE:
      break;
    case DriftDiffusionDefs::FERMIH:
      break;
  }
  
  return val;
}


