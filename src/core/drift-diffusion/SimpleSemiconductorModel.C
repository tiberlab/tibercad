// $Id$

#include "SimpleSemiconductorModel.h"


using namespace DriftDiffusionDefs;


SimpleSemiconductorModel::SimpleSemiconductorModel(void)
  : _is_prepared(false)
{
}


void
SimpleSemiconductorModel::prepare_element_data(void)
{
  if (!_is_prepared)
  {
    calculate_equilibrium_properties();

    _is_prepared =  true;
  }
}

void
SimpleSemiconductorModel::do_init(void)
{
  DriftDiffusionProperties::do_init();

  get_conduction_band().band_edge = get_options().get_option("Ec", 2.269);
  get_conduction_band().effective_mass =
    get_options().get_option("meff_n", 0.6);
  get_valence_band().band_edge = get_options().get_option("Ev", 1.1047);
  get_valence_band().effective_mass =
    get_options().get_option("meff_p", 1.2);
  
  permittivity = get_options().get_option("eps_r", 11.7);

  polarization(0) = get_options().get_option("Px", 0.0);
}
