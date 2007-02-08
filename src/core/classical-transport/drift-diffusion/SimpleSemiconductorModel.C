// $Id$

#include "SimpleSemiconductorModel.h"
#include "Material.h"


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

  // for the moment we read them from the materials section
  ModelOptions& opt = get_material()->get_options();

  get_conduction_band().band_edge = opt.get_option("Ec", 2.269);
  get_valence_band().band_edge = opt.get_option("Ev", 1.1047);
  get_conduction_band().effective_mass = opt.get_option("meff_n", 0.6);
  get_valence_band().effective_mass = opt.get_option("meff_p", 1.2);
  
  permittivity = opt.get_option("permittivity", 11.7);

  polarization(0) = opt.get_option("Px", 0.0);
  polarization(1) = opt.get_option("Py", 0.0);
  polarization(2) = opt.get_option("Pz", 0.0);
}
