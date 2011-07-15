// $Id$


#include "SemiconductorModel.h"
#include "DDsemiconductor.h"

#include "Material.h"
#include "Messages.h"


#include <sstream>


TIBER_MODULE(SemiconductorModel, ddbulk, default)



SemiconductorModel::~SemiconductorModel(void)
{
  reset();
}

SemiconductorModel::SemiconductorModel(const ModelOptions& options)
  : DriftDiffusionProperties(options),
    _is_prepared(false),
    _recompute_band_parameters(false)
{
}

void
SemiconductorModel::reset(void)
{
  DataMap::iterator begin = get_data_map().begin();
  DataMap::iterator end = get_data_map().end();
  _element_data.erase(begin, end);
}




void
SemiconductorModel::do_init(void)
{

  Parent::do_init();

  _recompute_band_parameters = get_option("recompute_band_parameters",
      _recompute_band_parameters);

}





void
SemiconductorModel::prepare_element_data(void)
{
  const Elem* elem = get_element();
  assert(elem != NULL);

  if (is_inhomogeneous())
  {
    const DataMap::const_iterator end = get_data_map().end();
    const DataMap::const_iterator it = get_data_map().find(elem);
    if ((it == end) || _recompute_band_parameters)
    {
      ElementData& elem_data = get_data_map()[elem];

      calculate_equilibrium_properties();

      elem_data.Ec = get_conduction_band_edge();
      elem_data.Ev = get_valence_band_edge();
      elem_data.mc = get_conduction_band().get_effective_mass();
      elem_data.mv = get_valence_band().get_effective_mass();
      elem_data.Ef0 = get_equilibrium_fermi_level();
      //elem_data.ni = get_intrinsic_density();

    }
    else
    {
      const ElementData& elem_data = it->second;

      get_conduction_band()._band_edge = elem_data.Ec;
      get_conduction_band()._effective_mass = elem_data.mc;
      get_valence_band()._band_edge = elem_data.Ev;
      get_valence_band()._effective_mass = elem_data.mv;

      // this sets the band edges and the effective DOS in the base class
      setup_band_edges();

      set_equilibrium_properties(elem_data.Ef0);
    }
  }
  else
  {
    if (!_is_prepared || _recompute_band_parameters)
    {
      calculate_equilibrium_properties();
      _is_prepared = true;
    }
  }
}
