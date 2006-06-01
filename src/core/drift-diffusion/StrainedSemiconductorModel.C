// $Id$

#include "StrainedSemiconductorModel.h"
#include "DDsemiconductor.h"

#include "elem.h"
#include "macrostrain.h"
#include "tensor.h"

#include <iostream>

using namespace DriftDiffusionDefs;

StrainedSemiconductorModel::StrainedSemiconductorModel(
    Macrostrain* strain)
  : SemiconductorModel(),
    _strain(strain),
    _ignore_strain(false)
{
}


StrainedSemiconductorModel::StrainedSemiconductorModel(
    const StrainedSemiconductorModel& model)
  : SemiconductorModel(model),
    _strain(model._strain),
    _ignore_strain(model._ignore_strain)
{
}

void
StrainedSemiconductorModel::prepare_element_data(void)
{
  if (!_ignore_strain)
  {
    assert(elem != NULL);

    const DataMap::const_iterator end = _element_data.end();
    const DataMap::const_iterator it = _element_data.find(elem);
    if (it == end)
    {
      // set strain
      get_physical_model()->set_strain(_strain->get_strain(elem, true));

      // call method of parent class
      double temp = electron_vt / Constants::k_B;
      SemiconductorModel::calculate_equilibrium_properties(BOTH,
          temp);

      // put them into _element_data
      ElementData& elem_data = _element_data[elem];
      elem_data.Ec = get_conduction_band_edge();
      elem_data.Ev = get_valence_band_edge();
      elem_data.mc = get_conduction_band_properties().effective_mass;
      elem_data.mv = get_valence_band_properties().effective_mass;
      elem_data.Ef0 = get_equilibrium_fermi_level();
      elem_data.n0 = get_equilibrium_electron_density();
      elem_data.p0 = get_equilibrium_hole_density();

      Tensor1 pol = _strain->get_piezopolarization(elem);
      elem_data.polarization(0) = pol(1); 
      elem_data.polarization(1) = pol(2); 
      elem_data.polarization(2) = pol(3); 
      polarization = elem_data.polarization;
    }
    else
    {
      const ElementData& elem_data = it->second;

      polarization = elem_data.polarization;

      get_conduction_band_properties().band_edge = elem_data.Ec; 
      get_conduction_band_properties().effective_mass = elem_data.mc; 
      get_valence_band_properties().band_edge = elem_data.Ev; 
      get_valence_band_properties().effective_mass = elem_data.mv; 

      equilibrium_fermi_level = elem_data.Ef0;
      equilibrium_electron_density = elem_data.n0;
      equilibrium_hole_density = elem_data.p0;

      // this sets the band edges and the effective DOS in the base class
      setup_band_edges();
    }

  }
  else
    SemiconductorModel::prepare_element_data();
}


void
StrainedSemiconductorModel::reset(void)
{
  DataMap::iterator begin = _element_data.begin();
  DataMap::iterator end = _element_data.end();
  _element_data.erase(begin, end);
  std::cerr << _element_data.size() << std::endl;
}
