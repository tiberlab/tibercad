// $Id$

#include "StrainedSemiconductorModel.h"
#include "DDsemiconductor.h"
#include "Macrostrain.h"
//#include "tensor.h"
#include "Utils.h"

#include "elem.h"

#include <iostream>


TIBER_MODULE(StrainedSemiconductorModel, strained)



using namespace DriftDiffusionDefs;

StrainedSemiconductorModel::StrainedSemiconductorModel(void)
  : strain_model_(NULL),
    ignore_strain_(false)
{
}



void
StrainedSemiconductorModel::prepare_element_data(void)
{
  if (!ignore_strain_)
  {
    const Elem* elem = get_element();
    assert(elem != NULL);

    const DataMap::const_iterator end = element_data_.end();
    const DataMap::const_iterator it = element_data_.find(elem);
    if (it == end)
    {
      // set strain
      set_strain(strain_model_->get_strain_crystal(elem));
      get_physical_model()->set_strain(get_strain());

      // call method of parent class
      SemiconductorModel::calculate_equilibrium_properties();

      // put them into element_data_
      ElementData& elem_data = element_data_[elem];
      elem_data.Ec = get_conduction_band_edge();
      elem_data.Ev = get_valence_band_edge();
      elem_data.mc = get_conduction_band().effective_mass;
      elem_data.mv = get_valence_band().effective_mass;
      elem_data.Ef0 = get_equilibrium_fermi_level();
      elem_data.ni = get_intrinsic_density();

      Tensor1 pol =
        strain_model_->get_built_in_polarization(elem, elem->centroid());
      elem_data.polarization(0) = pol(1); 
      elem_data.polarization(1) = pol(2); 
      elem_data.polarization(2) = pol(3); 
      polarization += elem_data.polarization;
    }
    else
    {
      const ElementData& elem_data = it->second;

      polarization += elem_data.polarization;

      get_conduction_band().band_edge = elem_data.Ec; 
      get_conduction_band().effective_mass = elem_data.mc; 
      get_valence_band().band_edge = elem_data.Ev; 
      get_valence_band().effective_mass = elem_data.mv; 

      equilibrium_fermi_level = elem_data.Ef0;
      intrinsic_density = elem_data.ni;

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
  DataMap::iterator begin = element_data_.begin();
  DataMap::iterator end = element_data_.end();
  element_data_.erase(begin, end);
}

void
StrainedSemiconductorModel::copy_from(const PhysicalModelInterface* rhs)
{
  SemiconductorModel::copy_from(rhs);

  const StrainedSemiconductorModel* sc =
    dynamic_cast<const StrainedSemiconductorModel*>(rhs);

  strain_model_ = sc->strain_model_;
  ignore_strain_ = sc->ignore_strain_;
}


void
StrainedSemiconductorModel::do_init(void)
{
  SemiconductorModel::do_init();

  
  std::string strain_sim =
    get_options().get_option("strain_simulation",
        Utils::extract_typename(typeid(strain_model_)));
  
  // find the strain calculation to use
  strain_model_ = dynamic_cast<Macrostrain*>(
      SimulationInterface::find_simulation(strain_sim));

  if (strain_model_ == NULL)
  {
    std::string msg("Simulation "+std::string(strain_sim)+" not found");
    throw InitFailedException(msg);
  }

  ignore_strain_ = get_options().get_option("ignore_strain", false);
}
