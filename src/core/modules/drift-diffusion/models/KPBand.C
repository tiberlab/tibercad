// $Id$

#include "KPBand.h"
#include "DDsemiconductor.h"
#include "DriftDiffusionProperties.h"
#include "Material.h"
#include "Messages.h"


TIBER_MODULE(KPBand, band_properties, kp)


KPBand::KPBand(const ModelOptions& options) :
  BandProperties(options),
  _bulk_model(NULL)
{

}



KPBand::~KPBand(void)
{

}


void
KPBand::create_submodels(void)
{
  ModelOptions kpopts;
  if (get_options().has_submodel("kp"))
  {
    kpopts = get_options().submodels_begin("kp")->second;
  }

  _bulk_model = DDsemiconductor::create(get_material()->get_structure(), kpopts);

  if (_bulk_model == NULL)
    throw InitFailedException("Cannot create KdotP bulk model for material "
        + get_material()->get_name());

  add_submodel("kp", _bulk_model);
  get_options().delete_submodels("kp");
}



void
KPBand::do_calculate(void)
{
  assert(_bulk_model != NULL);

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  // it wants temperature in K
  _bulk_model->set_temperature(get_lattice_temperature() / Constants::k_B);
  _bulk_model->set_strain(dd.get_strain());

  // calculate conduction and valence band data
  if (get_option("particle", "el") == "el")
  {
    _bulk_model->calculate_conduction_band_extremum();

    const std::vector<DDsemiconductor::band_extremum>& cbs =
        _bulk_model->get_conduction_band_energy_mass();

    // get minimum
    int id = 0;

    for (unsigned int i = 1; i < cbs.size(); i++)
    {
      if (cbs[i].energy < cbs[id].energy)
        id = i;
    }
    band_edge() = cbs[id].energy;
    effective_mass() = cbs[id].mass_DOS
      * std::pow(cbs[id].degeneracy, 2.0 / 3.0);
    degeneracy() = cbs[id].degeneracy;

  }
  else if (get_option("particle", "el") == "hl")
  {
    _bulk_model->calculate_valence_band_extremum();

    const std::vector<DDsemiconductor::band_extremum>& vbs =
      _bulk_model->get_valence_band_energy_mass();

    // get maximum
    int id = 0;

    //double kT = SimulationOptions::T * Constants::k_B;
    double kT = dd.get_lattice_temperature();
    double delta_max = 4.0 * kT;
    for (unsigned int i = 1; i < vbs.size(); i++)
    {
      if (vbs[i].energy > vbs[id].energy)
        id = i;
    }
    band_edge() = vbs[id].energy;
    degeneracy() = vbs[id].degeneracy;

    double tmp = 0;
    // include other bands
    for (unsigned int i = 0; i < vbs.size(); i++)
    {
      double delta = get_band_edge() - vbs[i].energy;
      if (delta < delta_max)
        tmp += vbs[i].degeneracy * std::pow(vbs[i].mass_DOS, 1.5)
          * std::exp(-delta / kT);
    }
    effective_mass() = std::pow(tmp, 2.0 / 3.0);
  }
}


void
KPBand::do_print_info(void)
{
  Messages::info("(band parameters obtained from bulk k.p)");
}
