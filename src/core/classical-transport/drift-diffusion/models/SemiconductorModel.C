// $Id$


#include "SemiconductorModel.h"
#include "ZbDDsemiconductor.h"
#include "WzDDsemiconductor.h"

#include "Alloy.h"

//#include "Database.h"

//#include "getpot.h"

#include <iostream>


TIBER_MODULE(SemiconductorModel, unstrained)


using namespace DriftDiffusionDefs;
using namespace std;


SemiconductorModel::~SemiconductorModel(void)
{
  PhysicalModelInterface::destroy(bulk_model_);
}

SemiconductorModel::SemiconductorModel(void)
  : bulk_model_(NULL),
    is_prepared_(false)
{
}


void
SemiconductorModel::do_init(void)
{
  Parent::do_init();

  const ModelOptions& opt = get_options();

  PhysicalModelInterface::destroy(bulk_model_);

  bulk_model_ = DDsemiconductor::create(get_material()->get_structure(), opt);

  if (bulk_model_ == NULL)
    throw InitFailedException("Unknown structure for DDsemiconductor");

  bulk_model_->set_material(get_material());

  bulk_model_->init();
  
}



void
SemiconductorModel::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  DriftDiffusionProperties::calculate_VCA(comp_A, comp_B, xa);

  //double (*alloy)(double, double, double, double) =
  //  Alloy::calculate_VCA_parameter;

  const SemiconductorModel* scA =
    dynamic_cast<const SemiconductorModel*>(comp_A);
  const SemiconductorModel* scB =
    dynamic_cast<const SemiconductorModel*>(comp_B);

  bulk_model_->build_alloy(scA->bulk_model_, scB->bulk_model_, xa);

  permittivity = alloy(scA->permittivity, scB->permittivity,
        xa, 0.0);

}



void
SemiconductorModel::prepare_element_data(void)
{
  if (!is_prepared_)
  {
    try
    {
      calculate_equilibrium_properties();
    }
    catch (...)
    {
    }
  
    is_prepared_ = true;
  }
}


//
// TODO
// Very crude implementation at the moment
// 
void
SemiconductorModel::extract_band_properties(void)
{
  // treat conduction band
  const std::vector<DDsemiconductor::band_extremum>& cbs =
    bulk_model_->get_conduction_band_energy_mass();
  // get minimum
  int id = 0;
  for (int i = 1; i < cbs.size(); i++)
  {
    if (cbs[i].energy < cbs[id].energy)
      id = i;
  }
  get_conduction_band().band_edge = cbs[id].energy;
  get_conduction_band().effective_mass = cbs[id].mass_DOS
    * std::pow(cbs[id].degeneracy, 2.0 / 3.0);
  
  // treat valence band
  const std::vector<DDsemiconductor::band_extremum>& vbs =
    bulk_model_->get_valence_band_energy_mass();
  cerr << "extr" << vbs.size() << "\n";

  // get maximum
  id = 0;
  // TODO should be local temperature
  double kT = SimulationOptions::T * Constants::k_B;
  double delta_max = 4.0 * kT;
  for (int i = 1; i < vbs.size(); i++)
  {
    cerr << i << "    " << vbs.size() << "         " << vbs[i].energy << "\n";

    if (vbs[i].energy > vbs[id].energy)
      id = i;
  }
  get_valence_band().band_edge = vbs[id].energy;
  
  double tmp = 0;
  // include other bands
  for (int i = 0; i < vbs.size(); i++)
  {
    double delta = get_valence_band().band_edge - vbs[i].energy;
    if (delta < delta_max)
      tmp += vbs[i].degeneracy * std::pow(vbs[i].mass_DOS, 1.5)
        * std::exp(-delta / kT);
  }
  get_valence_band().effective_mass = std::pow(tmp, 2.0 / 3.0);
  
}



void
SemiconductorModel::calculate_equilibrium_properties(void)
{

  assert(bulk_model_ != NULL);

  // calculate conduction and valence band data
  bulk_model_->calculate_conduction_band_extremum();
  bulk_model_->calculate_valence_band_extremum();

  // get the band properties from bulk_model_
  extract_band_properties();

  // call the method of the parent class
  Parent::calculate_equilibrium_properties();

}


void
SemiconductorModel::do_print_info(void)
{
  string space("    ");

  cout << space << "unstrained semiconductor model" << endl;
  if (SimulationOptions::verbose() > 1)
  {
    cout << endl;
    set_lattice_temperature(SimulationOptions::T);
    bulk_model_->calculate_conduction_band_extremum();
    bulk_model_->calculate_valence_band_extremum();
    extract_band_properties();
    setup_band_edges();
    calculate_equilibrium_properties();

    double deg = std::pow(2.0, 2.0 / 3.0);

    const std::vector<DDsemiconductor::band_extremum>& cbs =
      bulk_model_->get_conduction_band_energy_mass();
    cout << space << " - conduction bands:\n";
    for (int i = 0 ; i < cbs.size(); i++)
    {
      cout << space << "   Ec = " << cbs[i].energy
        << ", m = " << cbs[i].mass_DOS
        << ", d = " << cbs[i].degeneracy << endl;
    }
    cout << space << "   Nc = " << get_conduction_band().effective_DOS << " cm^-3"
      << "  m_dos = " << get_conduction_band().effective_mass / deg << "\n";

    //_bulk_model->calculate_valence_band_extremum();
    const std::vector<DDsemiconductor::band_extremum>& vbs =
      bulk_model_->get_valence_band_energy_mass();
    cout << space << " - valence bands:\n";
    for (int i = 0 ; i < vbs.size(); i++)
    {
      cout << space << "   Ev = " << vbs[i].energy
        << ", m = " << vbs[i].mass_DOS
        << ", d = " << vbs[i].degeneracy << endl;
    }
    cout << space << "   Nv = " << get_valence_band().effective_DOS << " cm^-3"
      << "  m_dos = " << get_valence_band().effective_mass / deg << "\n";

    cout << space << " - Eg = " <<
      get_conduction_band().band_edge - get_valence_band().band_edge <<
      ", Ef0 = " << get_equilibrium_fermi_level()
      << ", ni = " << std::sqrt(get_intrinsic_density_squared());
    cout << endl;
  }
}
