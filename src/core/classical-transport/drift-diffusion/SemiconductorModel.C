// $Id$

#include "Alloy.h"

#include "SemiconductorModel.h"
#include "ZbDDsemiconductor.h"
#include "WzDDsemiconductor.h"

#include "Database.h"

#include "getpot.h"

#include <iostream>

using namespace DriftDiffusionDefs;
using namespace std;


SemiconductorModel::~SemiconductorModel(void)
{
  PhysicalModelInterface::destroy(_bulk_model);
}

SemiconductorModel::SemiconductorModel(void)
  : _bulk_model(NULL),
    _is_prepared(false)
{
}


void
SemiconductorModel::do_init(void)
{
  Parent::do_init();

  const ModelOptions& opt =  get_options ();

  PhysicalModelInterface::destroy(_bulk_model);

  _bulk_model = DDsemiconductor::create(get_material()->get_structure(), opt);

  if (_bulk_model == NULL)
    throw InitFailedException("Unknown structure for DDsemiconductor");

  _bulk_model->set_material(get_material());

  _bulk_model->init();
  
}

/*
void
SemiconductorModel::read_bowing_parameters(void)
{
}
*/

void
SemiconductorModel::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  permittivity = data("permittivity", 12.93);
 

  if (get_material()->get_structure() == "zb")
  {
    permittivity = data("permittivity", 12.93);
  }
  else
  {
    permittivity = data("permittivity", 9.5);
  }
}


void
SemiconductorModel::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  DriftDiffusionProperties::calculate_VCA(comp_A, comp_B, xa);

  double (*alloy)(double, double, double, double) =
    Alloy::calculate_VCA_parameter;

  const SemiconductorModel* scA =
    dynamic_cast<const SemiconductorModel*>(comp_A);
  const SemiconductorModel* scB =
    dynamic_cast<const SemiconductorModel*>(comp_B);

  _bulk_model->build_alloy(scA->_bulk_model, scB->_bulk_model, xa);

  permittivity = alloy(scA->permittivity, scB->permittivity,
        xa, 0.0);

}



void
SemiconductorModel::prepare_element_data(void)
{
  if (!_is_prepared)
  {
    try
    {
      calculate_equilibrium_properties();
    }
    catch (...)
    {
    }
  
    _is_prepared = true;
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
    _bulk_model->get_conduction_band_energy_mass();
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
    _bulk_model->get_valence_band_energy_mass();
  // get maximum
  id = 0;
  double kT = SimulationOptions::T * Constants::k_B;
  double delta_max = 4.0 * kT;
  for (int i = 1; i < vbs.size(); i++)
  {
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
SemiconductorModel::print_info(void) const
{
  _bulk_model->calculate_conduction_band_extremum();
  const std::vector<DDsemiconductor::band_extremum>& cbs =
    _bulk_model->get_conduction_band_energy_mass();
  cout << " - conduction bands:\n";
  for (int i = 0 ; i < cbs.size(); i++)
  {
    cout << "   Ec = " << cbs[i].energy
      << ", m = " << cbs[i].mass_DOS
      << ", d = " << cbs[i].degeneracy << endl;
  }
  cout << "   Nc = " << get_conduction_band().effective_DOS << " cm^-3"
    << "  m_dos = " << get_conduction_band().effective_mass << "\n";

  _bulk_model->calculate_valence_band_extremum();
  const std::vector<DDsemiconductor::band_extremum>& vbs =
    _bulk_model->get_valence_band_energy_mass();
  cout << " - valence bands:\n";
  for (int i = 0 ; i < vbs.size(); i++)
  {
    cout << "   Ev = " << vbs[i].energy
      << ", m = " << vbs[i].mass_DOS
      << ", d = " << vbs[i].degeneracy << endl;
  }
  cout << "   Nv = " << get_valence_band().effective_DOS << " cm^-3"
    << "  m_dos = " << get_valence_band().effective_mass << "\n";

  cout << " - Ef0 = " << get_equilibrium_fermi_level()
    << ", ni^2 = " << get_intrinsic_density_squared();
  cout << endl;
}

void
SemiconductorModel::calculate_equilibrium_properties(void)
{

  assert(_bulk_model != NULL);

  // calculate conduction and valence band data
  _bulk_model->calculate_conduction_band_extremum();
  _bulk_model->calculate_valence_band_extremum();

  // get the band properties from _bulk_model
  extract_band_properties();

  // call the method of the parent class
  Parent::calculate_equilibrium_properties();

}

