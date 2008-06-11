// $Id$


#include "SemiconductorModel.h"
#include "ZbDDsemiconductor.h"
#include "WzDDsemiconductor.h"

#include "Alloy.h"

//#include "Database.h"

//#include "getpot.h"

#include <iostream>


TIBER_MODULE(SemiconductorModel, default)


using namespace DriftDiffusionDefs;
using namespace std;


SemiconductorModel::~SemiconductorModel(void)
{
  reset();
  PhysicalModelInterface::destroy(bulk_model_);
}

SemiconductorModel::SemiconductorModel(void)
  : bulk_model_(NULL),
    is_prepared_(false),
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

  const ModelOptions& opt = get_options();

  _recompute_band_parameters = get_parameter("recompute_band_parameters",
      _recompute_band_parameters);

  PhysicalModelInterface::destroy(bulk_model_);

  bulk_model_ = DDsemiconductor::create(get_material()->get_structure(), opt);

  if (bulk_model_ == NULL)
    throw InitFailedException("Unknown structure for DDsemiconductor");

  bulk_model_->set_material(get_material());
  bulk_model_->set_simulator_id(get_simulator_id());

  bulk_model_->init();
  
}



void
SemiconductorModel::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  Parent::calculate_VCA(comp_A, comp_B, xa);

  const SemiconductorModel* scA =
    dynamic_cast<const SemiconductorModel*>(comp_A);
  const SemiconductorModel* scB =
    dynamic_cast<const SemiconductorModel*>(comp_B);

  bulk_model_->build_alloy(scA->bulk_model_, scB->bulk_model_, xa);

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
      elem_data.mc = get_conduction_band().effective_mass;
      elem_data.mv = get_valence_band().effective_mass;
      elem_data.Ef0 = get_equilibrium_fermi_level();
      elem_data.ni = get_intrinsic_density();

    }
    else
    {
      const ElementData& elem_data = it->second;

      //set_polarization(elem_data.polarization);

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
  {
    if (!is_prepared_ || _recompute_band_parameters)
    {
      calculate_equilibrium_properties();
      is_prepared_ = true;
    }
  }
}


//
// Very crude implementation at the moment
// 
void
SemiconductorModel::extract_band_properties(void)
{
  // treat conduction band
  const std::vector<DDsemiconductor::band_extremum>& cbs =
    bulk_model_->get_conduction_band_energy_mass();

  get_conduction_band().band_edges.resize(1);

  // get minimum
  int id = 0;
  get_conduction_band().band_edges[0] = cbs[0].energy;

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
  
  get_valence_band().band_edges.resize(vbs.size());

  // get maximum
  id = 0;
  get_valence_band().band_edges[0] = vbs[0].energy;
  
  //double kT = SimulationOptions::T * Constants::k_B;
  double kT = get_lattice_temperature();
  double delta_max = 4.0 * kT;
  for (int i = 1; i < vbs.size(); i++)
  {
    get_valence_band().band_edges[i] = vbs[i].energy;
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

    //cerr << i << ", Ev = " << vbs[i].energy << ", m = " << vbs[i].mass_DOS << ", d = " << vbs[i].degeneracy << endl;
  }
  get_valence_band().effective_mass = std::pow(tmp, 2.0 / 3.0);
  //cerr << "DOS mass = " << get_valence_band().effective_mass  << "\n***\n";

  
}



void
SemiconductorModel::calculate_equilibrium_properties(void)
{

  assert(bulk_model_ != NULL);

  // it wants temperature in K
  bulk_model_->set_temperature(get_lattice_temperature() / Constants::k_B);
  bulk_model_->set_strain(get_strain());

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

  cout << space << "default semiconductor model (using k.p)" << endl;
  Parent::do_print_info();

  if (SimulationOptions::verbose() > 1)
  {
    cout << endl;
    set_lattice_temperature(SimulationOptions::T);

    calculate_equilibrium_properties();
    setup_band_edges();

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
