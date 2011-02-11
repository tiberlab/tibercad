// $Id$


#include "SemiconductorModel.h"
#include "DDsemiconductor.h"

#include "Material.h"
#include "Messages.h"


#include <sstream>


TIBER_MODULE(SemiconductorModel, ddbulk, default)


using namespace DriftDiffusionDefs;
using namespace std;


SemiconductorModel::~SemiconductorModel(void)
{
  reset();
}

SemiconductorModel::SemiconductorModel(const ModelOptions& options)
  : DriftDiffusionProperties(options),
    _is_prepared(false),
    _bulk_model(NULL),
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
SemiconductorModel::create_submodels(void)
{

  DriftDiffusionProperties::create_submodels();

  // this is dirty as for now the band parameters are not given
  // by a physical_model
  ModelOptions opt(get_options());
  opt.delete_option("relax_polarization");
  opt.delete_option("strain_simulation");
  opt.delete_option("statistics");
  opt.delete_all_submodels();

  if (_bulk_model == NULL)
  {
    _bulk_model = DDsemiconductor::create(get_material()->get_structure(), opt);

    if (_bulk_model == NULL)
      throw InitFailedException("Unknown structure for DDsemiconductor");

    add_submodel("bandstructure", _bulk_model);

  }
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
      elem_data.mc = get_conduction_band().effective_mass;
      elem_data.mv = get_valence_band().effective_mass;
      elem_data.Ef0 = get_equilibrium_fermi_level();
      //elem_data.ni = get_intrinsic_density();

    }
    else
    {
      const ElementData& elem_data = it->second;

      get_conduction_band().band_edge = elem_data.Ec;
      get_conduction_band().effective_mass = elem_data.mc;
      get_valence_band().band_edge = elem_data.Ev;
      get_valence_band().effective_mass = elem_data.mv;

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


//
// Very crude implementation at the moment
//
void
SemiconductorModel::extract_band_properties(void)
{
  // treat conduction band
  const std::vector<DDsemiconductor::band_extremum>& cbs =
    _bulk_model->get_conduction_band_energy_mass();

  get_conduction_band().band_edges.resize(1);

  // get minimum
  int id = 0;
  get_conduction_band().band_edges[0] = cbs[0].energy;

  for (unsigned int i = 1; i < cbs.size(); i++)
  {
    if (cbs[i].energy < cbs[id].energy)
      id = i;
  }
  get_conduction_band().band_edge = cbs[id].energy;
  get_conduction_band().effective_mass = cbs[id].mass_DOS
    * std::pow(cbs[id].degeneracy, 2.0 / 3.0);
  get_conduction_band().degeneracy = cbs[id].degeneracy;

  // treat valence band
  const std::vector<DDsemiconductor::band_extremum>& vbs =
    _bulk_model->get_valence_band_energy_mass();

  get_valence_band().band_edges.resize(vbs.size());

  // get maximum
  id = 0;
  get_valence_band().band_edges[0] = vbs[0].energy;

  //double kT = SimulationOptions::T * Constants::k_B;
  double kT = get_lattice_temperature();
  double delta_max = 4.0 * kT;
  for (unsigned int i = 1; i < vbs.size(); i++)
  {
    get_valence_band().band_edges[i] = vbs[i].energy;
    if (vbs[i].energy > vbs[id].energy)
      id = i;
  }
  get_valence_band().band_edge = vbs[id].energy;
  get_valence_band().degeneracy = vbs[id].degeneracy;

  double tmp = 0;
  // include other bands
  for (unsigned int i = 0; i < vbs.size(); i++)
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

  assert(_bulk_model != NULL);

  // it wants temperature in K
  _bulk_model->set_temperature(get_lattice_temperature() / Constants::k_B);
  _bulk_model->set_strain(get_strain());

  // calculate conduction and valence band data
  _bulk_model->calculate_conduction_band_extremum();
  _bulk_model->calculate_valence_band_extremum();

  // get the band properties from _bulk_model
  extract_band_properties();

  // call the method of the parent class
  Parent::calculate_equilibrium_properties();

}



void
SemiconductorModel::do_print_info(void)
{

  Messages::info("default semiconductor model (using k.p)");
  Parent::do_print_info();

  if (SimulationOptions::verbose() > 1)
  {
    Messages::newline();
    set_lattice_temperature(SimulationOptions::T);

    calculate_equilibrium_properties();
    setup_band_edges();

    double deg = std::pow(2.0, 2.0 / 3.0);

    const std::vector<DDsemiconductor::band_extremum>& cbs =
      _bulk_model->get_conduction_band_energy_mass();

    const std::vector<DDsemiconductor::band_extremum>& vbs =
      _bulk_model->get_valence_band_energy_mass();

    ostringstream os;

    os << " - Eg = " <<
      get_conduction_band().band_edge - get_valence_band().band_edge <<
      ", Ef0 = " << get_equilibrium_fermi_level()
      << ", ni = " << std::sqrt(get_intrinsic_density_squared())
      << Messages::endl;

    os << " - conduction bands:" << Messages::endl;
    for (int i = 0 ; i < cbs.size(); i++)
    {
      os << "   Ec = " << cbs[i].energy
        << ", m = " << cbs[i].mass_DOS
        << ", d = " << cbs[i].degeneracy << Messages::endl;
    }
    os << "   Nc = " << get_conduction_band().effective_DOS << " cm^-3"
      << "  m_dos = " << get_conduction_band().effective_mass / deg
      << Messages::endl;

    os << " - valence bands:" << Messages::endl;
    for (int i = 0 ; i < vbs.size(); i++)
    {
      os << "   Ev = " << vbs[i].energy
        << ", m = " << vbs[i].mass_DOS
        << ", d = " << vbs[i].degeneracy << Messages::endl;
    }
    os << "   Nv = " << get_valence_band().effective_DOS << " cm^-3"
      << "  m_dos = " << get_valence_band().effective_mass / deg;

    Messages::info(os.str());
    Messages::newline();
  }
}
