// $Id$

#include "SimpleSemiconductorModel.h"
#include "Material.h"


TIBER_MODULE(SimpleSemiconductorModel, simple)



using namespace std;
using namespace DriftDiffusionDefs;


SimpleSemiconductorModel::SimpleSemiconductorModel(void)
  : is_prepared_(false)
{
}


void
SimpleSemiconductorModel::prepare_element_data(void)
{
  if (!is_prepared_)
  {
    calculate_equilibrium_properties();

    is_prepared_ =  true;
  }
}

void
SimpleSemiconductorModel::do_init(void)
{
  DriftDiffusionProperties::do_init();

  // for the moment we read them from the materials section
  ModelOptions& opt = get_material()->get_options();

  get_conduction_band().band_edge = opt.get_option("Ec", 2.2288);
  get_valence_band().band_edge = opt.get_option("Ev", 1.1047);
  double deg = std::pow(2.0, 2.0 / 3.0);
  get_conduction_band().effective_mass = deg * opt.get_option("meff_n", 1.082);
  get_valence_band().effective_mass = deg * opt.get_option("meff_p", 1.1432);
  
}

void
SimpleSemiconductorModel::do_print_info(void)
{
  string space("    ");
  
  set_lattice_temperature(SimulationOptions::T);
  setup_band_edges();
  calculate_equilibrium_properties();

  double deg = std::pow(2.0, 2.0 / 3.0);
  cout << space << "simple semiconductor model" << endl;
  cout << space << "Ec = " << get_conduction_band().band_edge <<
    ", m_DOS = " << get_conduction_band().effective_mass / deg <<
    ", Nc = " << get_conduction_band().effective_DOS << endl;
  cout << space << "Ev = " << get_valence_band().band_edge <<
    ", m_DOS = " << get_valence_band().effective_mass / deg <<
    ", Nv = " << get_valence_band().effective_DOS << endl;
  cout << space << "Eg = " <<
    get_conduction_band().band_edge - get_valence_band().band_edge <<
    ", Ef0 = " << get_equilibrium_fermi_level()
      << ", ni = " << std::sqrt(get_intrinsic_density_squared());
  cout << endl;
}
