// $Id$

#include "SimpleSemiconductorModel.h"
#include "Material.h"
#include "Messages.h"


TIBER_MODULE(SimpleSemiconductorModel, ddbulk, simple)



using namespace std;
using namespace DriftDiffusionDefs;


SimpleSemiconductorModel::SimpleSemiconductorModel(const ModelOptions& options)
  : DriftDiffusionProperties(options),
    is_prepared_(false)
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
  get_conduction_band()._band_edge = get_option("Ec", 2.2288);
  //get_conduction_band().band_edges.resize(1);
  //get_conduction_band().band_edges[0] = get_conduction_band().band_edge;
  get_valence_band()._band_edge = get_option("Ev", 1.1047);
  //get_valence_band().band_edges.resize(1);
  //get_valence_band().band_edges[0] = get_valence_band().band_edge;

  double deg = std::pow(2.0, 2.0 / 3.0);
  get_conduction_band()._effective_mass = deg * get_option("m_dos_e", 1.082);
  get_valence_band()._effective_mass = deg * get_option("m_dos_h", 1.1432);
  get_conduction_band()._degeneracy = 2;
  get_valence_band()._degeneracy = 2;

  if (has_option("Nc"))
  {
    double Nc = get_option("Nc", 1e20);
    double kT = Constants::k_B * SimulationOptions::T;
    get_conduction_band()._effective_mass =
      std::pow(Nc / get_DOS_factor(), 2.0/3.0) / kT;
  }

  if (has_option("Nv"))
  {
    double Nv = get_option("Nv", 1e20);
    double kT = Constants::k_B * SimulationOptions::T;
    get_valence_band()._effective_mass =
      std::pow(Nv / get_DOS_factor(), 2.0/3.0) / kT;
  }

  DriftDiffusionProperties::do_init();
}

void
SimpleSemiconductorModel::do_print_info(void)
{

  set_lattice_temperature(SimulationOptions::T);
  setup_band_edges();
  calculate_equilibrium_properties();

  double deg = std::pow(2.0, 2.0 / 3.0);
  Messages::info("simple semiconductor model " 
    "(with constant band parameters)");
  DriftDiffusionProperties::do_print_info();

  if (SimulationOptions::verbose() > 1)
  {
    ostringstream os;
    os << "Ec = " << get_conduction_band().get_band_edge() <<
      ", m_DOS = " << get_conduction_band().get_effective_mass() / deg <<
      ", Nc = " << get_conduction_band().get_effective_DOS() << Messages::endl;
    os << "Ev = " << get_valence_band().get_band_edge() <<
      ", m_DOS = " << get_valence_band().get_effective_mass() / deg <<
      ", Nv = " << get_valence_band().get_effective_DOS() << Messages::endl;
    os << "Eg = " <<
      get_conduction_band().get_band_edge() - get_valence_band().get_band_edge() <<
      ", Ef0 = " << get_equilibrium_fermi_level()
      << ", ni = " << std::sqrt(get_intrinsic_density_squared());

    Messages::info(os.str());
    Messages::newline();
  }
}
