
#include "MEBulkModel.h"
#include "ParticleDensity.h"
//#include "RecombinationModelInterface.h"
//#include "MobilityModelInterface.h"
//#include "ThermoelectricPower.h"
#include "SimulationInterface.h"
#include "Material.h"
#include "Database.h"
//#include "Dopant.h"
//#include "Trap.h"
#include "Constants.h"
#include "InitFailedException.h"
//#include "RotatedCrystal.h"
#include "Embracing.h"
#include "Messages.h"
#include "TypeDefs.h"
//#include "PolarizationModel.h"
//#include "PermittivityModel.h"

#include "elem.h"


#include <cmath>


using namespace std;




MEBulkModel::~MEBulkModel(void)
{

}



MEBulkModel::MEBulkModel(const ModelOptions& options)
  : MasterEquationsProperties(options)

{
}



MEBulkModel*
MEBulkModel::create(const std::string& name, const Material* mat,
    const ModelOptions& options)
{
  return PhysicalModel::create<MEBulkModel>("me_bulk_" + name, mat, options);
}



void
MEBulkModel::read_database(void)
{
  const Database& db = get_database();

  //db.set_section("");
  //is_dielectric() = db.get("dielectric", is_dielectric());

  //bool diel_as_sc = get_option("dielectric_as_semiconductor", false);
  //is_dielectric() &= !diel_as_sc;

}





void
MEBulkModel::parse_options(void)
{

//  is_dielectric() = get_option("dielectric", is_dielectric());

//  bool diel_as_sc = get_option("dielectric_as_semiconductor", false);
//  is_dielectric() &= !diel_as_sc;


}




void
MEBulkModel::do_init(void)
{

  parse_options();

  MasterEquationsProperties::do_init();


  // calculate the equilibrium
  set_lattice_temperature(SimulationOptions::T);
  //set_equilibrium_properties();


}



bool
MEBulkModel::has_solution(void) const
{
  return SimulationInterface::get_simulation(get_simulator_id())->is_solved();
}



void
MEBulkModel::do_reinit(const Elem* elem)
{

  if (get_element() != elem)
  {
    set_element(elem);
    set_coordinates(elem->vertex_average());

    // get the nodal temperatures
    get_temperature_interface().get_temperature(elem, _nodal_lattice_vt);

    // get the mean temperature on the element
    set_lattice_temperature(
        get_temperature_interface().get_temperature(elem, elem->vertex_average()));

  }

  // here we assume thermal equilibrium
  electron_vt = hole_vt = get_lattice_temperature();

}





void
MEBulkModel::calculate_equilibrium_properties(void)
{
  // remember the coupling va in set equilibrium properties
  //int coupling_bkp = get_coupling_type();
  //set_coupling_type(MasterEquationsDefs::BOTH);



  // call this method to properly set conduction and valence band DOS
  // and energy
  double kT = get_lattice_temperature();
  set_carrier_temperatures(kT, kT);


  MEBandProperties& cb = get_conduction_band();
  cb.calculate(kT);
  double Ec = cb.get_band_edge();


  MEBandProperties& vb = get_valence_band();
  vb.calculate(kT);
  double Ev = vb.get_band_edge();


  double ni2 = cb.get_effective_DOS() * vb.get_effective_DOS()
    * exp(-get_band_gap() / kT);
  double ni = sqrt(ni2);
  set_intrinsic_density(ni);

  set_equilibrium_fermi_level(0.5 * (Ec + Ev));

  calculate_densities();

  set_intrinsic_density(sqrt(get_electron_density()) * sqrt(get_hole_density()));
  set_equilibrium_n(get_electron_density());
  set_equilibrium_p(get_hole_density());

}



void
MEBulkModel::set_equilibrium_properties(double Ef)
{
  // remember the coupling
  int coupling_bkp = get_coupling_type();
  set_coupling_type(MasterEquationsDefs::BOTH);

  set_equilibrium_fermi_level(Ef);
  //set_potentials(Ef);
  calculate_densities();

  set_intrinsic_density(sqrt(get_electron_density()) * sqrt(get_hole_density()));
  set_equilibrium_n(get_electron_density());
  set_equilibrium_p(get_hole_density());

  // restore original coupling
  set_coupling_type(coupling_bkp);
}



std::vector<double>&
MEBulkModel::get_temperature_at_nodes()
{
  return _nodal_lattice_vt;
}



void
MEBulkModel::do_print_info(void)
{
  Messages m;

  if (get_temperature_interface().has_simulation())
    m.info("using lattice temperature from: " +
        get_temperature_interface().get_simulation()->get_name());


  set_lattice_temperature(SimulationOptions::T);
  //calculate_equilibrium_properties();

  //double deg = std::pow(2.0, 2.0 / 3.0);

  m.info("Conduction band:");
  m.indent();
  get_conduction_band().print_info();

  ostringstream os;
  os << "Ec = " << get_conduction_band().get_band_edge()
      << ", Nc = " << get_conduction_band().get_effective_DOS() << " cm^-3\n"
      << "m_dos = " << get_conduction_band().get_effective_mass()
      << ", v_th = " << get_conduction_band().get_thermal_velocity(
          Constants::k_B * SimulationOptions::T)
      << " cm/s\n";
  Messages::info(os.str());

  m.unindent();

  os.str("");

  m.info("Valence band:");
  m.indent();
  get_valence_band().print_info();

  os << "Ev = " << get_valence_band().get_band_edge()
      << ", Nv = " << get_valence_band().get_effective_DOS() << " cm^-3\n"
      << "m_dos = " << get_valence_band().get_effective_mass()
      << ", v_th = " << get_valence_band().get_thermal_velocity(
          Constants::k_B *SimulationOptions::T)
      << " cm/s\n";
  Messages::info(os.str());


  m.unindent();
  os.str("");

  os << "Eg = "
      << get_conduction_band().get_band_edge() - get_valence_band().get_band_edge()
      << ", Ef0 = " << get_equilibrium_fermi_level()
      << ", ni = " << std::sqrt(get_intrinsic_density_squared())
      << Messages::endl;
  Messages::info(os.str());
}

