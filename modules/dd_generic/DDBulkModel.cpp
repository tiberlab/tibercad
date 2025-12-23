/*  
 * This file is part of the tiberCAD module dd_generic.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file DDBulkModel.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "DDBulkModel.h"
#include "tibercad/physics/misc/ParticleDensity.h"
#include "RecombinationModelInterface.h"
#include "ThermoelectricPower.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/misc/Dopant.h"
#include "tibercad/physics/misc/Trap.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/base/InitFailedException.h"
#include "tibercad/embracing/Embracing.h"
#include "tibercad/io/Messages.h"
#include "tibercad/base/TypeDefs.h"
#include "tibercad/physics/misc/PolarizationModel.h"
#include "tibercad/physics/misc/PermittivityModel.h"
#include "tibercad/math/TensorOperators.h"

#include "libmesh/elem.h"


#include <cmath>


using namespace std;




DDBulkModel::~DDBulkModel(void)
{

}





DDBulkModel::DDBulkModel(const ModelOptions& options)
  : DriftDiffusionProperties(options),
    _is_inhomogeneous(false),
    _use_predictor(true),
    _eTEpowerGrad(0.0),
    _hTEpowerGrad(0.0),
    _eTEpower(0),
    _hTEpower(0),
    _polarization(0),
    _thermoelectric_power(NULL),
    _relax_polariz(1.0)
{
}



DDBulkModel*
DDBulkModel::create(const std::string& name, const Material* mat,
    const ModelOptions& options)
{
  return PhysicalModel::create<DDBulkModel>("ddbulk_" + name, mat, options);
}



void
DDBulkModel::read_database(void)
{
  const Database& db = get_database();
  db.set_section("");

}





void
DDBulkModel::parse_options(void)
{

  get_parameter("relax_polarization", _relax_polariz);

  _use_predictor = get_option("use_density_predictor", _use_predictor);





  // the temperature simulation
  string temp_simul = get_option("thermal_simulation", "");
  _is_inhomogeneous |= get_temperature_interface().set_simulation(temp_simul);


  // the strain simulation
  string strain_simul = get_option("strain_simulation", "");
  _is_inhomogeneous |= _strain_if.set_simulation(strain_simul);

}




void
DDBulkModel::prepare_submodels(void)
{

  DriftDiffusionProperties::prepare_submodels();

  //if (get_carrier_properties().size() == 0)
  //  throw InitFailedException("At least one carrier MUST be provided");

  // permittivity Default
  {
    PermittivityModel* pm = NULL;
    ModelOptions opts;
    opts.set_option("type", "constant");
    create_submodel(pm, "permittivity", opts);
  }


  // polarization models
  create_submodels(_pm, "polarization");


  {

    //
    // Recombinations
    //
    create_recombination_models();

    //
    // Thermoelectric power
    //
    // TODO bring into CarrierProperties

    auto it = get_options().submodels_begin("thermoelectric_power");
    auto end = get_options().submodels_end("thermoelectric_power");
    if (it != end)
    {
      _thermoelectric_power =
          ThermoelectricPower::create_model("default", get_material(), it->second);

      if (_thermoelectric_power == NULL)
        throw InitFailedException("Could not create thermoelectric power model");

      add_submodel("thermoelectricpower", _thermoelectric_power);
    }
  }

}




void
DDBulkModel::do_init(void)
{
  parse_options();

  DriftDiffusionProperties::do_init();


  // calculate the equilibrium
  set_lattice_temperature(SimulationOptions::T);

  calculate_equilibrium_properties();


  // permittivity
  SubmodelIterator it = submodels_begin("permittivity");
  PermittivityModel* pm =  dynamic_cast<PermittivityModel*>(it->second);
  set_relative_permittivity(pm->get_permittivity());


}





bool
DDBulkModel::has_solution(void) const
{
  return SimulationInterface::get_simulation(get_simulator_id())->is_solved();
}








void
DDBulkModel::do_reinit(const Elem* elem)
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

     _strain_if.get_crystal_strain(elem, elem->vertex_average(), get_strain());


      _polarization = 0;
      for (size_t n = 0; n < _pm.size(); n++)
      {
        _pm[n]->set_strain(get_strain());
        _pm[n]->calculate(get_element(), get_coordinates());
        _polarization += _pm[n]->get_polarization();
      }
      set_polarization(_polarization);

      //this->prepare_element_data();
  }

  // here we assume thermal equilibrium
  this->set_carrier_temperatures(get_lattice_temperature());

}








void
DDBulkModel::calculate_mobilities(void)
{
  PointData& pd = get_pd();
  pd.q_mobility.resize(this->n_known_carriers());
  pd.q_conductivity.resize(this->n_known_carriers());

  for (auto&& cp : get_carrier_properties())
  {
    pd.q_conductivity[cp.first] = cp.second->get_conductivity();
    pd.q_mobility[cp.first] = cp.second->get_mobility();
  }

}








void
DDBulkModel::calculate_equilibrium_properties(void)
{

  // call this method to properly set conduction and valence band DOS
  // and energy
  double kT = get_lattice_temperature();

  this->set_carrier_temperatures(kT);


  auto& carriers = get_carrier_properties();

  double Ec = 0.0, Ev = 0.0;
  int ncb = 0, nvb = 0;
  for (auto&& cp: carriers)
  {
    if (!cp.second->is_dopant())
    {
      cp.second->calculate(kT);
      if (cp.second->get_charge() < 0)
      {
        Ec += cp.second->get_band_edge();
        ncb++;
      }
      else if (cp.second->get_charge() > 0)
      {
        Ev += cp.second->get_band_edge();
        nvb++;
      }
    }
  }
  if (ncb > 0)
    Ec /= ncb;
  if (nvb > 0)
    Ev /= nvb;


  for (auto&& cp: carriers)
  {
    cp.second->use_quantum(false);
  }
  
/*
  double Nd = get_material()->get_total_donor_density();
  double Na = get_material()->get_total_acceptor_density();


  double ni2 = cb.get_effective_DOS() * vb.get_effective_DOS()
    * exp(-get_band_gap() / kT);
  double ni = sqrt(ni2);
  set_intrinsic_density(ni);

  double guess;
  // Hmm... Is there a better guess?
  if (Nd > Na)
  {
    guess = Ec - kT * log(cb.get_effective_DOS() / (Nd + ni));
  }
  else
  {
    guess = Ev + kT * log(vb.get_effective_DOS() / (Na + ni));
  }
*/

  // In some cases guess can be Inf or NaN. Then we set it to midband energy
  //if (std::isinf(guess) || std::isnan(guess))
  double guess = 0.5 * (Ec + Ev);

  /*
   * We use standard Newton. This should work always, as the density
   * is a strictly monotone function of the electric potential with
   * lim_{+-infty} = +-infty .
   */

  double x = guess;
  // 1e-4 V error seems to be good enough...
  double eps = 1e-4, dens_max = 1e10;
  double y;

  //set_carrier_temperatures(kT, kT);

  /* for testing
  ofstream of("charge.dat");
  double I = Ec + 0.2 - (Ev -0.2);
  double h = I / 10000;
  of << "# " << Ec << " " << Ev << "\n";
  for (unsigned int i = 0; i < 10000; i++)
  {
    double x = Ev - 0.2 + i * h;
    set_potentials(x);
    calculate_densities();
    calculate_traps();
    calculate_ionized_dopants();
    double f = get_charge_density();
    double df_fermi[2];
    get_charge_density_derivatives(df_fermi);
    double df = -(df_fermi[0] + df_fermi[1]);
    of << (x-Ec) << " " << f << " " << df <<  "\n";
  }
  of.close();
  */

  // is set to true when calculation is successful
  bool success = false;

  // the minimum x (below the zero, f > 0)
  double xmin = Ev - 0.5;

  // the maximum x (above the zero, f < 0)
  double xmax = Ec + 0.5;

  
  //cerr << "***" << get_owner()->get_name() << endl;;
  for (unsigned int i = 0; i < 10000; ++i)
  {
    set_potentials(x);
    calculate_densities();
    //calculate_traps();
    calculate_ionized_dopants();
    double f = get_charge_density();
    double df = 0;
    for (auto&& c : this->get_carrier_properties())
      df -= get_charge_density_derivative(c.first);

    //cout << "x = " << x << " f = " << f << " df = " << df << endl;

    if (f > 0) xmin = x;
    else if (f < 0) xmax = x;

    double residual_dens = fabs(f);

    double dx = 0.0;
    if (residual_dens > ParticleDensity::MINDENSITY)
    {
      // At low temperatures everything is very sensitive on dx, so we don't
      // allow it to be bigger than k*T. At high temperatures this should not
      // have any impact
      dx = - f / df;

      y = x + dx;
      // we limit Ef to (xmin, xmax)
      if (y > xmax)
        dx = 0.5 * (xmax - x);
      else if (y < xmin)
        dx = 0.5 * (xmin - x);

    }

    y = x + dx;

    double error = fabs(dx);
    //cout << "x = " << y << " error = " << dx << " res. dens. = "
    //  << residual_dens << " Ec = " << Ec << " Ev = " << Ev << endl;

    x = y;

    if ((error < eps) && (residual_dens < dens_max))
    {
      success = true;
      break;
    }
  }

  if (!success)
  {
    ostringstream os;
    os << "Could not find equilibrium properties for material "
        << get_material()->get_name() << " around point "
        << get_coordinates();
    //throw SolveFailedException(os.str());
  }

  //set_intrinsic_density(sqrt(get_electron_density()) * sqrt(get_hole_density()));

  //set_equilibrium_n(get_electron_density());
  //set_equilibrium_p(get_hole_density());

  // for a dielectric we don't need much...
  set_equilibrium_fermi_level(y);


  for (auto&& cp: carriers)
  {
     bool quantum = cp.second->has_quantum();
     cp.second->use_quantum(quantum);

     set_bulk_equilibrium_density(cp.first, get_q_density(cp.first));

  }

}





void
DDBulkModel::set_equilibrium_properties(double Ef)
{

  set_equilibrium_fermi_level(Ef);
  set_potentials(Ef);
  calculate_densities();

  //set_intrinsic_density(sqrt(get_electron_density()) * sqrt(get_hole_density()));
  //set_equilibrium_n(get_electron_density());
  //set_equilibrium_p(get_hole_density());

}








void
DDBulkModel::compute_thermoelectric_powers(void)
{
  _eTEpower = 0;
  _hTEpower = 0;

  if (_thermoelectric_power != NULL)
  {
    /*
    _thermoelectric_power->set_potentials(
        get_electron_electro_chemical_potential(),
        get_hole_electro_chemical_potential(),
        get_electric_potential());

    double cb = get_conduction_band_edge();

    double vb = get_valence_band_edge();

    _thermoelectric_power->set_band_edges(cb, vb);

    _thermoelectric_power->set_temperature(get_lattice_temperature());

    _thermoelectric_power->calculate();

    _eTEpower = _thermoelectric_power->get_electrons_thermoelectric_power();

    _hTEpower = _thermoelectric_power->get_holes_thermoelectric_power();
    */
  }
}

void
DDBulkModel::compute_thermoelectric_power_gradient(void)
{
  _eTEpowerGrad = 0;
  _hTEpowerGrad = 0;
   if (_thermoelectric_power != NULL)
   {
     /*
     _thermoelectric_power->set_potential_gradients(get_grad_fermi_e(),
         get_grad_fermi_h(),get_electric_field());

     _thermoelectric_power->set_temperature(get_lattice_temperature());

     _thermoelectric_power->calculate_derivatives();

     _eTEpowerGrad =  _thermoelectric_power->get_electron_thermoelectric_power_gradient();

     _hTEpowerGrad =  _thermoelectric_power->get_hole_thermoelectric_power_gradient();
      */
    }
}

std::vector<double>&
DDBulkModel::get_temperature_at_nodes()
{
  return _nodal_lattice_vt;
}



void
DDBulkModel::do_print_info(void)
{
  Messages m;

  if (_strain_if.has_simulation())
    m.info("using strain simulation: " +
      _strain_if.get_simulation()->get_name());
  else if (trace(get_strain()) != 0.0)
    m.info("using strain from input file");

  if (get_temperature_interface().has_simulation())
    m.info("using lattice temperature from: " +
        get_temperature_interface().get_simulation()->get_name());


  set_lattice_temperature(SimulationOptions::T);
  calculate_equilibrium_properties();

  ostringstream os;

  for (auto&& cp: get_carrier_properties())
  {
    m.info("Carrier:");
    m.indent();
    cp.second->print_info();

    os << "Band edge = " << cp.second->get_band_edge()
       << ", N0 = " << cp.second->get_effective_DOS() << " cm^-3\n"
       << "m_dos = " << cp.second->get_effective_mass()
       << ", v_th = " << cp.second->get_thermal_velocity(Constants::k_B * SimulationOptions::T) << " cm/s\n\n";

    Messages::info(os.str());
    m.unindent();
    os.str("");
  }

  Messages::endl;
  Messages::info(os.str());

  /*
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
  */
/*
  os << "Ec = " << get_conduction_band().get_band_edge()
      << ", Nc = " << get_conduction_band().get_effective_DOS() << " cm^-3\n"
      << "m_dos = " << get_conduction_band().get_effective_mass() / deg
      << ", v_th = " << get_conduction_band().get_thermal_velocity(
          Constants::k_B * SimulationOptions::T)
      << " cm/s\n";
  Messages::info(os.str());
*/
  /*
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
  */
/*
  os << "Ev = " << get_valence_band().get_band_edge()
      << ", Nv = " << get_valence_band().get_effective_DOS() << " cm^-3\n"
      << "m_dos = " << get_valence_band().get_effective_mass() / deg
      << ", v_th = " << get_valence_band().get_thermal_velocity(
          Constants::k_B *SimulationOptions::T)
      << " cm/s\n";
  Messages::info(os.str());
*/
  /*
  m.unindent();
  os.str("");

  os << "Eg = "
      << get_conduction_band().get_band_edge() - get_valence_band().get_band_edge()
      << ", Ef0 = " << get_equilibrium_fermi_level()
      << ", ni = " << std::sqrt(get_intrinsic_density_squared())
      << Messages::endl;
  Messages::info(os.str());
  */

}

