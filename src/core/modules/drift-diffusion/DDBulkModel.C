// $Id$

#include "DDBulkModel.h"
#include "ParticleDensity.h"
#include "RecombinationModelInterface.h"
#include "MobilityModelInterface.h"
#include "ThermoelectricPower.h"
#include "SimulationInterface.h"
#include "Material.h"
#include "Database.h"
#include "Dopant.h"
#include "Trap.h"
#include "Constants.h"
#include "InitFailedException.h"
#include "RotatedCrystal.h"
#include "Embracing.h"
#include "Messages.h"
#include "TypeDefs.h"
#include "PolarizationModel.h"
#include "PermittivityModel.h"

#include "elem.h"


#include <cmath>


using namespace std;




DDBulkModel::~DDBulkModel(void)
{

}





DDBulkModel::DDBulkModel(const ModelOptions& options)
  : DriftDiffusionProperties(options),
    _is_inhomogeneous(false),
    _use_predictor(true),
    _electron_mobility(NULL),
    _hole_mobility(NULL),
    _eTEpowerGrad(0.0),
    _hTEpowerGrad(0.0),
    _eTEpower(0),
    _hTEpower(0),
    _polarization(0),
    _background_conductivity(0.0),
    _thermoelectric_power(NULL),
    _is_dielectric(false),
    _relax_polariz(1.0)
{
}



DDBulkModel*
DDBulkModel::create(const std::string& name, const Material* mat,
    const ModelOptions& options)
{
  return PhysicalModelInterface::create<DDBulkModel>("ddbulk_" + name, mat, options);
}



void
DDBulkModel::read_database(void)
{
  const Database& db = get_database();
  db.set_section("");
  _is_dielectric = db.get("dielectric", _is_dielectric);

  bool diel_as_sc = get_option("dielectric_as_semiconductor", false);
  _is_dielectric &= !diel_as_sc;

}





void
DDBulkModel::parse_options(void)
{

  get_parameter("relax_polarization", _relax_polariz);

  _use_predictor = get_option("use_density_predictor", _use_predictor);

  _is_dielectric = get_option("dielectric", _is_dielectric);

  bool diel_as_sc = get_option("dielectric_as_semiconductor", false);
  _is_dielectric &= !diel_as_sc;



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

  // permittivity Default
  {
    PermittivityModel* pm = NULL;
    ModelOptions opts;
    opts.set_option("type", "constant");
    create_submodel(pm, "permittivity", opts);
  }


  if (&get_conduction_band() == NULL)
  {
    BandProperties* cb;
    ModelOptions opts;
    opts.set_option("particle", "el");
    opts.set_option("type", "kp");
    create_submodel(cb, "band_properties", opts);
    set_conduction_band(cb);
  }


  if (&get_valence_band() == NULL)
  {
    BandProperties* vb;
    ModelOptions opts;
    opts.set_option("particle", "hl");
    opts.set_option("type", "kp");
    create_submodel(vb, "band_properties", opts);
    set_valence_band(vb);
  }


/*
  // particle densities
  {
    ModelOptions opts;
    opts.set_option("statistics", "fermidirac");

    bool e_done = false;
    bool h_done = false;

    ModelOptions::submodel_iterator
      it(get_options().submodels_begin("particle_density"));
    ModelOptions::submodel_iterator
      end(get_options().submodels_end("particle_density"));



    while (it != end)
    {
      ModelOptions& o = it->second;
      ++it;

      const string& particle = o.get_option("particle", "");
      if (!o.find_option("statistics"))
        o.set_option("statistics", "fermidirac");

      if (particle == "electron")
      {
        if (e_done)
          throw InitFailedException("Only one particle_density model per "
              "particle and region is allowed");

        e_done = true;
      }
      else if (particle == "hole")
      {
        if (h_done)
          throw InitFailedException("Only one particle_density model per "
              "particle and region is allowed");

        h_done = true;
      }
      else
      {
        if (e_done || h_done)
          throw InitFailedException("Only one particle_density model per "
              "particle and region is allowed");

        // this case is valid for both

        o.set_option("particle", "electron");

        opts = o;
        opts.set_option("particle", "hole");
        get_options().add_submodel("particle_density", opts);

        e_done = h_done = true;
      }
    }

    if (!e_done)
    {
      opts.set_option("particle", "electron");
      get_options().add_submodel("particle_density", opts);
    }

    if (!h_done)
    {
      opts.set_option("particle", "hole");
      get_options().add_submodel("particle_density", opts);
    }

    vector<PhysicalModelInterface*> pd;
    create_submodels(pd, "particle_density");
  }
*/

  // polarization models
  create_submodels(_pm, "polarization");


  if (!is_dielectric())
  {

    //
    // mobilities
    //
    ModelOptions::submodel_iterator
      mobit(get_options().submodels_begin("mobility"));

    ModelOptions::submodel_iterator
      it(get_options().submodels_begin("electron_mobility"));
    ModelOptions::submodel_iterator
      end(get_options().submodels_end("electron_mobility"));

    // create electron mobility model
    if (it != end)
    {
      (it->second).set_option("particle", "electron");
      _electron_mobility = create_mobility_model(it->second);
    }
    else if (mobit != get_options().submodels_end("mobility"))
    {
      (mobit->second).set_option("particle", "electron");
      (mobit->second).set_option("name", string("electron_mobility"));
      _electron_mobility = create_mobility_model(mobit->second);
    }
    else
    {
      ModelOptions opts;
      opts.set_option("particle", "electron");
      opts.set_option("name", string("electron_mobility"));
      _electron_mobility = create_mobility_model(opts);
    }



    // create hole mobility model
    it = get_options().submodels_begin("hole_mobility");
    end = get_options().submodels_end("hole_mobility");

    if (it != end)
    {
      (it->second).set_option("particle", "hole");
      _hole_mobility = create_mobility_model(it->second);
    }
    else if (mobit != get_options().submodels_end("mobility"))
    {
      (mobit->second).set_option("particle", "hole");
      (mobit->second).set_option("name", string("hole_mobility"));
      _hole_mobility = create_mobility_model(mobit->second);
    }
    else
    {
      ModelOptions opts;
      opts.set_option("particle", "hole");
      opts.set_option("name", string("hole_mobility"));
      _hole_mobility = create_mobility_model(opts);
    }



    //
    // Recombinations
    //
    create_recombination_models();

    //
    // Thermoelectric power
    //

    it = get_options().submodels_begin("thermoelectric_power");
    end = get_options().submodels_end("thermoelectric_power");
    if (it != end)
    {
      _thermoelectric_power =
          ThermoelectricPower::create_model("default", get_material(), it->second);

      if (_thermoelectric_power == NULL)
        throw InitFailedException("Could not create thermoelectric power model");

      add_submodel("thermoelectricpower", _thermoelectric_power);
    }
  }
  else
  {
    // a dielectric does not need all this models
    delete_submodel("thermoelectricpower");
    delete_submodel("hole_mobility");
    delete_submodel("electron_mobility");
    delete_submodel("generation");
  }
}




void
DDBulkModel::do_init(void)
{

  parse_options();

  DriftDiffusionProperties::do_init();


  //if (_is_dielectric)
  //  _background_conductivity =
  //      0.5 * get_option("background_conductivity", 1e-3 * Constants::e) / Constants::e;
  //else
    _background_conductivity =
        0.5 * get_option("background_conductivity", 1e-3 * Constants::e) / Constants::e;


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








MobilityModelInterface*
DDBulkModel::create_mobility_model(const ModelOptions& options)
{
  string model_name = options.get_name();
  model_name = options.get_option("type", model_name);
  if (model_name.empty())
    model_name = "constant";

  MobilityModelInterface* mobility_model =
    MobilityModelInterface::create(model_name, get_material(), options);

  if (mobility_model == NULL)
    throw InitFailedException("No such mobility model: " + model_name);

  add_submodel(options.get_option("name",""), mobility_model);

  return mobility_model;
}









void
DDBulkModel::reinit(const Elem* elem)
{

  if (get_element() != elem)
  {
    set_element(elem);
    set_coordinates(elem->centroid());

    // get the nodal temperatures
    get_temperature_interface().get_temperature(elem, _nodal_lattice_vt);

    // get the mean temperature on the element
    set_lattice_temperature(
        get_temperature_interface().get_temperature(elem, elem->centroid()));

     _strain_if.get_crystal_strain(elem, elem->centroid(), get_strain());

      _polarization = 0;
      for (size_t n = 0; n < _pm.size(); n++)
      {
        _pm[n]->set_strain(get_strain());
        _pm[n]->calculate(get_element(), get_coordinates());
        _polarization += _pm[n]->get_polarization();
      }
      set_polarization(_polarization);

      this->prepare_element_data();
  }

  // here we assume thermal equilibrium
  get_pd().electron_vt = get_pd().hole_vt = get_lattice_temperature();

}








void
DDBulkModel::calculate_mobilities(void)
{
  PointData& pd = get_pd();
  if (is_dielectric())
  {
    pd.electron_mobility = 0.0;
    pd.hole_mobility = 0.0;
  }
  else
  {
    pd.electron_mobility = _electron_mobility->get_mobility();
    _electron_mobility->get_derivative_grad_fermi(pd.electron_mobility_derivatives);
    pd.hole_mobility = _hole_mobility->get_mobility();
    _hole_mobility->get_derivative_grad_fermi(pd.hole_mobility_derivatives);
  }

  pd.electron_conductivity =
    pd.electron_mobility * pd.electron_density + _background_conductivity;
  pd.hole_conductivity =
    pd.hole_mobility * pd.hole_density + _background_conductivity;
}








void
DDBulkModel::calculate_equilibrium_properties(void)
{

  // call this method to properly set conduction and valence band DOS
  // and energy
  double kT = get_lattice_temperature();
  set_carrier_temperatures(kT, kT);

  BandProperties& cb = get_conduction_band();
  cb.calculate(kT);
  double Ec = cb.get_band_edge();

  BandProperties& vb = get_valence_band();
  vb.calculate(kT);
  double Ev = vb.get_band_edge();


  // for a dielectric we don't need much...
  if (is_dielectric())
  {
    set_equilibrium_fermi_level(0.5 * (Ec + Ev));
    set_intrinsic_density(sqrt(cb.get_effective_DOS() * vb.get_effective_DOS())
        * exp(-0.5 * get_band_gap() / kT));
    return;
  }

  // remember the coupling
  int coupling_bkp = get_coupling_type();
  set_coupling_type(DriftDiffusionDefs::BOTH);

  bool quantum_el = get_electrons().has_quantum_density();
  bool quantum_hl = get_holes().has_quantum_density();
  get_electrons().use_quantum_density(false);
  get_holes().use_quantum_density(false);


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


  // In some cases guess can be Inf or NaN. Then we set it to midband energy
  //if (std::isinf(guess) || std::isnan(guess))
    guess = 0.5 * (Ec + Ev);

  /*
   * We use standard Newton. This should work always, as the density
   * is a strictly monotone function of the electric potential with
   * lim_{+-infty} = +-infty .
   */

  double x = guess;
  // 1e-4 V error seems to be good enough...
  double eps = 1e-4, dens_max = 1e6;
  double error, residual_dens, y;

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

  //cerr << "***\n";
  for (unsigned int i = 0; i < 200; ++i)
  {
    set_potentials(x);
    calculate_densities();
    calculate_traps();
    calculate_ionized_dopants();

    double f = get_charge_density();
    double df_fermi[2];
    get_charge_density_derivatives(df_fermi);
    double df = -(df_fermi[0] + df_fermi[1]);

    if (f > 0) xmin = x;
    else if (f < 0) xmax = x;

    residual_dens = fabs(f);

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

    error = fabs(dx);
    //cerr << "x = " << y << " error = " << dx << " res. dens. = "
    //  << residual_dens << endl;

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
    throw SolveFailedException(os.str());
  }

  set_intrinsic_density(sqrt(get_electron_density()) * sqrt(get_hole_density()));
  set_equilibrium_n(get_electron_density());
  set_equilibrium_p(get_hole_density());

  // for a dielectric we don't need much...
  if (is_dielectric())
    set_equilibrium_fermi_level(0.0);
  else
    set_equilibrium_fermi_level(y);

  // restore original coupling
  set_coupling_type(coupling_bkp);

  get_electrons().use_quantum_density(quantum_el);
  get_holes().use_quantum_density(quantum_hl);
}





void
DDBulkModel::set_equilibrium_properties(double Ef)
{
  // remember the coupling
  int coupling_bkp = get_coupling_type();
  set_coupling_type(DriftDiffusionDefs::BOTH);

  set_equilibrium_fermi_level(Ef);
  set_potentials(Ef);
  calculate_densities();

  set_intrinsic_density(sqrt(get_electron_density()) * sqrt(get_hole_density()));
  set_equilibrium_n(get_electron_density());
  set_equilibrium_p(get_hole_density());

  // restore original coupling
  set_coupling_type(coupling_bkp);
}








void
DDBulkModel::compute_thermoelectric_powers(void)
{
  if (_thermoelectric_power != NULL)
  {
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

  }
  else
  {
    _eTEpower = 0;
    _hTEpower = 0;
  }
}

void
DDBulkModel::compute_thermoelectric_power_gradient(void)
{
   if (_thermoelectric_power != NULL)
   {
     _thermoelectric_power->set_potential_gradients(get_grad_fermi_e(),
         get_grad_fermi_h(),get_electric_field());

     _thermoelectric_power->set_temperature(get_lattice_temperature());

     _thermoelectric_power->calculate_derivatives();

     _eTEpowerGrad =  _thermoelectric_power->get_electron_thermoelectric_power_gradient();

     _hTEpowerGrad =  _thermoelectric_power->get_hole_thermoelectric_power_gradient();

    }
   else
   {
     _eTEpowerGrad = 0;
     _hTEpowerGrad = 0;
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

  double deg = std::pow(2.0, 2.0 / 3.0);

  m.info("Conduction band:");
  m.indent();
  get_conduction_band().print_info();

  ostringstream os;
  os << "Ec = " << get_conduction_band().get_band_edge()
      << ", Nc = " << get_conduction_band().get_effective_DOS() << " cm^-3\n"
      << "m_dos = " << get_conduction_band().get_effective_mass() / deg
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
      << "m_dos = " << get_valence_band().get_effective_mass() / deg
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

