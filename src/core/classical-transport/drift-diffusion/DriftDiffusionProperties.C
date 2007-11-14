// $Id$

#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "MobilityModelInterface.h"
#include "ThermoelectricPower.h"
#include "Material.h"
#include "Database.h"
#include "Dopant.h"
#include "Constants.h"
#include "InitFailedException.h"
#include "RotatedCrystal.h"

#include "elem.h"
#include "getpot.h"

#include <cmath>


using namespace std;


// we calculate in cm, therefore the factor 1e6
// the electron charge enters because we take k*T in electron volts
const double
DriftDiffusionProperties::_DOS_factor = pow(2.0 * M_PI * Constants::me /
      (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;




DriftDiffusionProperties::DriftDiffusionProperties(void)
  : electron_conductivity_derivatives(3, 0.0),
    hole_conductivity_derivatives(3, 0.0),
    electron_recombination_rate_derivatives(3, 0.0),
    hole_recombination_rate_derivatives(3, 0.0),
    pyro_polarization(3, 0.0),
    bow_pyro(0.0),
    _elem(NULL),
    _statistics(TiberCad::BOLTZMANN),
    _coupling(DriftDiffusionDefs::BOTH),
    _strain(0),
    _electron_mobility(NULL),
    _hole_mobility(NULL),
    _eTEpower(0),
    _hTEpower(0),
    _polarization(3, 0.0),
    _thermoelectric_power(NULL),
    _is_dielectric(false),
    _relax_polariz(1.0)
{
}


void
DriftDiffusionProperties::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  _is_dielectric = data("dielectric", _is_dielectric);
  permittivity = data("permittivity", 1.0);

  // pyropolarization
  Tensor1 pol;
  pol(1) = pol(2) = 0.0;
  pol(3) = data("Pz", 0.0);
  pol = (mat->get_rotated_crystal()).RotMatrix * pol;
  pyro_polarization(0) = pol(1);
  pyro_polarization(1) = pol(2);
  pyro_polarization(2) = pol(3);
  bow_pyro = data("bow_Pz", 0.0);

}


void
DriftDiffusionProperties::set_variable_value(double value, ID id)
{
  //if (id == RELAXPOLARIZ)
    _relax_polariz = value;
}


double
DriftDiffusionProperties::get_variable_value(ID id)
{
  return _relax_polariz;
}



void
DriftDiffusionProperties::do_init(void)
{

  const string stat("B");
  if (get_options().get_option("statistics", stat) == "FD")
    set_statistics(TiberCad::FERMIDIRAC);

  
  std::string s(get_parameter("relax_polarization", ""));
  _relax_polariz = check_and_register(s, _relax_polariz);


  //
  // setup holes and electrons
  //
  _holes.set_particle_charge(1.0);
  _holes.set_statistics(_statistics);
  _electrons.set_statistics(_statistics);

  // we could have quantum density simulations for them
  {
    vector<string> qd;
    get_parameter("electron_quantum_density", qd);
    for (int i = 0; i < qd.size(); i++)
      _electrons.add_quantum_density(qd[i]);

    qd.resize(0);
    get_parameter("hole_quantum_density", qd);
    for (int i = 0; i < qd.size(); i++)
      _holes.add_quantum_density(qd[i]);
  }


  _is_dielectric = get_parameter("dielectric", _is_dielectric);
  permittivity = get_parameter("permittivity", permittivity);


  pyro_polarization(0) = get_parameter("Px", pyro_polarization(0));
  pyro_polarization(1) = get_parameter("Py", pyro_polarization(1));
  pyro_polarization(2) = get_parameter("Pz", pyro_polarization(2));


  // the temperature simulation
  string temp_simul = get_options().get_option("thermal_simulation", "");
  _lattice_temp.set_simulation(temp_simul);


  //
  // mobilities
  // 

  // create electron mobility model
  PhysicalModelInterface::destroy(_electron_mobility);

  ModelOptions::const_submodel_iterator
    it(get_options().submodels_begin("electron_mobility"));
  ModelOptions::const_submodel_iterator
    end(get_options().submodels_end("electron_mobility"));

  if (it != end)
    _electron_mobility = create_mobility_model(it->second);
  else
  {
    ModelOptions opts;
    opts.set_option("name", string("electron_mobility"));
    _electron_mobility = create_mobility_model(opts);
  }
  _electron_mobility->set_carrier_type('e');
  _electron_mobility->init();


  // create hole mobility model
  PhysicalModelInterface::destroy(_hole_mobility);

  it = get_options().submodels_begin("hole_mobility");
  end = get_options().submodels_end("hole_mobility");

  if (it != end)
    _hole_mobility = create_mobility_model(it->second);
  else
  {
    ModelOptions opts;
    opts.set_option("name", string("hole_mobility"));
    _hole_mobility = create_mobility_model(opts);
  }
  _hole_mobility->set_carrier_type('h');
  _hole_mobility->init();


  
  //
  // Recombinations (we can have several models!)
  //

  // create recombination models
  it = get_options().submodels_begin("recombination");
  end = get_options().submodels_end("recombination");
  for ( ; it != end; ++it)
  {
    const std::string& name = (it->second).get_option("model", "");
    add_recombination_model(name, it->second);
  }
 
   

   // create a pointer to thermoelectric power
   PhysicalModelInterface::destroy(_thermoelectric_power);
   
   it = get_options().submodels_begin("thermoelectric_power");
   end = get_options().submodels_end("thermoelectric_power");

   if (it != end)
   {
     _thermoelectric_power = dynamic_cast<ThermoelectricPower*>(
      PhysicalModelInterface::create("thermoelectric_power", it->second));

    if (_thermoelectric_power == NULL)
      throw InitFailedException("Could not create thermoelectric power model");

    _thermoelectric_power->set_material(get_material());
    _thermoelectric_power->init();  
   }

}




DriftDiffusionProperties::~DriftDiffusionProperties(void)
{
  clear_recombination();
  PhysicalModelInterface::destroy(_electron_mobility);
  PhysicalModelInterface::destroy(_hole_mobility);
  PhysicalModelInterface::destroy(_thermoelectric_power);
}




void
DriftDiffusionProperties::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  const DriftDiffusionProperties* scA =
    dynamic_cast<const DriftDiffusionProperties*>(comp_A);
  const DriftDiffusionProperties* scB =
    dynamic_cast<const DriftDiffusionProperties*>(comp_B);


  permittivity = alloy(scA->permittivity, scB->permittivity, xa);


  conduction_band.effective_mass =
    alloy(scA->conduction_band.effective_mass,
        scB->conduction_band.effective_mass, xa);
  conduction_band.effective_DOS =
    alloy(scA->conduction_band.effective_DOS,
        scB->conduction_band.effective_DOS, xa);
  conduction_band.band_edge =
    alloy(scA->conduction_band.band_edge,
        scB->conduction_band.band_edge, xa);

  valence_band.effective_mass =
    alloy(scA->valence_band.effective_mass,
        scB->valence_band.effective_mass, xa);
  valence_band.effective_DOS =
    alloy(scA->valence_band.effective_DOS,
        scB->valence_band.effective_DOS, xa);
  valence_band.band_edge =
    alloy(scA->valence_band.band_edge,
        scB->valence_band.band_edge, xa);

  _electron_mobility->build_alloy(scA->_electron_mobility,
      scB->_electron_mobility, xa);
  _hole_mobility->build_alloy(scA->_hole_mobility,
      scB->_hole_mobility, xa);

  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
  {
    ID id = it->first;
    (it->second)->build_alloy(scA->get_recombination_model(id),
                              scB->get_recombination_model(id), xa);
  }

  // pyropolarization
  pyro_polarization(0) = alloy(scA->pyro_polarization(0),
      scB->pyro_polarization(0), xa, bow_pyro);
  pyro_polarization(1) = alloy(scA->pyro_polarization(1),
      scB->pyro_polarization(1), xa, bow_pyro);
  pyro_polarization(2) = alloy(scA->pyro_polarization(2),
      scB->pyro_polarization(2), xa, bow_pyro);
}




void
DriftDiffusionProperties::add_recombination_model(
    const string& model_name, const ModelOptions& options)
{
  RecombinationModelInterface* model =
    RecombinationModelInterface::create(model_name, options);

  if (model == NULL)
    throw InitFailedException("No such recombination model: " + model_name);

  ID id = model->get_id();
  _recombination_models[id] = model;
  model->set_driftdiffusionproperties(this);
  model->set_material(get_material());
  model->set_simulator_id(get_simulator_id());
  model->init();
}

 

 
MobilityModelInterface*
DriftDiffusionProperties::create_mobility_model(const ModelOptions& options)
{
  const string& model_name = options.get_option("model", "constant");

  MobilityModelInterface* mobility_model =
    MobilityModelInterface::create(model_name, options);

  if (mobility_model == NULL)
    throw InitFailedException("No such mobility model: " + model_name);
 
  mobility_model->set_driftdiffusionproperties(this);
  mobility_model->set_material(get_material());
  mobility_model->set_simulator_id(get_simulator_id());

  return mobility_model;
}





void
DriftDiffusionProperties::clear_recombination(void)
{
  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
    PhysicalModelInterface::destroy(it->second);

  _recombination_models.clear();
}



void
DriftDiffusionProperties::reinit(const Elem* elem)
{
 
  if  ( _elem != elem) 
  {
    _elem = elem;
    _coord = elem->centroid();

    // get the nodal temperatures
    _lattice_temp.get_temperature(elem, _nodal_lattice_vt);
    //for (int i = 0; i < _nodal_lattice_vt.size(); i++)
    //  _nodal_lattice_vt[i] *= Constants::k_B;

    // get the mean temperature on the element
    lattice_vt = Constants::k_B *
      _lattice_temp.get_temperature(elem, elem->centroid());

    // here we assume thermal equilibrium
    electron_vt = hole_vt = lattice_vt;

    _polarization = 0.0;

    this->prepare_element_data();
    
    _polarization += _relax_polariz * pyro_polarization;
  }

}




void
DriftDiffusionProperties::calculate_densities(void)
{
  //double kT = lattice_vt;
  double kTe = electron_vt;
  double kTh = hole_vt;
  
  const BandProperties& cb = conduction_band;
  const BandProperties& vb = valence_band;

  double Ec = get_conduction_band_edge();
  double Ev = get_valence_band_edge();
 
  _electrons.set_element_and_point(_elem, _coord);
  _electrons.set_classical_parameters(cb.effective_DOS,
      Ec - electric_potential, -fermi_e, kTe);
  electron_density = _electrons.get_particle_density();
  electron_density_derivative = _electrons.get_particle_density_derivative();
  
  _holes.set_element_and_point(_elem, _coord);
  _holes.set_classical_parameters(vb.effective_DOS,
      -Ev + electric_potential, fermi_h, kTh);
  hole_density = _holes.get_particle_density();
  // TODO where to put the sign?
  hole_density_derivative = -_holes.get_particle_density_derivative();
 
}




void
DriftDiffusionProperties::calculate_ionized_dopants(void)
{
  double kT = lattice_vt;

  double Ec = get_conduction_band_edge();
  double Ev = get_valence_band_edge();

  double arg_e = -fermi_e + electric_potential - Ec;
  double arg_h = fermi_h - electric_potential + Ev;

  double Nd = 0, dNd = 0;
  double Na = 0, dNa = 0;
  
  Material::dopant_iterator it(get_material()->donors_begin());
  Material::dopant_iterator end(get_material()->donors_end());
  for ( ; it != end; ++it)
  {
    (*it)->calculate_doping_density(_elem, _coord);
    Nd += (*it)->get_ionized_dopant_density(arg_e, kT);
    dNd += (*it)->get_ionized_dopant_density_derivative(arg_e, kT);
  }
  ionized_donor_density = Nd;
  ionized_donor_density_derivative = dNd;

  it = get_material()->acceptors_begin();
  end = get_material()->acceptors_end();
  for ( ; it != end; ++it)
  {
    (*it)->calculate_doping_density(_elem, _coord);
    Na += (*it)->get_ionized_dopant_density(arg_h, kT);
    dNa -= (*it)->get_ionized_dopant_density_derivative(arg_h, kT);
  }
  ionized_acceptor_density = Na;
  ionized_acceptor_density_derivative = dNa;

}




void
DriftDiffusionProperties::calculate_net_recombination_rates(void)
{
  electron_recombination_rate = 0;
  electron_recombination_rate_derivatives[0] = 0;
  electron_recombination_rate_derivatives[1] = 0;
  hole_recombination_rate = 0;
  hole_recombination_rate_derivatives[0] = 0;
  hole_recombination_rate_derivatives[1] = 0;

  double Re, Rh;
  vector<double> dRe(3), dRh(3);

  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
  {
    (it->second)->get_net_recombination_rates(Re, Rh);
    (it->second)->get_net_recombination_rate_derivatives(dRe, dRh);

    electron_recombination_rate += Re;
    electron_recombination_rate_derivatives[0] += dRe[0];
    electron_recombination_rate_derivatives[1] += dRe[1];
    hole_recombination_rate += Rh;
    hole_recombination_rate_derivatives[0] += dRh[0];
    hole_recombination_rate_derivatives[1] += dRh[1];
  }
}




void
DriftDiffusionProperties::calculate_mobilities(void)
{
  //double kT = lattice_vt;

  double mue = _electron_mobility->get_mobility();
  //double electron_diffusivity = kT * mue;
  electron_mobility = mue;
  //electron_mobility = electron_diffusivity * dn_over_n;
  //electron_conductivity = electron_diffusivity * dn;
  //electron_conductivity_derivatives[0] = electron_diffusivity * dn2;
  //electron_conductivity_derivatives[1] = electron_diffusivity * dn2;
  double muh = _hole_mobility->get_mobility();
  //double hole_diffusivity = kT * muh;
  hole_mobility = muh;
  //hole_mobility = -hole_diffusivity * dp_over_p;
  //hole_conductivity = -hole_diffusivity * dp;
  //hole_conductivity_derivatives[0] = -hole_diffusivity * dp2;
  //hole_conductivity_derivatives[2] = -hole_diffusivity * dp2;
}




void
DriftDiffusionProperties::calculate_electro_chemical_potentials(void)
{

  if (_coupling & DriftDiffusionDefs::ELECTRONS)
  {
    double kTe = electron_vt;

    const BandProperties& cb = conduction_band;
    double Ec = get_conduction_band_edge();
    double Nc = cb.effective_DOS;

    if (electron_density > 0.0)
      fermi_e = -kTe * log(electron_density / Nc) - Ec + electric_potential;
    else
      fermi_e = -10.0;

    if (! _coupling & DriftDiffusionDefs::HOLES)
      fermi_h = fermi_e;
  }
  
  if (_coupling & DriftDiffusionDefs::HOLES)
  {
    double kTh = hole_vt;
  
    const BandProperties& vb = valence_band;
    double Ev = get_valence_band_edge();
    double Nv = vb.effective_DOS;
  
    if (hole_density > 0.0)
      fermi_h = kTh * log(hole_density / Nv) - Ev + electric_potential;
    else
      fermi_h = -10.0;

    if (! _coupling & DriftDiffusionDefs::ELECTRONS)
      fermi_e = fermi_h;
  }
}



// TODO
void
DriftDiffusionProperties::get_net_recombination_rates(
    vector<double>& rates)
{
}



int
DriftDiffusionProperties::get_net_recombination_rate_IDs(
    vector<ID>& ids)
{
  int n = _recombination_models.size();

  ids.resize(n);

  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  int ctr = 0;
  for ( ; it != end; ++it, ctr++)
    ids[ctr] = (it->first);

  return n;
}



double
DriftDiffusionProperties::get_net_recombination_rate(ID id)
{
  double r = 0.0, dummy;
  
  RecombinationModelInterface* rec =
    get_recombination_model(id);
  if (rec != NULL)
    rec->get_net_recombination_rates(r, dummy);

  return r;
}




void
DriftDiffusionProperties::calculate_equilibrium_properties(void)
{
  
  // call this method to properly set conduction and valence band DOS
  // and energy
  setup_band_edges();

  // remember the coupling
  int coupling_bkp = _coupling;
  _coupling = DriftDiffusionDefs::BOTH;


  double kT = get_lattice_temperature();

  const BandProperties& cb = conduction_band;
  const BandProperties& vb = valence_band;
  double Nd = get_material()->get_total_donor_density();
  double Na = get_material()->get_total_acceptor_density();

  double ni2 = cb.effective_DOS * vb.effective_DOS
    * exp(-get_band_gap() / kT);
  double ni = sqrt(ni2);
  intrinsic_density = ni;

  double guess;
  // Hmm... Is there a better guess?
  if (Nd > Na)
  {
    guess = cb.band_edge - kT *
      log(cb.effective_DOS / (Nd + ni));
  }
  else
  {
    guess = vb.band_edge + kT *
      log(vb.effective_DOS / (Na + ni));
  }

  // In some cases guess can be Inf or NaN. Then we set it to midband energy
  if (isinf(guess) || isnan(guess))
    guess = 0.5 * (cb.band_edge + vb.band_edge);

  /*
   * We use standard Newton. This should work always, as the density
   * is a strictly monotone function of the electric potential with
   * lim_{+-infty} = +-infty .
   */

  double x = guess;
  // 1e-6 V error seems to be good enough...
  double eps = 1e-6, dens_max = 1e6;
  double error, residual_dens, y;

  //set_carrier_temperatures(kT, kT);

  do
  {
    set_potentials(x);
    calculate_densities();
    calculate_ionized_dopants();

    double f  = hole_density - electron_density + ionized_donor_density -
      ionized_acceptor_density;
    double df = hole_density_derivative - electron_density_derivative +
      ionized_donor_density_derivative - ionized_acceptor_density_derivative;


    // At low temperatures everything is very sensitive on dx, so we don't
    // allow it to be bigger than k*T. At high temperatures this should not
    // have any impact
    double dx = - f / df;
    if (fabs(dx) > kT)
      if (dx > 0)
        dx = kT;
      else
        dx = -kT;

    y = x + dx;

    error = fabs(y - x);
    residual_dens = fabs(f);
    //cerr << "x = " << y << " error = " << error << " res. dens. = "
    //  << residual_dens << endl;
    
    x = y;
  }
  while ((error > eps) || (residual_dens > dens_max));

  intrinsic_density = sqrt(electron_density) * sqrt(hole_density);

  equilibrium_fermi_level =  y;
  
  // restore original coupling
  _coupling = coupling_bkp;

}


void
DriftDiffusionProperties::copy_from(const PhysicalModelInterface* rhs)
{

  const DriftDiffusionProperties* mod =
    dynamic_cast<const DriftDiffusionProperties*>(rhs);

  equilibrium_fermi_level = mod->get_equilibrium_fermi_level();
  intrinsic_density = mod->get_intrinsic_density();
  _coupling = mod->_coupling;
  _statistics = mod->_statistics;
  _strain = mod->_strain;
  conduction_band = mod->conduction_band;
  valence_band = mod->valence_band;
  
}





void  DriftDiffusionProperties::compute_thermoelectric_powers()
{

  if (_thermoelectric_power != NULL)
  {

    _thermoelectric_power->set_fermi_potential(fermi_e,fermi_h);
    
    double cb = get_conduction_band_edge()-electric_potential;
    
    double vb = get_valence_band_edge()-electric_potential;
    
    _thermoelectric_power->set_band_edges(cb,vb);
 
    _thermoelectric_power->set_mobility_term(5.0,5.0);

    
    _thermoelectric_power->set_temperature(lattice_vt / Constants::k_B);
    
    _thermoelectric_power->re_init();
    
    _eTEpower = _thermoelectric_power->get_electrons_thermoelectric_power();

    _hTEpower = _thermoelectric_power->get_holes_thermoelectric_power();
    
  }
 

}

std::vector<double> DriftDiffusionProperties::get_temperature_node()
{
  return _nodal_lattice_vt;
}
