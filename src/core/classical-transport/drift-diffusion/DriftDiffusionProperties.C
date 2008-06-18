// $Id$

#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "MobilityModelInterface.h"
#include "ThermoelectricPower.h"
#include "SimulationInterface.h"
#include "Material.h"
#include "Database.h"
#include "Dopant.h"
#include "Constants.h"
#include "InitFailedException.h"
#include "RotatedCrystal.h"
#include "PyroPolarization.h"

#include "elem.h"
#include "getpot.h"


#include <cmath>


using namespace std;


DriftDiffusionProperties::PointData::PointData(void)
  : //electron_conductivity_derivatives(3, 0.0),
    //hole_conductivity_derivatives(3, 0.0),
    electron_recombination_rate_derivatives(3, 0.0),
    hole_recombination_rate_derivatives(3, 0.0)
{
}



// we calculate in cm, therefore the factor 1e6
// the electron charge enters because we take k*T in electron volts
const double
DriftDiffusionProperties::_DOS_factor = pow(2.0 * M_PI * Constants::me /
      (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;




DriftDiffusionProperties::DriftDiffusionProperties(void)
  : equilibrium_fermi_level(0.0),
    _is_inhomogeneous(false),
    _pd(NULL),
    _elem(NULL),
    _statistics(TiberCad::BOLTZMANN),
    _coupling(DriftDiffusionDefs::BOTH),
    _strain(0),
    _electron_mobility(NULL),
    _hole_mobility(NULL),
    _eTEpower(0),
    _hTEpower(0),
    _eTEpowerGrad(0.0),
    _hTEpowerGrad(0.0),
    _pyropolarization(NULL),
    _polarization(3, 0.0),
    _thermoelectric_power(NULL),
    _is_dielectric(false),
    _relax_polariz(1.0),
    _generalized_einstein_relation(true)
{
  _pd = new PointData();
  _pd_stack.push(pair<PointData*, bool>(_pd, 1));
}


void
DriftDiffusionProperties::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  _is_dielectric = data("dielectric", _is_dielectric);
  permittivity = data("permittivity", 1.0);

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


  _generalized_einstein_relation =
    get_options().get_option("generalized_einstein_relation",
        _generalized_einstein_relation);

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

  PhysicalModelInterface::destroy(_pyropolarization);
  _pyropolarization = PyroPolarization::create(get_material());
  _pyropolarization->set_material(get_material());
  _pyropolarization->set_simulator_id(get_simulator_id());
  _pyropolarization->init();  

  // TODO read pyropolarization from input
  //pyro_polarization(0) = get_parameter("Px", pyro_polarization(0));
  //pyro_polarization(1) = get_parameter("Py", pyro_polarization(1));
  //pyro_polarization(2) = get_parameter("Pz", pyro_polarization(2));


  // the temperature simulation
  string temp_simul = get_options().get_option("thermal_simulation", "");
  _is_inhomogeneous |= _lattice_temp.set_simulation(temp_simul);


  // the strain simulation
  string strain_simul = get_options().get_option("strain_simulation", "");
  _is_inhomogeneous |= _strain_if.set_simulation(strain_simul);
  


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
 
   
  //
  // Thermoelectric power
  // 

  it = get_options().submodels_begin("thermoelectric_power");
  end = get_options().submodels_end("thermoelectric_power");
  if (it != end)
  {
    PhysicalModelInterface::destroy(_thermoelectric_power);

    _thermoelectric_power = ThermoelectricPower::create_model("default", it->second);

    if (_thermoelectric_power == NULL)
      throw InitFailedException("Could not create thermoelectric power model");

    _thermoelectric_power->set_material(get_material());
    _thermoelectric_power->set_driftdiffusionproperties(this);
    _thermoelectric_power->set_simulator_id(get_simulator_id());
    _thermoelectric_power->init();  
  }

}




DriftDiffusionProperties::~DriftDiffusionProperties(void)
{
  clear_recombination();
  PhysicalModelInterface::destroy(_pyropolarization);
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

  _pyropolarization->build_alloy(scA->_pyropolarization,
      scB->_pyropolarization, xa);


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
DriftDiffusionProperties::lock(PointData* pd)
{
  bool del_after_use = false;
  if (pd == NULL)
  {
    pd = new PointData();
    del_after_use = true;
  }
  
  _pd_stack.push(pair<PointData*, bool>(pd, del_after_use));
  _pd = pd;
}



void
DriftDiffusionProperties::unlock(void)
{
  pair<PointData*, bool>& top = _pd_stack.top();
  if (top.second)
    delete top.first;

  _pd_stack.pop();

  _pd = (_pd_stack.top()).first;
}





void
DriftDiffusionProperties::reinit(const Elem* elem)
{
 
  if (_elem != elem) 
  {
    _elem = elem;
    _coord = elem->centroid();

    // get the nodal temperatures
    _lattice_temp.get_temperature(elem, _nodal_lattice_vt);
    //for (int i = 0; i < _nodal_lattice_vt.size(); i++)
    //  _nodal_lattice_vt[i] *= Constants::k_B;

    // get the mean temperature on the element
    _lattice_vt = Constants::k_B *
      _lattice_temp.get_temperature(elem, elem->centroid());

    _strain_if.get_strain_data(elem, _strain, _polarization);
    
    // pyropolarization is Tensor1
    _pyropolarization->calculate_polarization(_elem, _coord, _lattice_vt);
    _polarization(0) += _pyropolarization->get_polarization()(1);
    _polarization(1) += _pyropolarization->get_polarization()(2);
    _polarization(2) += _pyropolarization->get_polarization()(3);
    _polarization *= _relax_polariz;
    this->prepare_element_data();
  }

  // here we assume thermal equilibrium
  _pd->electron_vt = _pd->hole_vt = _lattice_vt;

}




void
DriftDiffusionProperties::calculate_densities(void)
{
  //double kT = _lattice_vt;
  double kTe = _pd->electron_vt;
  double kTh = _pd->hole_vt;
  
  const BandProperties& cb = conduction_band;
  const BandProperties& vb = valence_band;

  double Ec = get_conduction_band_edge();
  double Ev = get_valence_band_edge();
 
  _electrons.set_element_and_point(_elem, _coord);
  _electrons.set_classical_parameters(cb.effective_DOS,
      Ec - _pd->electric_potential, -_pd->fermi_e, kTe);
  _pd->electron_density = _electrons.get_particle_density();
  _pd->electron_density_derivative = _electrons.get_particle_density_derivative();
  
  _holes.set_element_and_point(_elem, _coord);
  _holes.set_classical_parameters(vb.effective_DOS,
      -Ev + _pd->electric_potential, _pd->fermi_h, kTh);
  _pd->hole_density = _holes.get_particle_density();
  // TODO where to put the sign?
  _pd->hole_density_derivative = -_holes.get_particle_density_derivative();
 
}




void
DriftDiffusionProperties::calculate_ionized_dopants(void)
{
  double kT = _lattice_vt;

  double Ec = get_conduction_band_edge();
  double Ev = get_valence_band_edge();

  double arg_e = -_pd->fermi_e + _pd->electric_potential - Ec;
  double arg_h = _pd->fermi_h - _pd->electric_potential + Ev;

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
  _pd->ionized_donor_density = Nd;
  _pd->ionized_donor_density_derivative = dNd;

  it = get_material()->acceptors_begin();
  end = get_material()->acceptors_end();
  for ( ; it != end; ++it)
  {
    (*it)->calculate_doping_density(_elem, _coord);
    Na += (*it)->get_ionized_dopant_density(arg_h, kT);
    dNa -= (*it)->get_ionized_dopant_density_derivative(arg_h, kT);
  }
  _pd->ionized_acceptor_density = Na;
  _pd->ionized_acceptor_density_derivative = dNa;

}




void
DriftDiffusionProperties::calculate_net_recombination_rates(void)
{
  _pd->electron_recombination_rate = 0;
  _pd->electron_recombination_rate_derivatives[0] = 0;
  _pd->electron_recombination_rate_derivatives[1] = 0;
  _pd->hole_recombination_rate = 0;
  _pd->hole_recombination_rate_derivatives[0] = 0;
  _pd->hole_recombination_rate_derivatives[1] = 0;

  double Re, Rh;
  vector<double> dRe(3), dRh(3);

  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
  {
    (it->second)->get_net_recombination_rates(Re, Rh);
    (it->second)->get_net_recombination_rate_derivatives(dRe, dRh);

    _pd->electron_recombination_rate += Re;
    _pd->electron_recombination_rate_derivatives[0] += dRe[0];
    _pd->electron_recombination_rate_derivatives[1] += dRe[1];
    _pd->hole_recombination_rate += Rh;
    _pd->hole_recombination_rate_derivatives[0] += dRh[0];
    _pd->hole_recombination_rate_derivatives[1] += dRh[1];
  }
}




void
DriftDiffusionProperties::calculate_mobilities(void)
{
  double kT = _lattice_vt;

  // electrons
  double mue = _electron_mobility->get_mobility();
  double arg_e = -_pd->fermi_e +
    _pd->electric_potential - get_conduction_band_edge();
    
  // holes
  double muh = _hole_mobility->get_mobility();
  double arg_h = _pd->fermi_h -
    _pd->electric_potential + get_valence_band_edge();

  if (_generalized_einstein_relation)
  {
    // use generalized Einstein relations!

    if (arg_e > -10 * _lattice_vt)
    {
      double Dn = kT * mue;
      mue = Dn * get_electron_density_derivative() / get_electron_density();
    }

    if (arg_h > -10 * _lattice_vt)
    {
      double Dp = kT * muh;
      muh = -Dp * get_hole_density_derivative() / get_hole_density();
    }
  }

  _pd->electron_mobility = mue;
  _pd->hole_mobility = muh;
}




void
DriftDiffusionProperties::calculate_electro_chemical_potentials(void)
{

  if (_coupling & DriftDiffusionDefs::ELECTRONS)
  {
    double kTe = _pd->electron_vt;

    const BandProperties& cb = conduction_band;
    double Ec = get_conduction_band_edge();
    double Nc = cb.effective_DOS;

    if (_pd->electron_density > 0.0)
      _pd->fermi_e = -kTe * log(_pd->electron_density / Nc) - Ec +
        _pd->electric_potential;
    else
      _pd->fermi_e = -10.0;

    if (! _coupling & DriftDiffusionDefs::HOLES)
      _pd->fermi_h = _pd->fermi_e;
  }
  
  if (_coupling & DriftDiffusionDefs::HOLES)
  {
    double kTh = _pd->hole_vt;
  
    const BandProperties& vb = valence_band;
    double Ev = get_valence_band_edge();
    double Nv = vb.effective_DOS;
  
    if (_pd->hole_density > 0.0)
      _pd->fermi_h = kTh * log(_pd->hole_density / Nv) - Ev +
        _pd->electric_potential;
    else
      _pd->fermi_h = -10.0;

    if (! _coupling & DriftDiffusionDefs::ELECTRONS)
      _pd->fermi_e = _pd->fermi_h;
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

    double f  = _pd->hole_density - _pd->electron_density +
      _pd->ionized_donor_density - _pd->ionized_acceptor_density;
    double df = _pd->hole_density_derivative - _pd->electron_density_derivative +
      _pd->ionized_donor_density_derivative - 
      _pd->ionized_acceptor_density_derivative;

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

  intrinsic_density = sqrt(_pd->electron_density) * sqrt(_pd->hole_density);

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




void
DriftDiffusionProperties::compute_thermoelectric_powers(void)
{
  if (_thermoelectric_power != NULL)
  {
    _thermoelectric_power->set_potentials(_pd->fermi_e, _pd->fermi_h,_pd->electric_potential);
    
    double cb = get_conduction_band_edge();
    
    double vb = get_valence_band_edge();
    
    _thermoelectric_power->set_band_edges(cb, vb);
    
    _thermoelectric_power->set_temperature(_lattice_vt);
    
    _thermoelectric_power->calculate();
    
    _eTEpower = _thermoelectric_power->get_electrons_thermoelectric_power();

    _hTEpower = _thermoelectric_power->get_holes_thermoelectric_power();
    
  }
}

void
DriftDiffusionProperties::compute_thermoelectric_power_gradient(void)
{
   if (_thermoelectric_power != NULL)
   {
     _thermoelectric_power->set_potential_gradients(_grad_fermi_e,_grad_fermi_h,_electric_field);

     _thermoelectric_power->set_temperature(_lattice_vt);

     _thermoelectric_power->calculate_derivatives();

     _eTEpowerGrad =  _thermoelectric_power->get_electron_thermoelectric_power_gradient();

     _hTEpowerGrad =  _thermoelectric_power->get_hole_thermoelectric_power_gradient();

    }
}

std::vector<double>&
DriftDiffusionProperties::get_temperature_at_nodes()
{
  return _nodal_lattice_vt;
}
  


void
DriftDiffusionProperties::do_print_info(void)
{
  string space("    ");
  if (_strain_if.has_simulation())
    cout << space << "using strain simulation: " <<
      _strain_if.get_simulation()->get_name() << endl;
  else if (trace(_strain) == 0.0)
    cout << space << "unstrained" << endl;
  else
    cout << space << "using strain from input file";

  if (_lattice_temp.has_simulation())
    cout << space << "using temperature simulation: " <<
      _lattice_temp.get_simulation()->get_name() << endl;

}
