// $Id$

#include "DriftDiffusionProperties.h"
#include "DriftDiffusion.h"
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
#include "getpot.h"


#include <cmath>


using namespace std;




DriftDiffusionProperties::PointData::PointData(void)
  : //electron_conductivity_derivatives(3, 0.0),
    //hole_conductivity_derivatives(3, 0.0),
    electron_recombination_rate(0.0),
    electron_recombination_rate_derivatives(3, 0.0),
    hole_recombination_rate(0.0),
    hole_recombination_rate_derivatives(3, 0.0)
{
}



// we calculate in cm, therefore the factor 1e6
// the electron charge enters because we take k*T in electron volts
const double
DriftDiffusionProperties::_DOS_factor = pow(2.0 * M_PI *
    Constants::me / (Constants::h * Constants::h) *
    Constants::e, 1.5) / 1e6;




DriftDiffusionProperties::DriftDiffusionProperties(const ModelOptions& options)
  : PhysicalModel(options),
    _is_inhomogeneous(false),
    _use_predictor(false),
    _driftdiffusion(NULL),
    _pd(NULL),
    _elem(NULL),
    _statistics(TiberCad::BOLTZMANN),
    _coupling(DriftDiffusionDefs::BOTH),
    _strain(0),
    _equilibrium_fermi_level(0.0),
    _intrinsic_density(1e10),
    _equilibrium_n(0.0),
    _equilibrium_p(0.0),
    _electron_mobility(NULL),
    _hole_mobility(NULL),
    _eTEpowerGrad(0.0),
    _hTEpowerGrad(0.0),
    _eTEpower(0),
    _hTEpower(0),
    _polarization(0),
    _permittivity(0),
    _background_conductivity(0.0),
    _thermoelectric_power(NULL),
    _is_dielectric(false),
    _electrons("electron"),
    _holes("hole"),
    _relax_polariz(1.0)
{
  _pd = new PointData();
}


void
DriftDiffusionProperties::read_database(void)
{
  Database& db = get_database();
  db.set_section("");
  _is_dielectric = db.get("dielectric", _is_dielectric);

}





void
DriftDiffusionProperties::parse_options(void)
{

  const string stat("B");
  if (get_options().get_option("statistics", stat) == "FD")
    set_statistics(TiberCad::FERMIDIRAC);


  get_parameter("relax_polarization", _relax_polariz);

  _use_predictor = get_option("use_density_predictor", _use_predictor);

  _is_dielectric = get_option("dielectric", _is_dielectric);


  // the temperature simulation
  string temp_simul = get_option("thermal_simulation", "");
  _is_inhomogeneous |= _lattice_temp.set_simulation(temp_simul);


  // the strain simulation
  string strain_simul = get_option("strain_simulation", "");
  _is_inhomogeneous |= _strain_if.set_simulation(strain_simul);

}



void
DriftDiffusionProperties::setup_electrons_and_holes(void)
{
  _holes.set_statistics(_statistics);
  _electrons.set_statistics(_statistics);

  // we could have quantum density simulations for them
  {
    vector<string> qd;
    get_option("electron_quantum_density", qd);
    for (size_t i = 0; i < qd.size(); i++)
      _electrons.add_quantum_density(qd[i]);

    if (_electrons.get_quantum_simulation() != NULL)
    {
      Embracing* emb =
        _driftdiffusion->create_embracing_region(
            _electrons.get_quantum_simulation(), get_options(), true);
      _electrons.set_embracing(emb);
    }

    qd.resize(0);
    get_option("hole_quantum_density", qd);
    for (size_t i = 0; i < qd.size(); i++)
      _holes.add_quantum_density(qd[i]);
    if (_holes.get_quantum_simulation() != NULL)
    {
      Embracing* emb =
        _driftdiffusion->create_embracing_region(
            _holes.get_quantum_simulation(), get_options(), true);
      _holes.set_embracing(emb);
    }
  }
}


void
DriftDiffusionProperties::create_submodels(void)
{

  // permittivity Default
  if (!get_options().has_submodel("permittivity"))
  {
    ModelOptions opts;
    opts.set_option("type", "constant");
    get_options().add_submodel("permittivity", opts);
  }
  
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
      (it->second).set_option("particle", "e");
      _electron_mobility = create_mobility_model(it->second);
    }
    else if (mobit != get_options().submodels_end("mobility"))
    {
      (mobit->second).set_option("particle", "e");
      (mobit->second).set_option("name", string("electron_mobility"));
      _electron_mobility = create_mobility_model(mobit->second);
    }
    else
    {
      ModelOptions opts;
      opts.set_option("particle", "e");
      opts.set_option("name", string("electron_mobility"));
      _electron_mobility = create_mobility_model(opts);
    }



    // create hole mobility model
    it = get_options().submodels_begin("hole_mobility");
    end = get_options().submodels_end("hole_mobility");

    if (it != end)
    {
      (it->second).set_option("particle", "h");
      _hole_mobility = create_mobility_model(it->second);
    }
    else if (mobit != get_options().submodels_end("mobility"))
    {
      (mobit->second).set_option("particle", "h");
      (mobit->second).set_option("name", string("hole_mobility"));
      _hole_mobility = create_mobility_model(mobit->second);
    }
    else
    {
      ModelOptions opts;
      opts.set_option("particle", "h");
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
          ThermoelectricPower::create_model("default", it->second);

      if (_thermoelectric_power == NULL)
        throw InitFailedException("Could not create thermoelectric power model");

      _thermoelectric_power->set_material(get_material());
      _thermoelectric_power->set_driftdiffusionproperties(this);
      add_submodel("thermoelectricpower", _thermoelectric_power);
    }
  }
  else
  {
    // a dielectric does not need all this models
    delete_submodel("thermoelectricpower");
    delete_submodel("hole_mobility");
    delete_submodel("electron_mobility");
    delete_submodel("recombination");
    delete_submodel("generation");
    _recombination_models.clear();
  }

  // eliminate them from the submodel list
  get_options().delete_submodels("mobility");
  get_options().delete_submodels("electron_mobility");
  get_options().delete_submodels("hole_mobility");
  get_options().delete_submodels("recombination");
  get_options().delete_submodels("generation");
  get_options().delete_submodels("thermoelectric_power");
}


void
DriftDiffusionProperties::create_recombination_models(void)
{

  // for each trap we add an SRH recombination model
  ModelOptions::submodel_iterator it(get_options().submodels_begin("trap"));
  ModelOptions::submodel_iterator end(get_options().submodels_end("trap"));
  for (; it != end; ++it)
  {
    ModelOptions opts(it->second);
    opts.set_option("trap", true);
    opts.set_option("type", "srh");
    opts.set_option("name", "recombination");
    get_options().add_submodel("recombination", opts);
  }


  // we create them only if they do not exist yet
  vector<ID> ids;
  get_net_recombination_rate_IDs(ids);
  if (!is_dielectric() && !ids.size())
  {
    // we can have several models!

    // create recombination models
    it = get_options().submodels_begin("recombination");
    end = get_options().submodels_end("recombination");
    for ( ; it != end; ++it)
    {
      std::string name = (it->second).get_option("model", "");
      name = (it->second).get_option("type", name);
      add_recombination_model(name, it->second);
    }

    // this is only a logical division
    it = get_options().submodels_begin("generation");
    end = get_options().submodels_end("generation");
    for ( ; it != end; ++it)
    {
      std::string name = (it->second).get_option("model", "");
      name = (it->second).get_option("type", name);
      add_recombination_model(name, it->second);
    }
  }
}


void
DriftDiffusionProperties::do_init(void)
{
  assert(_driftdiffusion != NULL);

  parse_options();

  setup_electrons_and_holes();

  SubmodelIterator it = submodels_begin("trap");
  const SubmodelIterator end = submodels_end("trap");
  for ( ; it != end; ++it)
  {
    Trap* t = static_cast<Trap*>(it->second);

    if (t->get_particle() == 'e')
      _etraps.insert(t);
    else if (t->get_particle() == 'h')
      _htraps.insert(t);
  }


  _background_conductivity =
      0.5 * get_option("background_conductivity", 2e4 * Constants::e) / Constants::e;

  // calculate the equilibrium
  set_lattice_temperature(SimulationOptions::T);
  calculate_equilibrium_properties();
  setup_band_edges();


  // Polarization 
  it = submodels_begin("polarization");
  const PhysicalModelInterface::SubmodelIterator  it_end(submodels_end("polarization"));
  for ( ; it != it_end ; ++it)
    _pm.push_back(dynamic_cast<PolarizationModel*>(it->second));
  

  // Permittivity 
  it = submodels_begin("permittivity");
  PermittivityModel* pm =  dynamic_cast<PermittivityModel*>(it->second);
  assert(pm != NULL);
  _permittivity = pm->get_permittivity();


}




DriftDiffusionProperties::~DriftDiffusionProperties(void)
{
  delete _pd;
}



bool
DriftDiffusionProperties::has_solution(void) const
{
  return _driftdiffusion->is_solved();
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
  _recombination_models.insert(make_pair(id, model));

  model->set_driftdiffusionproperties(this);
  model->set_material(get_material());

  add_submodel("recombination", model);
}




MobilityModelInterface*
DriftDiffusionProperties::create_mobility_model(const ModelOptions& options)
{
  string model_name = options.get_option("model", "constant");
  model_name = options.get_option("type", model_name);

  MobilityModelInterface* mobility_model =
    MobilityModelInterface::create(model_name, options);

  if (mobility_model == NULL)
    throw InitFailedException("No such mobility model: " + model_name);

  mobility_model->set_driftdiffusionproperties(this);
  mobility_model->set_material(get_material());

  add_submodel(options.get_option("name",""), mobility_model);

  return mobility_model;
}





void
DriftDiffusionProperties::clear_recombination(void)
{
  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
    destroy(it->second);

  _recombination_models.clear();
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

    // get the mean temperature on the element
    _lattice_vt = Constants::k_B *
      _lattice_temp.get_temperature(elem, elem->centroid());

     _strain_if.get_crystal_strain(elem, elem->centroid(), _strain);

      _polarization = 0;
      for (size_t n = 0; n < _pm.size(); n++)
      {
        _pm[n]->set_strain(_strain);
        _pm[n]->calculate(_elem, _coord);
        _polarization += _pm[n]->get_polarization();
      }
      set_polarization(_polarization);

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

  double relax = _driftdiffusion->_relaxation;

  _electrons.set_element_and_point(_elem, _coord);
  _electrons.set_classical_parameters(cb.effective_DOS,
      Ec - _pd->electric_potential, -_pd->fermi_e, kTe);
  _pd->electron_density = _electrons.get_particle_density();
  _pd->electron_density_derivative = _electrons.get_particle_density_derivative();
  _pd->gamma_n = _electrons.get_gamma();
  if (_electrons.is_quantum_density() && has_solution() && use_predictor())
  {
    _electrons.use_quantum_density(false);
    double dens = _electrons.get_particle_density();
    _pd->electron_density_derivative = _electrons.get_particle_density_derivative();
    _electrons.set_classical_parameters(cb.effective_DOS,
        Ec - _pd->old_electric_potential, -_pd->old_fermi_e, kTe);
    double old_dens = _electrons.get_particle_density();
    //if (old_dens > 1e-6)
    {
      //double fac = max(0.1, min(10.0, dens / old_dens));
      //_pd->electron_density_derivative *= _pd->electron_density / old_dens;
      //_pd->electron_density *= fac;
      double qdens = relax * _pd->electron_density + (1.0 - relax) * old_dens;
      _pd->electron_density_derivative *= qdens / old_dens;
      _pd->electron_density = qdens / old_dens * dens;
    }

    /*
    double arg_cl = std::log(dens);
    double arg_q = std::log(_pd->electron_density);
    double arg = relax * arg_q + (1.0 - relax) * arg_cl;
    _pd->electron_density = std::exp(arg);
    _pd->electron_density_derivative = _pd->electron_density / kTe;
    */

    _electrons.use_quantum_density(true);
    /* simpler but slower convergence
    {
     double arg = (_pd->electric_potential - _pd->old_electric_potential
         - _pd->fermi_e + _pd->old_fermi_e) / kTe;
     double fac = max(0.1, min(10.0, exp(arg)));
     _pd->electron_density *= fac;
     _pd->electron_density_derivative = _pd->electron_density / kTe;
    }
    */
  }

  _holes.set_element_and_point(_elem, _coord);
  _holes.set_classical_parameters(vb.effective_DOS,
      -Ev + _pd->electric_potential, _pd->fermi_h, kTh);
  _pd->hole_density = _holes.get_particle_density();
  // TODO where to put the sign?
  _pd->hole_density_derivative = -_holes.get_particle_density_derivative();
  _pd->gamma_p = _holes.get_gamma();
  if (_holes.is_quantum_density() && has_solution() && use_predictor())
  {
    _holes.use_quantum_density(false);
    double dens = _holes.get_particle_density();
    _pd->hole_density_derivative = -_holes.get_particle_density_derivative();
    _holes.set_classical_parameters(vb.effective_DOS,
        -Ev + _pd->old_electric_potential, _pd->old_fermi_h, kTh);
    double old_dens = _holes.get_particle_density();
    //if (old_dens > 1e-6)
    {
      //double fac = max(0.1, min(10.0, dens / old_dens));
      //_pd->hole_density_derivative *= _pd->hole_density / old_dens;
      //_pd->hole_density *= fac;
      double qdens = relax * _pd->hole_density + (1.0 - relax) * old_dens;
      _pd->hole_density_derivative *= qdens / old_dens;
      _pd->hole_density = qdens / old_dens * dens;
    }
    _holes.use_quantum_density(true);
    /* simpler but slower convergence
    {
     double arg = (_pd->electric_potential - _pd->old_electric_potential) / kTe;
     double fac = max(0.1, min(10.0, exp(arg)));
     _pd->electron_density *= fac;
     _pd->electron_density_derivative = _pd->electron_density / kTe;
    }
    */
  }

}


void
DriftDiffusionProperties::calculate_traps(void)
{
  double Ec = get_conduction_band_edge() - _pd->electric_potential;
  double Ev = get_valence_band_edge() - _pd->electric_potential;

  _pd->ionized_electron_traps = 0.0;
  _pd->ionized_electron_traps_derivative = 0.0;
  if (_etraps.size() > 0)
  {
    double nt = 0, dnt = 0;
    double kT = _pd->electron_vt;
    set<Trap*>::iterator it(_etraps.begin());
    const set<Trap*>::iterator end(_etraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, -_pd->fermi_e, kT);
      nt += (*it)->get_ionized_density();
      dnt += (*it)->get_ionized_density_derivative();
    }

    _pd->ionized_electron_traps = nt;
    _pd->ionized_electron_traps_derivative = dnt;
  }

  _pd->ionized_hole_traps = 0;
  _pd->ionized_hole_traps_derivative = 0;
  if (_htraps.size() > 0)
  {
    double nt = 0, dnt = 0;
    double kT = _pd->hole_vt;
    set<Trap*>::iterator it(_htraps.begin());
    const set<Trap*>::iterator end(_htraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, -_pd->fermi_h, kT);
      nt += (*it)->get_ionized_density();
      dnt += (*it)->get_ionized_density_derivative();
    }

    _pd->ionized_hole_traps = nt;
    _pd->ionized_hole_traps_derivative = dnt;
  }
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
  if (is_dielectric())
  {
    _pd->electron_mobility = 0.0;
    _pd->hole_mobility = 0.0;
  }
  else
  {
    _pd->electron_mobility = _electron_mobility->get_mobility();
    _pd->hole_mobility = _hole_mobility->get_mobility();
  }

  _pd->electron_conductivity =
    _pd->electron_mobility * _pd->electron_density + _background_conductivity;
  _pd->hole_conductivity =
    _pd->hole_mobility * _pd->hole_density + _background_conductivity;
}



void
DriftDiffusionProperties::get_electron_mobility_derivatives(RealGradient& dmu) const
{
  if (!is_dielectric())
    _electron_mobility->get_derivative_grad_fermi(dmu);
}



void
DriftDiffusionProperties::get_hole_mobility_derivatives(RealGradient& dmu) const
{
  if (!is_dielectric())
    _hole_mobility->get_derivative_grad_fermi(dmu);
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
  ignore_unused_variable(rates);
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
  double rec = 0.0;

  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
    if (it->first == id)
    {
      double r, dummy;
      it->second->get_net_recombination_rates(r, dummy);
      rec += r;
    }

  return rec;
}




void
DriftDiffusionProperties::calculate_equilibrium_properties(void)
{

  // call this method to properly set conduction and valence band DOS
  // and energy
  setup_band_edges();

  const BandProperties& cb = conduction_band;
  const BandProperties& vb = valence_band;
  double Ec = cb.band_edge;
  double Ev = vb.band_edge;

  double kT = get_lattice_temperature();

  // for a dielectric we don't need much...
  if (is_dielectric())
  {
    _equilibrium_fermi_level = 0.5 * (Ec + Ev);
    double ni2 = cb.effective_DOS * vb.effective_DOS
        * exp(-get_band_gap() / kT);
    double ni = sqrt(ni2);
    _intrinsic_density = ni;
    return;
  }

  // remember the coupling
  int coupling_bkp = _coupling;
  _coupling = DriftDiffusionDefs::BOTH;

  bool quantum_el = _electrons.has_quantum_density();
  bool quantum_hl = _holes.has_quantum_density();
  _electrons.use_quantum_density(false);
  _holes.use_quantum_density(false);



  double Nd = get_material()->get_total_donor_density();
  double Na = get_material()->get_total_acceptor_density();


  double ni2 = cb.effective_DOS * vb.effective_DOS
    * exp(-get_band_gap() / kT);
  double ni = sqrt(ni2);
  _intrinsic_density = ni;

  double guess;
  // Hmm... Is there a better guess?
  if (Nd > Na)
  {
    guess = Ec - kT * log(cb.effective_DOS / (Nd + ni));
  }
  else
  {
    guess = Ev + kT * log(vb.effective_DOS / (Na + ni));
  }


  // In some cases guess can be Inf or NaN. Then we set it to midband energy
  if (isinf(guess) || isnan(guess))
    guess = 0.5 * (Ec + Ev);


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

  //cerr << "***\n";
  for (unsigned int i = 0; i < 100; ++i)
  {
    set_potentials(x);
    calculate_densities();
    calculate_traps();
    calculate_ionized_dopants();

    double f = get_charge_density();
    double df_fermi[2];
    get_charge_density_derivatives(df_fermi);
    double df = -(df_fermi[0] + df_fermi[1]);

    residual_dens = fabs(f);

    double dx = 0.0;
    if (residual_dens > ParticleDensity::MINDENSITY)
    {
      // At low temperatures everything is very sensitive on dx, so we don't
      // allow it to be bigger than k*T. At high temperatures this should not
      // have any impact
      dx = - f / df;

      y = x + dx;
      // we limit Ef to (Ec + 0.2, Ev - 0.2)
      if ((Ec - y) < -0.2)
        dx = Ec + 0.2 - x;
      else if ((Ev - y) > 0.2)
        dx = Ev - 0.2 - x;


      //if (residual_dens > kT)
      //  if (dx > 0)
      //    dx = kT;
      //  else
      //    dx = -kT;
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

  _intrinsic_density = sqrt(_pd->electron_density) * sqrt(_pd->hole_density);
  _equilibrium_n = _pd->electron_density;
  _equilibrium_p = _pd->hole_density;

  // for a dielectric we don't need much...
  if (is_dielectric())
    _equilibrium_fermi_level = 0.0;
  else
    _equilibrium_fermi_level =  y;

  // restore original coupling
  _coupling = coupling_bkp;

  _electrons.use_quantum_density(quantum_el);
  _holes.use_quantum_density(quantum_hl);
}





void
DriftDiffusionProperties::set_equilibrium_properties(double Ef)
{
  _equilibrium_fermi_level = Ef;
  set_potentials(Ef);
  calculate_densities();

  _intrinsic_density = sqrt(_pd->electron_density) * sqrt(_pd->hole_density);
  _equilibrium_n = _pd->electron_density;
  _equilibrium_p = _pd->hole_density;
}





void
DriftDiffusionProperties::copy_from(const PhysicalModelInterface* rhs)
{

  const DriftDiffusionProperties* mod =
    dynamic_cast<const DriftDiffusionProperties*>(rhs);

  //equilibrium_fermi_level = mod->get_equilibrium_fermi_level();
  //intrinsic_density = mod->get_intrinsic_density();
  _driftdiffusion = mod->_driftdiffusion;
  //_coupling = mod->_coupling;
  //_statistics = mod->_statistics;
  //_strain = mod->_strain;
  //conduction_band = mod->conduction_band;
  //valence_band = mod->valence_band;

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
  else
  {
    _eTEpower = 0;
    _hTEpower = 0;
  }
}

void
DriftDiffusionProperties::compute_thermoelectric_power_gradient(void)
{
   if (_thermoelectric_power != NULL)
   {
     _thermoelectric_power->set_potential_gradients(_grad_fermi_e,
         _grad_fermi_h,_electric_field);

     _thermoelectric_power->set_temperature(_lattice_vt);

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
DriftDiffusionProperties::get_temperature_at_nodes()
{
  return _nodal_lattice_vt;
}



void
DriftDiffusionProperties::do_print_info(void)
{
  if (_strain_if.has_simulation())
    Messages::info("using strain simulation: " +
      _strain_if.get_simulation()->get_name());
  else if (trace(_strain) == 0.0)
    Messages::info("unstrained model");
  else
    Messages::info("using strain from input file");

  if (_lattice_temp.has_simulation())
    Messages::info("using lattice temperature from: " +
      _lattice_temp.get_simulation()->get_name());

}

