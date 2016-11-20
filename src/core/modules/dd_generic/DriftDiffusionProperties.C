// $Id: DriftDiffusionProperties.C 4238 2016-04-22 19:39:02Z maufder $

#include "DriftDiffusionProperties.h"
#include "ParticleDensity.h"
#include "RecombinationModelInterface.h"
#include "MobilityModelInterface.h"
#include "SimulationInterface.h"
#include "Material.h"
#include "Database.h"
#include "Dopant.h"
#include "Trap.h"
#include "Particle.h"
#include "Constants.h"
#include "InitFailedException.h"
#include "RotatedCrystal.h"
#include "Embracing.h"
#include "Messages.h"
#include "TypeDefs.h"

#include "elem.h"


#include <cmath>


using namespace std;




DriftDiffusionProperties::PointData::PointData(void)
  : //electron_conductivity_derivatives(3, 0.0),
    //hole_conductivity_derivatives(3, 0.0),
    electron_recombination_rate(0.0),
    electron_recombination_rate_derivatives(4, 0.0),
    hole_recombination_rate(0.0),
    hole_recombination_rate_derivatives(4, 0.0)
{
}






DriftDiffusionProperties::DriftDiffusionProperties(const ModelOptions& options)
  : PhysicalModel(options),
    _is_inhomogeneous(false),

    _pd(NULL),
    _elem(NULL),
    _coupling(DriftDiffusionDefs::BOTH),
    _strain(0),
    _conduction_band(NULL),
    _valence_band(NULL),
    _equilibrium_fermi_level(0.0),
    _intrinsic_density(1e10),
    _equilibrium_n(0.0),
    _equilibrium_p(0.0),
    //_eTEpowerGrad(0.0),
    //_hTEpowerGrad(0.0),
    //_eTEpower(0),
    //_hTEpower(0),
    //_polarization(0),
    _permittivity(0),
    //_thermoelectric_power(NULL),


    _is_dielectric(false)
    //_relax_polariz(1)
{
  _pd = new PointData();
}



DriftDiffusionProperties*
DriftDiffusionProperties::create(const std::string& name, const Material* mat,
    const ModelOptions& options)
{
  return PhysicalModelInterface::create<DriftDiffusionProperties>("ddbulk_" + name, mat, options);
}








void
DriftDiffusionProperties::prepare_submodels(void)
{

  /*
  if (get_options().has_submodel("carrier"))
  {
    ModelOptions::const_submodel_iterator it(get_options().submodels_begin("carrier"));

    const ModelOptions& opts = it->second;
    if (!opts.find_option("particle"))
    {
      ModelOptions elopts(opts);
      elopts["particle"] = "el";
      ModelOptions hlopts(elopts);
      hlopts["particle"] = "hl";
      if (elopts.find_option("Ec"))
      {
        elopts["band_edge"] = elopts["Ec"];
        elopts.delete_option("Ec");
        hlopts.delete_option("Ec");
      }
      if (elopts.find_option("Nc"))
      {
        elopts["effective_DOS"] = elopts["Nc"];
        elopts.delete_option("Nc");
        hlopts.delete_option("Nc");
      }
      if (elopts.find_option("m_dos_e"))
      {
        elopts["DOS_mass"] = elopts["m_dos_e"];
        elopts.delete_option("m_dos_e");
        hlopts.delete_option("m_dos_e");
      }
      if (elopts.find_option("band_gap"))
        if (elopts.find_option("Ev"))
          elopts["reference_energy"] = elopts["Ev"];

      if (hlopts.find_option("Ev"))
      {
        hlopts["band_edge"] = hlopts["Ev"];
        elopts.delete_option("Ev");
        hlopts.delete_option("Ev");
      }
      if (hlopts.find_option("Nv"))
      {
        hlopts["effective_DOS"] = hlopts["Nv"];
        elopts.delete_option("Nv");
        hlopts.delete_option("Nv");
      }
      if (hlopts.find_option("m_dos_h"))
      {
        hlopts["DOS_mass"] = hlopts["m_dos_h"];
        elopts.delete_option("m_dos_h");
        hlopts.delete_option("m_dos_h");
      }

      if (++it != get_options().submodels_end("carrier"))
        throw InitFailedException("Multiple definition of band properties for material "
            + get_material()->get_name());

      get_options().delete_submodels("carrier");
      get_options().add_submodel("carrier", elopts);
      get_options().add_submodel("carrier", hlopts);
    }
  }

  vector<CarrierProperties*> bp;
  PhysicalModelInterface::create_submodels(bp, "carrier");

  if (bp.size() > 2)
    throw InitFailedException("Multiple definition of band properties for material "
        + get_material()->get_name());

  for (size_t i = 0; i < bp.size(); i++)
  {
    if (string("el") == bp[i]->get_options().get_option("particle", ""))
      _conduction_band = bp[i];
    else if (string("hl") == bp[i]->get_options().get_option("particle", ""))
      _valence_band = bp[i];
    else
      throw InitFailedException("Unknown particle for Drift-Diffusion model");
  }

  if (get_options().has_submodel("conduction_band"))
  {
    if (_conduction_band != NULL)
    throw InitFailedException("Multiple definition of conduction band properties for material "
        + get_material()->get_name());

    ModelOptions opts(get_options().submodels_begin("conduction_band")->second);
    opts.set_option("particle", "el");
    create_submodel(_conduction_band, "carrier", opts);
  }




  if (get_options().has_submodel("valence_band"))
  {
    if (_valence_band != NULL)
    throw InitFailedException("Multiple definition of valence band properties for material "
        + get_material()->get_name());

    ModelOptions opts(get_options().submodels_begin("valence_band")->second);
    opts.set_option("particle", "hl");
    create_submodel(_valence_band, "carrier", opts);
  }

  */

  if (get_options().has_submodel("carrier"))
  { 
    vector<CarrierProperties*> cp;
    PhysicalModelInterface::create_submodels(cp, "carrier");

    for (auto icp: cp)
    {
      std::string cp_name = icp->get_name();

      if (_carrier_properties.find(cp_name) == _carrier_properties.end())
        _carrier_properties.insert(make_pair(cp_name, icp));
      else
        throw ModelErrorException("Carrier models MUST have unique names. Name: '" + cp_name + "' is not unique" );
    }
  }

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


  // polarization models
  //create_submodels(_pm, "polarization");


  // traps
  vector<PhysicalModelInterface*> pd;
  create_submodels(pd, "trap");


  if (!is_dielectric())
  {


    //
    // Recombinations
    //
    create_recombination_models();

    //
    // Thermoelectric power
    //
    /*
    ModelOptions::submodel_iterator
      it(get_options().submodels_begin("thermoelectric_power"));
    ModelOptions::submodel_iterator
      end(get_options().submodels_end("thermoelectric_power"));
    if (it != end)
    {
      _thermoelectric_power =
          ThermoelectricPower::create_model("default", get_material(), it->second);

      if (_thermoelectric_power == NULL)
        throw InitFailedException("Could not create thermoelectric power model");

      add_submodel("thermoelectricpower", _thermoelectric_power);
    }
    */
  }
  else
  {
    // a dielectric does not need all this models
    delete_submodel("recombination");
    delete_submodel("generation");
    recomb_iterator it = _recombination_models.begin();
    recomb_iterator end = _recombination_models.end();
    for ( ; it != end; ++it)
      PhysicalModelInterface::destroy(it->second);
    _recombination_models.clear();
  }


}


void
DriftDiffusionProperties::create_recombination_models(void)
{
  // we create them only if they do not exist yet
  vector<ID> ids;
  get_net_recombination_rate_IDs(ids);
  if (!is_dielectric() && !ids.size())
  {
    //
    // we can have several models!
    //

    list<ModelOptions> newopts;

    // for each trap we add an SRH recombination model
    ModelOptions::submodel_iterator it(get_options().submodels_begin("trap"));
    ModelOptions::submodel_iterator end(get_options().submodels_end("trap"));
    for (; it != end; ++it)
    {
      ModelOptions opts(it->second);
      if (opts.get_option("recombination_center", false))
      {
        it->second.delete_option("recombination_center");
        opts.set_option("trap", true);
        opts.set_option("type", "srh");
        opts.set_option("name", "recombination");
        opts.set_key("recombination");
        newopts.insert(newopts.end(), opts);
      }
    }

    it = get_options().submodels_begin("generation");
    end = get_options().submodels_end("generation");
    for ( ; it != end; ++it)
    {
      ModelOptions opts(it->second);
      opts.set_option("name", "recombination");
      opts.set_key("recombination");
      newopts.insert(newopts.end(), opts);
    }
    // we must get rid of it in the options, otherwise we get
    // trouble when copying
    get_options().delete_submodels("generation");

    list<ModelOptions>::iterator lit(newopts.begin());
    list<ModelOptions>::iterator lend(newopts.end());
    for ( ; lit != lend; ++lit)
      get_options().add_submodel((*lit).get_key(), *lit);


    vector<RecombinationModelInterface*> pd;
    create_submodels(pd, "recombination");

    for (int i = 0; i < pd.size(); ++i)
      _recombination_models.insert(make_pair(pd[i]->get_id(), pd[i]));

  }
}


void
DriftDiffusionProperties::do_init(void)
{
  // hand the temperature interface over to the band parameter models
  //_conduction_band->set_temperature_interface(_lattice_temp);
  //_valence_band->set_temperature_interface(_lattice_temp);

  SubmodelIterator it = submodels_begin("trap");
  SubmodelIterator end = submodels_end("trap");
  for ( ; it != end; ++it)
  {
    Trap* t = static_cast<Trap*>(it->second);

    if (t->get_particle() == 'e')
      _etraps.insert(t);
    else if (t->get_particle() == 'h')
      _htraps.insert(t);
  }


  
  // calculate the equilibrium
  set_lattice_temperature(SimulationOptions::T);

}



void
DriftDiffusionProperties::set_conduction_band(CarrierProperties* cb)
{
  _conduction_band = cb;
}


void
DriftDiffusionProperties::set_valence_band(CarrierProperties* vb)
{
  _valence_band = vb;
}


DriftDiffusionProperties::~DriftDiffusionProperties(void)
{
  delete _pd;
}



bool
DriftDiffusionProperties::has_solution(void) const
{
  return SimulationInterface::get_simulation(get_simulator_id())->is_solved();
}






void
DriftDiffusionProperties::add_recombination_model(
    const string& model_name, const ModelOptions& options)
{
  RecombinationModelInterface* model =
    RecombinationModelInterface::create(model_name, get_material(), options);

  if (model == NULL)
    throw InitFailedException("No such recombination model: " + model_name);

  ID id = model->get_id();
  _recombination_models.insert(make_pair(id, model));

  add_submodel("recombination", model);
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






/*
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
*/



void
DriftDiffusionProperties::calculate_densities(void)
{
  //const std::map<std::string, CarrierProperties*>& carriers = get_carrier_properties();
  for (auto&& cp: get_carrier_properties())
  {
    cp.second->set_temperature(_pd->carrier_vt[cp.first]);
    pair<double, double> dens_der(cp.second->get_density_and_derivative(_pd->fermi_potential[cp.first], _pd->electric_potential));
    _pd->q_density[cp.first] = dens_der.first;
    _pd->q_density_derivative[cp.first] = dens_der.second;
  }


  /*
  CarrierProperties& cb = get_conduction_band();
  CarrierProperties& vb = get_valence_band();
  double kTe = _pd->electron_vt;
  double kTh = _pd->hole_vt;
  cb.set_temperature(kTe);
  vb.set_temperature(kTh);

  pair<double, double> el(cb.get_density_and_derivative(_pd->fermi_e, _pd->electric_potential));
  _pd->electron_density = el.first;
  _pd->electron_density_derivative = el.second;

  pair<double, double> hl(vb.get_density_and_derivative(_pd->fermi_h, _pd->electric_potential));
  _pd->hole_density = hl.first;
  _pd->hole_density_derivative = hl.second;
  
  
  _pd->gamma_n = cb.get_gamma();
  _pd->gamma_p = vb.get_gamma();

  if (!(_coupling & DriftDiffusionDefs::ELECTRONS))
  {
    _pd->electron_density = 0.0;
    _pd->electron_density_derivative = 0.0;
  }

  if (!(_coupling & DriftDiffusionDefs::HOLES))
  {
    _pd->hole_density = 0.0;
    _pd->hole_density_derivative = 0.0;
  }
  */
}


void
DriftDiffusionProperties::calculate_traps(void)
{
  _pd->ionized_electron_traps = 0.0;
  _pd->ionized_hole_traps = 0.0;

  vector<double>& dntdEf = _pd->ionized_traps_derivative;
  dntdEf.resize(2);
  dntdEf[0] = dntdEf[1] = 0.0;
  /*
  double Ec = get_conduction_band_edge() - _pd->electric_potential;
  double Ev = get_valence_band_edge() - _pd->electric_potential;
  double phi =  _pd->electric_potential;

  if (_etraps.size() > 0)
  {
    double nt = 0;
    std::vector<double> derivatives;
    set<Trap*>::iterator it(_etraps.begin());
    const set<Trap*>::iterator end(_etraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, phi);
      Particle el(-1, _pd->electron_density, _pd->fermi_e, _pd->electron_vt);
      Particle hl(1, _pd->hole_density, _pd->fermi_h, _pd->hole_vt);
      nt += (*it)->get_ionized_density_and_derivative(_elem, _coord, el, hl, derivatives);
      // the negative sign is because the derivative is given with respect
      // to the quasi fermi level, not the electrochemical potential.
      dntdEf[0] -=
          derivatives[0] * _pd->electron_density_derivative +
          derivatives[2] + derivatives[4];
      dntdEf[1] -=
          derivatives[1] * _pd->hole_density_derivative +
          derivatives[3];
    }

    _pd->ionized_electron_traps = nt;
  }

  if (_htraps.size() > 0)
  {
    double nt = 0;
    std::vector<double> derivatives;
    set<Trap*>::iterator it(_htraps.begin());
    const set<Trap*>::iterator end(_htraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, phi);
      Particle el(-1, _pd->electron_density, _pd->fermi_e, _pd->electron_vt);
      Particle hl(1, _pd->hole_density, _pd->fermi_h, _pd->hole_vt);
      nt += (*it)->get_ionized_density_and_derivative(_elem, _coord, el, hl, derivatives);
      dntdEf[0] -=
          derivatives[0] * _pd->electron_density_derivative +
          derivatives[2];
      dntdEf[1] -=
          derivatives[1] * _pd->hole_density_derivative +
          derivatives[3] + derivatives[4];
    }

    _pd->ionized_hole_traps = nt;
  }
  */
}


void
DriftDiffusionProperties::calculate_ionized_dopants(void)
{
  /*
  double kT = _lattice_vt;

  double Ec = get_conduction_band_edge();
  double Ev = get_valence_band_edge();

  double arg_e = -_pd->fermi_e + _pd->electric_potential - Ec;
  double arg_h = _pd->fermi_h - _pd->electric_potential + Ev;
  */

  double Nd = 0, dNd = 0;
  double Na = 0, dNa = 0;

  /*
  Material::dopant_iterator it(get_material()->donors_begin());
  Material::dopant_iterator end(get_material()->donors_end());
  for ( ; it != end; ++it)
  {
    (*it)->calculate_doping_density(_elem, _coord);
    Nd += (*it)->get_ionized_dopant_density(arg_e, kT);
    dNd += (*it)->get_ionized_dopant_density_derivative(arg_e, kT);
  }
  */
  _pd->ionized_donor_density = Nd;
  _pd->ionized_donor_density_derivative = dNd;

  /*
  it = get_material()->acceptors_begin();
  end = get_material()->acceptors_end();
  for ( ; it != end; ++it)
  {
    (*it)->calculate_doping_density(_elem, _coord);
    Na += (*it)->get_ionized_dopant_density(arg_h, kT);
    dNa -= (*it)->get_ionized_dopant_density_derivative(arg_h, kT);
  }
  */
  _pd->ionized_acceptor_density = Na;
  _pd->ionized_acceptor_density_derivative = dNa;

}




void
DriftDiffusionProperties::calculate_net_recombination_rates(void)
{
  for (auto&& cpi: get_carrier_properties())
  {
    _pd->q_recombination_rate[cpi.first] = 0;
    for (auto && cpj : get_carrier_properties())
    { 
      vector<double> zero(2);
      _pd->q_recombination_rate_derivatives[cpi.first][cpj.first] = zero;
      _pd->q_recombination_rate_derivatives[cpi.first][cpj.first] = zero;
      _pd->q_recombination_rate_derivatives[cpi.first][cpj.first] = zero;
      _pd->q_recombination_rate_derivatives[cpi.first][cpj.first] = zero;
    }
  }

  double R0, R1;
  vector<double> dR0(4), dR1(4);

  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
  {
    (it->second)->get_net_recombination_rates(R0, R1);
    (it->second)->get_net_recombination_rate_derivatives(dR0, dR1);

    string cp0 = (it->second)->get_carriers()[0];
    string cp1 = (it->second)->get_carriers()[1];

    _pd->q_recombination_rate[cp0] += R0;
    _pd->q_recombination_rate_derivatives[cp0][cp0][0] += dR0[0];
    _pd->q_recombination_rate_derivatives[cp0][cp1][0] += dR0[1];
    _pd->q_recombination_rate_derivatives[cp0][cp0][1] += dR0[2];
    _pd->q_recombination_rate_derivatives[cp0][cp1][1] += dR0[3];
    _pd->q_recombination_rate[cp1] += R1;
    _pd->q_recombination_rate_derivatives[cp1][cp1][0] += dR1[0];
    _pd->q_recombination_rate_derivatives[cp1][cp0][0] += dR1[1];
    _pd->q_recombination_rate_derivatives[cp1][cp1][1] += dR1[2];
    _pd->q_recombination_rate_derivatives[cp1][cp0][1] += dR1[3];
  }
}






/*
double 
DriftDiffusionProperties::get_electron_mobility_derivative_potential(void) const
{
  if (!is_dielectric())
    _electron_mobility->get_derivative_potential();
}

double 
DriftDiffusionProperties::get_hole_mobility_derivative_potential(void) const
{
  if (!is_dielectric())
    _hole_mobility->get_derivative_potential();
}

void 
DriftDiffusionProperties::get_electron_mobility_derivative_grad_potential(libMesh::RealGradient& dmu) const
{
  if (!is_dielectric())
    _electron_mobility->get_derivative_grad_potential(dmu);
}

void 
DriftDiffusionProperties::get_hole_mobility_derivative_grad_potential(libMesh::RealGradient& dmu) const
{
  if (!is_dielectric())
    _hole_mobility->get_derivative_grad_potential(dmu);
}

void
DriftDiffusionProperties::get_electron_mobility_derivative_grad_fermi(libMesh::RealGradient& dmu) const
{
  if (!is_dielectric())
    _electron_mobility->get_derivative_grad_fermi(dmu);
}



void
DriftDiffusionProperties::get_hole_mobility_derivative_grad_fermi(libMesh::RealGradient& dmu) const
{
  if (!is_dielectric())
    _hole_mobility->get_derivative_grad_fermi(dmu);
}
*/


void
DriftDiffusionProperties::calculate_electro_chemical_potentials(void)
{

  if (_coupling & DriftDiffusionDefs::ELECTRONS)
  {
    double kTe = _pd->electron_vt;

    CarrierProperties& cb = get_conduction_band();
    cb.set_temperature(kTe);
    double Ec = get_conduction_band_edge();
    double Nc = cb.get_effective_DOS();

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

    CarrierProperties& vb = get_valence_band();
    vb.set_temperature(kTh);
    double Ev = get_valence_band_edge();
    double Nv = vb.get_effective_DOS();

    if (_pd->hole_density > 0.0)
      _pd->fermi_h = kTh * log(_pd->hole_density / Nv) - Ev +
        _pd->electric_potential;
    else
      _pd->fermi_h = -10.0;

    if (! _coupling & DriftDiffusionDefs::ELECTRONS)
      _pd->fermi_e = _pd->fermi_h;
  }
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



pair<double, double>
DriftDiffusionProperties::get_net_recombination_rate(ID id)
{
  double rece = 0.0;
  double rech = 0.0;

  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
    if (it->first == id)
    {
      double re, rh;
      it->second->get_net_recombination_rates(re, rh);
      rece += re;
      rech += rh;
    }

  return(make_pair(rece, rech));
}


void
DriftDiffusionProperties::set_equilibrium_properties(double Ef)
{
  // remember the coupling
  int coupling_bkp = _coupling;
  _coupling = DriftDiffusionDefs::BOTH;

  _equilibrium_fermi_level = Ef;
  set_potentials(Ef);
  calculate_densities();

  _intrinsic_density = sqrt(_pd->electron_density) * sqrt(_pd->hole_density);
  _equilibrium_n = _pd->electron_density;
  _equilibrium_p = _pd->hole_density;

  // restore original coupling
  _coupling = coupling_bkp;
}


void
DriftDiffusionProperties::set_carrier_properties(const std::map<std::string, CarrierProperties*>& cp)
{
  _carrier_properties = cp;
}

void
DriftDiffusionProperties::set_carrier_properties(std::string name, CarrierProperties* cp)
{
  if (_carrier_properties.find(name) != _carrier_properties.end() )
    _carrier_properties[name]=cp;
  else
    _carrier_properties.insert(std::make_pair(name, cp));
}



/*
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
*/



