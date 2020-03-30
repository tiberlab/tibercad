// $Id: DriftDiffusionProperties.C 4238 2016-04-22 19:39:02Z maufder $

#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
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

#include "libmesh/elem.h"


#include <cmath>


using namespace std;




DriftDiffusionProperties::PointData::PointData(void)
{
}






DriftDiffusionProperties::DriftDiffusionProperties(const ModelOptions& options)
  : PhysicalModel(options),
    _is_inhomogeneous(false),
    _pd(NULL),
    _elem(NULL),
    _strain(0),
    _equilibrium_fermi_level(0.0),
    //_polarization(0),
    _permittivity(0),
    _donor_reference_carrier(0),
    _acceptor_reference_carrier(0)
    //_thermoelectric_power(NULL),
{
  _pd = new PointData();
}



DriftDiffusionProperties*
DriftDiffusionProperties::create(const std::string& name, const Material* mat,
    const ModelOptions& options)
{
  return PhysicalModel::create<DriftDiffusionProperties>("ddbulk_" + name, mat, options);
}




void
DriftDiffusionProperties::copy_from(const PhysicalModel* rhs)
{
  const DriftDiffusionProperties* other =
      dynamic_cast<const DriftDiffusionProperties*>(rhs);

  this->set_known_carriers(other->get_known_carriers());
}




void
DriftDiffusionProperties::prepare_submodels(void)
{

  if (get_options().has_submodel("carrier"))
  { 
    // this should be empty here, actually
    _carrier_properties.clear();

    vector<CarrierProperties*> carriers;
    PhysicalModel::create_submodels(carriers, "carrier");

    // the check for duplicates
    set<string> used_names;

    for (auto cp : carriers)
    {
      string cp_name = cp->get_name();

      if (!used_names.count(cp_name))
      {
        int id = this->get_carrier_id(cp_name);
        // we should never have this case, by construction:
        if (id < 0)
          throw ModelErrorException("Carrier '" + cp_name + "' is unknown" );

        _carrier_properties[id] = cp;

        used_names.insert(cp_name);
      }
      else
        throw ModelErrorException("Carrier models must have unique names. "
            "Name: '" + cp_name + "' has already been used" );
    }
  }




  // traps
  vector<PhysicalModel*> pd;
  create_submodels(pd, "trap");


  {


    //
    // Recombinations
    //
    create_recombination_models();

  }

  _grad_fermi.resize(n_known_carriers());

}


void
DriftDiffusionProperties::create_recombination_models(void)
{
  // we create them only if they do not exist yet
  vector<ID> ids;
  get_net_recombination_rate_IDs(ids);

  if (!ids.size())
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
    {

        _recombination_models.insert(make_pair(pd[i]->get_id(), pd[i]));
    }

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

  auto rec = _recombination_models.begin();
  for ( ; rec != _recombination_models.end(); )
  {
    auto cur = rec;
    ++rec;

    const vector<ID>& ids = cur->second->get_carrier_ids();

      bool valid = true;
      for ( auto&& id : ids)
      {
        valid &= (get_carrier_properties(id) != nullptr);
      }

      if (!valid)
      {
        _recombination_models.erase(cur);
        auto it = this->submodels_begin();
        for ( ; it != this->submodels_end(); ++it)
        {
          if (it->second == cur->second)
          {
            this->delete_submodel(it);
          }
        }
        destroy(cur->second);
      }
  }

  _donor_reference_carrier = unknown_carrier_id;
  _acceptor_reference_carrier = unknown_carrier_id;
  string donor_ref = "electron";
  string acceptor_ref = "hole";

  if (get_options().has_submodel("doping"))
  {
    const ModelOptions& opts = get_options().submodels_begin("doping")->second;
    donor_ref = opts.get_option("donor_reference_carrier", donor_ref);
    acceptor_ref = opts.get_option("acceptor_reference_carrier", acceptor_ref);
  }
  donor_ref = get_option("donor_reference_carrier", donor_ref);
  acceptor_ref = get_option("acceptor_reference_carrier", acceptor_ref);

  _donor_reference_carrier = this->get_carrier_id(donor_ref);
  _acceptor_reference_carrier = this->get_carrier_id(acceptor_ref);

  if ((_donor_reference_carrier == unknown_carrier_id) &&
      (get_material()->donors_begin() != get_material()->donors_end()))
  {
    throw InitFailedException("No valid carrier specified as reference for donors.");
  }
  if ((_acceptor_reference_carrier == unknown_carrier_id) &&
      (get_material()->acceptors_begin() != get_material()->acceptors_end()))
  {
    throw InitFailedException("No valid carrier specified as reference for acceptors.");
  }

  
  // calculate the equilibrium
  set_lattice_temperature(SimulationOptions::T);

  _equilibrium_densities.resize(0);
  _equilibrium_densities.resize(n_known_carriers(), 0.0);

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
  _pd->q_density.resize(this->n_known_carriers(), 0.0);
  _pd->q_density_derivative.resize(this->n_known_carriers(), 0.0);
  _pd->charge_density_derivative.resize(this->n_known_carriers(), 0.0);

  for (auto&& cp: get_carrier_properties())
  {
    cp.second->set_temperature(_pd->carrier_vt[cp.first]);
    pair<double, double> dens_der(cp.second->get_density_and_derivative(
        _pd->fermi_potential[cp.first], _pd->electric_potential));

    if (cp.second->is_dopant())
    {
      double DOS = cp.second->get_effective_DOS();
      _pd->q_density[cp.first] = DOS - dens_der.first;
    }
    else
    {
      _pd->q_density[cp.first] = dens_der.first;
    }

    _pd->q_density_derivative[cp.first] = dens_der.second;

    double q = cp.second->get_charge();

    _pd->charge_density_derivative[cp.first] = - q * dens_der.second;
  }
}


void
DriftDiffusionProperties::calculate_traps(void)
{
  _pd->ionized_electron_traps = 0.0;
  _pd->ionized_hole_traps = 0.0;

  vector<double>& dntdEf = _pd->ionized_traps_derivative;
  dntdEf.resize(0);
  dntdEf.resize(this->n_known_carriers(), 0.0);

  if (_etraps.size() > 0)
  {
    double nt = 0;
    std::vector<double> derivatives;
    for (auto&& trap : _etraps)
    {
      // get coupled carriers
      //vector<string> carriers({"electrons", "holes"});
      vector<string> carriers;
      trap->get_options().get_option("carriers", carriers);
      if (carriers.size() != 2)
        throw InitFailedException("trap declaration needs two carriers");

      CarrierProperties* c1 = this->get_carrier_properties(carriers[0]);
      CarrierProperties* c2 = this->get_carrier_properties(carriers[1]);

      double edge1 = c1->get_band_edge();
      double edge2 = c2->get_band_edge();

      if (edge1 < edge2)
      {
        swap(c1, c2);
        swap(edge1, edge2);
      }

      ID id1 = c1->get_carrier_id();
      ID id2 = c2->get_carrier_id();

      double Ec = edge1 - c1->get_charge() * c1->get_charge_sign() * _pd->electric_potential;
      double Ev = edge2 - c2->get_charge() * c2->get_charge_sign() * _pd->electric_potential;
      double phi =  _pd->electric_potential;

      trap->set_energies(Ec, Ev, phi);
      Particle el(-1, _pd->q_density[id1], _pd->fermi_potential[id1], _pd->carrier_vt[id1]);
      Particle hl(1, _pd->q_density[id2], _pd->fermi_potential[id2], _pd->carrier_vt[id2]);
      nt += trap->get_ionized_density_and_derivative(_elem, _coord, el, hl, derivatives);

      // the negative sign is because the derivative is given with respect
      // to the quasi fermi level, not the electrochemical potential.
      dntdEf[id1] -=
          derivatives[0] * _pd->q_density_derivative[id1] +
          derivatives[2] + derivatives[4];
      dntdEf[id2] -=
          derivatives[1] * _pd->q_density_derivative[id2] +
          derivatives[3];

    }

    _pd->ionized_electron_traps = nt;

  }

  if (_htraps.size() > 0)
  {
    double nt = 0;
    std::vector<double> derivatives;
    for (auto&& trap : _htraps)
    {
      // get coupled carriers
      //vector<string> carriers({"electrons", "holes"});
      vector<string> carriers;
      trap->get_options().get_option("carriers", carriers);
      if (carriers.size() != 2)
        throw InitFailedException("trap declaration needs two carriers");

      CarrierProperties* c1 = this->get_carrier_properties(carriers[0]);
      CarrierProperties* c2 = this->get_carrier_properties(carriers[1]);

      double edge1 = c1->get_band_edge();
      double edge2 = c2->get_band_edge();

      if (edge1 < edge2)
      {
        swap(c1, c2);
        swap(edge1, edge2);
      }

      ID id1 = c1->get_carrier_id();
      ID id2 = c2->get_carrier_id();

      double Ec = edge1 - c1->get_charge() * c1->get_charge_sign() * _pd->electric_potential;
      double Ev = edge2 - c2->get_charge() * c2->get_charge_sign() * _pd->electric_potential;
      double phi =  _pd->electric_potential;
      trap->set_energies(Ec, Ev, phi);
      Particle el(-1, _pd->q_density[id1], _pd->fermi_potential[id1], _pd->carrier_vt[id1]);
      Particle hl(1, _pd->q_density[id2], _pd->fermi_potential[id2], _pd->carrier_vt[id2]);
      nt += trap->get_ionized_density_and_derivative(_elem, _coord, el, hl, derivatives);
      dntdEf[id1] -=
          derivatives[0] * _pd->q_density_derivative[id1] +
          derivatives[2];
      dntdEf[id2] -=
          derivatives[1] * _pd->q_density_derivative[id2] +
          derivatives[3] + derivatives[4];
    }

    _pd->ionized_hole_traps = nt;
  }

  for ( unsigned int i = 0; i < dntdEf.size(); ++i)
    _pd->charge_density_derivative[i] += dntdEf[i];

}


void
DriftDiffusionProperties::calculate_ionized_dopants(void)
{

  double kT = _lattice_vt;
  _pd->ionized_donor_density = 0.0;
  _pd->ionized_donor_density_derivative = 0.0;
  _pd->ionized_acceptor_density = 0.0;
  _pd->ionized_acceptor_density_derivative = 0.0;

  Material::dopant_iterator it(get_material()->donors_begin());
  Material::dopant_iterator end(get_material()->donors_end());
  for ( ; it != end; ++it)
  {
    Dopant& dop = **it;

    unsigned int ref_carrier = _donor_reference_carrier;

    if (dop.get_options().find_option("reference_carrier"))
    {
      ref_carrier = this->get_carrier_id(dop.get_options().get_option("reference_carrier", ""));
    }

    double Ec = get_carrier_band_edge(ref_carrier);
    double arg_e = -_pd->fermi_potential[ref_carrier] + _pd->electric_potential - Ec;
    double Nd = 0, dNd = 0;

    dop.calculate_doping_density(_elem, _coord);
    Nd += dop.get_ionized_dopant_density(arg_e, kT);
    dNd += dop.get_ionized_dopant_density_derivative(arg_e, kT);

    _pd->ionized_donor_density = Nd;
    _pd->ionized_donor_density_derivative = dNd;
    _pd->charge_density_derivative[ref_carrier] -= dNd;
  }

  it = get_material()->acceptors_begin();
  end = get_material()->acceptors_end();

  for ( ; it != end; ++it)
  {

    Dopant& dop = **it;

    unsigned int ref_carrier = _acceptor_reference_carrier;

    if (dop.get_options().find_option("reference_carrier"))
    {
      ref_carrier = this->get_carrier_id(dop.get_options().get_option("reference_carrier", ""));
    }

    double Ev = get_carrier_band_edge(ref_carrier);
    double arg_h = _pd->fermi_potential[ref_carrier] - _pd->electric_potential + Ev;
    double Na = 0, dNa = 0;

    dop.calculate_doping_density(_elem, _coord);
    Na += dop.get_ionized_dopant_density(arg_h, kT);
    dNa -= dop.get_ionized_dopant_density_derivative(arg_h, kT);

    _pd->ionized_acceptor_density = Na;
    _pd->ionized_acceptor_density_derivative = dNa;
    _pd->charge_density_derivative[ref_carrier] += dNa;
  }

}




void
DriftDiffusionProperties::calculate_net_recombination_rates(void)
{
  _pd->q_recombination_rate.resize(this->n_known_carriers(), 0.0);
  _pd->q_recombination_rate_derivatives.resize(this->n_known_carriers());

  for (unsigned int i = 0; i < this->n_known_carriers(); ++i)
  {
    _pd->q_recombination_rate[i] = 0;
    _pd->q_recombination_rate_derivatives[i] = vector<double>(this->n_known_carriers() + 1, 0.0);
  }

  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
  {
    vector<double> R(this->n_known_carriers(), 0.0);
    vector<vector<double>> dR( this->n_known_carriers(), vector<double>(this->n_known_carriers() + 1, 0.0) );

    (it->second)->get_net_rate_and_derivatives(R, dR);

    const vector<ID>& carriers = (it->second)->get_carrier_ids();

    for (ID i = 0; i < carriers.size(); ++i)
    {
      _pd->q_recombination_rate[carriers[i]] += R[carriers[i]];
      for (ID j = 0; j < carriers.size(); ++j)
      {
        _pd->q_recombination_rate_derivatives[carriers[i]][carriers[j]] += dR[carriers[i]][carriers[j]];
      }
      _pd->q_recombination_rate_derivatives[carriers[i]][this->n_known_carriers()] += dR[carriers[i]][this->n_known_carriers()];
    }
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
      // TODO
      double re, rh;
      //it->second->get_net_recombination_rates(re, rh);
      rece += re;
      rech += rh;
    }

  return(make_pair(rece, rech));
}

/*
void
DriftDiffusionProperties::set_equilibrium_properties(double Ef)
{

  _equilibrium_fermi_level = Ef;
  set_potentials(Ef);
  calculate_densities();

  //_intrinsic_density = sqrt(_pd->electron_density) * sqrt(_pd->hole_density);
  //_equilibrium_n = _pd->electron_density;
  //_equilibrium_p = _pd->hole_density;

}
*/


void
DriftDiffusionProperties::set_carrier_properties(const std::map<ID, CarrierProperties*>& cp)
{
  _carrier_properties = cp;
}

/*
void
DriftDiffusionProperties::set_carrier_properties(std::string name, CarrierProperties* cp)
{
  if (_carrier_properties.find(name) != _carrier_properties.end() )
    _carrier_properties[name]=cp;
  else
    _carrier_properties.insert(std::make_pair(name, cp));
}
*/




