/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file DriftDiffusionProperties.C
 * \brief tiberCAD driftdiffusion module implementation.
 *
 * \note This file is part of module driftdiffusion.
 */


#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "MobilityModelInterface.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/misc/Dopant.h"
#include "tibercad/physics/misc/Trap.h"
#include "tibercad/physics/Particle.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/base/InitFailedException.h"
#include "tibercad/io/Messages.h"
#include "tibercad/base/TypeDefs.h"

#include "libmesh/elem.h"


#include <cmath>


using namespace std;

class Elem;


DriftDiffusionProperties::PointData::PointData(void)
  : //electron_conductivity_derivatives(3, 0.0),
    //hole_conductivity_derivatives(3, 0.0),
    elem(NULL),
    electron_recombination_rate(0.0),
    electron_recombination_rate_derivatives(4, 0.0),
    hole_recombination_rate(0.0),
    hole_recombination_rate_derivatives(4, 0.0)
{
}






DriftDiffusionProperties::DriftDiffusionProperties(const ModelOptions& options)
  : PhysicalModel(options),
    _is_inhomogeneous(false),
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
  _pd.push(PointData());
}



DriftDiffusionProperties*
DriftDiffusionProperties::create(const std::string& name, const Material* mat,
    const ModelOptions& options)
{
  return PhysicalModel::create<DriftDiffusionProperties>("ddbulk_" + name, mat, options);
}








void
DriftDiffusionProperties::prepare_submodels(void)
{


  if (get_options().has_submodel("band_properties"))
  {
    ModelOptions::const_submodel_iterator it(get_options().submodels_begin("band_properties"));

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

      if (++it != get_options().submodels_end("band_properties"))
        throw InitFailedException("Multiple definition of band properties for material "
            + get_material()->get_name());

      get_options().delete_submodels("band_properties");
      get_options().add_submodel("band_properties", elopts);
      get_options().add_submodel("band_properties", hlopts);
    }
  }

  vector<BandProperties*> bp;
  PhysicalModel::create_submodels(bp, "band_properties");

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
    create_submodel(_conduction_band, "band_properties", opts);
  }




  if (get_options().has_submodel("valence_band"))
  {
    if (_valence_band != NULL)
    throw InitFailedException("Multiple definition of valence band properties for material "
        + get_material()->get_name());

    ModelOptions opts(get_options().submodels_begin("valence_band")->second);
    opts.set_option("particle", "hl");
    create_submodel(_valence_band, "band_properties", opts);
  }



  // traps
  vector<PhysicalModel*> pd;
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
      PhysicalModel::destroy(it->second);
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
    Trap* t = dynamic_cast<Trap*>(it->second);

    if (t->get_particle() == 'e')
      _etraps.insert(t);
    else if (t->get_particle() == 'h')
      _htraps.insert(t);
  }



  // calculate the equilibrium
  set_lattice_temperature(SimulationOptions::T);

}


void
DriftDiffusionProperties::push_point_data(const PointData& pd)
{
  _pd.push(pd);
}

void
DriftDiffusionProperties::pop_point_data(void)
{
  _pd.pop();
  if (_pd.empty())
    _pd.push(PointData());
}


void
DriftDiffusionProperties::set_conduction_band(BandProperties* cb)
{
  _conduction_band = cb;
}


void
DriftDiffusionProperties::set_valence_band(BandProperties* vb)
{
  _valence_band = vb;
}


DriftDiffusionProperties::~DriftDiffusionProperties(void)
{
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





void
DriftDiffusionProperties::calculate_densities(void)
{
  BandProperties& cb = get_conduction_band();
  BandProperties& vb = get_valence_band();
  double kTe = get_pd().electron_vt;
  double kTh = get_pd().hole_vt;
  cb.set_temperature(kTe);
  vb.set_temperature(kTh);

  pair<double, double> el(cb.get_density_and_derivative(get_pd().fermi_e, get_pd().electric_potential));
  get_pd().electron_density = el.first;
  get_pd().electron_density_derivative = el.second;

  pair<double, double> hl(vb.get_density_and_derivative(get_pd().fermi_h, get_pd().electric_potential));
  get_pd().hole_density = hl.first;
  get_pd().hole_density_derivative = hl.second;

  if (!(_coupling & DriftDiffusionDefs::ELECTRONS))
  {
    get_pd().electron_density = 0.0;
    get_pd().electron_density_derivative = 0.0;
  }

  if (!(_coupling & DriftDiffusionDefs::HOLES))
  {
    get_pd().hole_density = 0.0;
    get_pd().hole_density_derivative = 0.0;
  }

}


void
DriftDiffusionProperties::calculate_traps(void)
{
  double Ec = get_conduction_band_edge() - get_pd().electric_potential;
  double Ev = get_valence_band_edge() - get_pd().electric_potential;
  double phi =  get_pd().electric_potential;

  get_pd().ionized_electron_traps = 0.0;
  get_pd().ionized_hole_traps = 0.0;
  vector<double>& dntdEf = get_pd().ionized_traps_derivative;
  dntdEf.resize(2);
  dntdEf[0] = dntdEf[1] = 0.0;

  if (_etraps.size() > 0)
  {
    double nt = 0;
    std::vector<double> derivatives;
    set<Trap*>::iterator it(_etraps.begin());
    const set<Trap*>::iterator end(_etraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, phi);
      Particle el(-1, get_pd().electron_density, get_pd().fermi_e, get_pd().electron_vt);
      Particle hl(1, get_pd().hole_density, get_pd().fermi_h, get_pd().hole_vt);
      nt += (*it)->get_ionized_density_and_derivative(get_pd().elem, get_pd().coord, el, hl, derivatives);
      // the negative sign is because the derivative is given with respect
      // to the quasi fermi level, not the electrochemical potential.
      dntdEf[0] -=
          derivatives[0] * get_pd().electron_density_derivative +
          derivatives[2] + derivatives[4];
      dntdEf[1] -=
          derivatives[1] * get_pd().hole_density_derivative +
          derivatives[3];
    }

    get_pd().ionized_electron_traps = nt;
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
      Particle el(-1, get_pd().electron_density, get_pd().fermi_e, get_pd().electron_vt);
      Particle hl(1, get_pd().hole_density, get_pd().fermi_h, get_pd().hole_vt);
      nt += (*it)->get_ionized_density_and_derivative(get_pd().elem, get_pd().coord, el, hl, derivatives);
      dntdEf[0] -=
          derivatives[0] * get_pd().electron_density_derivative +
          derivatives[2];
      dntdEf[1] -=
          derivatives[1] * get_pd().hole_density_derivative +
          derivatives[3] + derivatives[4];
    }

    get_pd().ionized_hole_traps = nt;
  }
}


void
DriftDiffusionProperties::calculate_ionized_dopants(void)
{
  double kT = _lattice_vt;

  double Ec = get_conduction_band_edge();
  double Ev = get_valence_band_edge();

  double arg_e = -get_pd().fermi_e + get_pd().electric_potential - Ec;
  double arg_h = get_pd().fermi_h - get_pd().electric_potential + Ev;

  double Nd = 0, dNd = 0;
  double Na = 0, dNa = 0;

  Material::dopant_iterator it(get_material()->donors_begin());
  Material::dopant_iterator end(get_material()->donors_end());
  for ( ; it != end; ++it)
  {
    (*it)->calculate_doping_density(get_pd().elem, get_pd().coord);
    Nd += (*it)->get_ionized_dopant_density(arg_e, kT);
    dNd += (*it)->get_ionized_dopant_density_derivative(arg_e, kT);
  }
  get_pd().ionized_donor_density = Nd;
  get_pd().ionized_donor_density_derivative = dNd;

  it = get_material()->acceptors_begin();
  end = get_material()->acceptors_end();
  for ( ; it != end; ++it)
  {
    (*it)->calculate_doping_density(get_pd().elem, get_pd().coord);
    Na += (*it)->get_ionized_dopant_density(arg_h, kT);
    dNa -= (*it)->get_ionized_dopant_density_derivative(arg_h, kT);
  }
  get_pd().ionized_acceptor_density = Na;
  get_pd().ionized_acceptor_density_derivative = dNa;

}




void
DriftDiffusionProperties::calculate_net_recombination_rates(void)
{
  get_pd().electron_recombination_rate = 0;
  get_pd().electron_recombination_rate_derivatives[0] = 0;
  get_pd().electron_recombination_rate_derivatives[1] = 0;
  get_pd().electron_recombination_rate_derivatives[2] = 0;
  get_pd().electron_recombination_rate_derivatives[3] = 0;
  get_pd().hole_recombination_rate = 0;
  get_pd().hole_recombination_rate_derivatives[0] = 0;
  get_pd().hole_recombination_rate_derivatives[1] = 0;
  get_pd().hole_recombination_rate_derivatives[2] = 0;
  get_pd().hole_recombination_rate_derivatives[3] = 0;

  double Re, Rh;
  vector<double> dRe(4), dRh(4);

  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it)
  {
    (it->second)->get_net_recombination_rates(Re, Rh);
    (it->second)->get_net_recombination_rate_derivatives(dRe, dRh);

    get_pd().electron_recombination_rate += Re;
    get_pd().electron_recombination_rate_derivatives[0] += dRe[0];
    get_pd().electron_recombination_rate_derivatives[1] += dRe[1];
    get_pd().electron_recombination_rate_derivatives[2] += dRe[2];
    get_pd().electron_recombination_rate_derivatives[3] += dRe[3];
    get_pd().hole_recombination_rate += Rh;
    get_pd().hole_recombination_rate_derivatives[0] += dRh[0];
    get_pd().hole_recombination_rate_derivatives[1] += dRh[1];
    get_pd().hole_recombination_rate_derivatives[2] += dRh[2];
    get_pd().hole_recombination_rate_derivatives[3] += dRh[3];
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
DriftDiffusionProperties::get_electron_mobility_derivative_grad_potential(RealGradient& dmu) const
{
  if (!is_dielectric())
    _electron_mobility->get_derivative_grad_potential(dmu);
}

void 
DriftDiffusionProperties::get_hole_mobility_derivative_grad_potential(RealGradient& dmu) const
{
  if (!is_dielectric())
    _hole_mobility->get_derivative_grad_potential(dmu);
}

void
DriftDiffusionProperties::get_electron_mobility_derivative_grad_fermi(RealGradient& dmu) const
{
  if (!is_dielectric())
    _electron_mobility->get_derivative_grad_fermi(dmu);
}



void
DriftDiffusionProperties::get_hole_mobility_derivative_grad_fermi(RealGradient& dmu) const
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
    double kTe = get_pd().electron_vt;

    BandProperties& cb = get_conduction_band();
    cb.set_temperature(kTe);
    double Ec = get_conduction_band_edge();
    double Nc = cb.get_effective_DOS();

    if (get_pd().electron_density > 0.0)
      get_pd().fermi_e = -kTe * log(get_pd().electron_density / Nc) - Ec +
        get_pd().electric_potential;
    else
      get_pd().fermi_e = -10.0;

    if (! _coupling & DriftDiffusionDefs::HOLES)
      get_pd().fermi_h = get_pd().fermi_e;
  }

  if (_coupling & DriftDiffusionDefs::HOLES)
  {
    double kTh = get_pd().hole_vt;

    BandProperties& vb = get_valence_band();
    vb.set_temperature(kTh);
    double Ev = get_valence_band_edge();
    double Nv = vb.get_effective_DOS();

    if (get_pd().hole_density > 0.0)
      get_pd().fermi_h = kTh * log(get_pd().hole_density / Nv) - Ev +
        get_pd().electric_potential;
    else
      get_pd().fermi_h = -10.0;

    if (! _coupling & DriftDiffusionDefs::ELECTRONS)
      get_pd().fermi_e = get_pd().fermi_h;
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

  _intrinsic_density = sqrt(get_pd().electron_density) * sqrt(get_pd().hole_density);
  _equilibrium_n = get_pd().electron_density;
  _equilibrium_p = get_pd().hole_density;

  // restore original coupling
  _coupling = coupling_bkp;
}







/*
void
DriftDiffusionProperties::compute_thermoelectric_powers(void)
{
  if (_thermoelectric_power != NULL)
  {
    _thermoelectric_power->set_potentials(get_pd().fermi_e, get_pd().fermi_h,get_pd().electric_potential);

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



