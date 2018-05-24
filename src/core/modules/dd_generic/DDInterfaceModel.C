// $Id: DDInterfaceModel.C 4135 2015-09-25 10:19:38Z maufder $

#include "DDInterfaceModel.h"
#include "DDBulkModel.h"
#include "Material.h"
#include "Alloy.h"
#include "MaterialBoundary.h"
#include "Trap.h"
#include "FowlerNordheim.h"
#include "RecombinationModelInterface.h"
#include "SimulationInterface.h"
#include "ModelErrorException.h"
#include "Variable.h"
#include "elem.h"

using namespace std;


DDInterfaceModel::DDInterfaceModel(const ModelOptions& options) :
  DriftDiffusionProperties(options),
  _internal_bd(false),
  _has_current(false),
  _ref_fermi_e(0),
  _ref_fermi_h(0),
  _emission(NULL),
  _eflux(0.0),
  _eflux_sim(NULL),
  _flux_predictor(false),
  _eflux_controlled(false),
  _ddprop_A(NULL),
  _ddprop_B(NULL)
{

}


DriftDiffusionProperties*
DDInterfaceModel::get_bulk_dd_properties(void) const
{
  DriftDiffusionProperties* ddprop = NULL;

  if (get_material() != NULL)
    ddprop = static_cast<DriftDiffusionProperties*>(
        get_material()->get_model(get_simulator_id()));

  return ddprop;

}

/*
void
DDInterfaceModel::prepare_submodels(void)
{
  // first prepare our own models
  DriftDiffusionProperties::prepare_submodels();

  // then create the correct model for the two adjacent materials
  // (we can be sure the owner is an interface, for now)


}
*/


DDInterfaceModel*
DDInterfaceModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{
  DDInterfaceModel* model = NULL;

  string name("interface");
  if (options.get_key() == "Contact")
    name = "ohmic";

  name = options.get_option("type", name);

  if (name == "interface")
    model = dynamic_cast<DDInterfaceModel*>(
        PhysicalModelInterface::create(_create, _destroy, boundary, options));
  else
    model = dynamic_cast<DDInterfaceModel*>(
        PhysicalModelInterface::create("ddbnd_" + name, boundary, options));

  if (model == NULL)
    throw ModelErrorException("Unknown drift-diffusion "
        "interface model: " + name);

  return model;
}


void
DDInterfaceModel::do_init(void)
{
  unsigned int n_carriers = this->n_known_carriers();

  _coeff_a.resize(n_carriers + 1, 0);
  _coeff_b.resize(n_carriers + 1, NEUMANN);
  _coeff_g.resize(n_carriers + 1, 0);

  _jacobian.resize(n_carriers + 1);
  for (unsigned int i = 0; i < n_carriers + 1; i++)
    _jacobian[i].resize(n_carriers + 1, 0);

  vector<string> zeroflux;
  get_option("isolated_carriers", zeroflux);

  for (auto&& carrier : zeroflux)
  {
    _zero_flux.insert(this->get_carrier_id(carrier));
  }


  const MaterialBoundary* bnd = dynamic_cast<const MaterialBoundary*>(get_owner());
  if (bnd == NULL)
    throw ModelErrorException("DriftDiffusion boundary models can "
        "be used only on region boundaries");

  // create the bulk models for the two adjacent materials

  // NOTE: in some cases there are both materials, but one does not contain the
  //       model, becuase the associated region is not inside the simulator's regions
  //       In that case, we assure that _ddprop_A always is non-NULL.

  DDBulkModel* ddprop = NULL;
  const Material* mat = bnd->get_material_B();
  if (mat != NULL)
  {
    const PhysicalModel* model = mat->get_model(get_simulator_id());
    if (model != NULL)
    {
      ddprop = static_cast<DDBulkModel*>(model->copy());

      if (mat->is_alloy())
      {
        const Alloy* alloy = static_cast<const Alloy*>(mat);
        ddprop->init_alloy(alloy->get_component_A()->get_model(get_simulator_id()),
            model, alloy->get_molar_fraction());
      }
      else
        ddprop->init();
    }
  }
  _ddprop_B = ddprop;
  
  // this one is always present
  mat = bnd->get_material_A();
  const PhysicalModel* model = bnd->get_material_A()->get_model(get_simulator_id());
  if (model != NULL)
  {
    ddprop = static_cast<DDBulkModel*>(model->copy());
    if (mat->is_alloy())
    {
      const Alloy* alloy = static_cast<const Alloy*>(mat);
      ddprop->init_alloy(alloy->get_component_A()->get_model(get_simulator_id()),
          alloy->get_component_B()->get_model(get_simulator_id()), alloy->get_molar_fraction());
    }
    else
      ddprop->init();

    _ddprop_A = ddprop;
  }
  else
  {
    assert(_ddprop_B != NULL);
    _ddprop_A = _ddprop_B;
    _ddprop_B = NULL;
  }

  SimulationInterface* si = SimulationInterface::get_simulation(get_simulator_id());
  if (!si->includes_region(bnd->get_id_A()))
    mat = bnd->get_material_B();

  assert(mat != NULL);

  // we set a bulk material, just in case a submodel needs it
  set_material(mat);

  //set_conduction_band(&_ddprop_A->get_conduction_band());
  //set_valence_band(&_ddprop_A->get_valence_band());

  set_carrier_properties(_ddprop_A->get_carrier_properties());

  // to setup common submodels
  DriftDiffusionProperties::do_init();

  /*
  // get surface trap models
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

  // get surface recombination models
  it = submodels_begin("recombination");
  end = submodels_end("recombination");
  for ( ; it != end; ++it)
  {
    RecombinationModelInterface* rec =
        static_cast<RecombinationModelInterface*>(it->second);
    _recombination_models.insert(rec);
  }
  */


  if (get_number_of_recombination_models() > 0)
  {
    for (unsigned int i = 0; i < n_carriers; i++)
      set_type(i, NEUMANN);
  }



  if (get_option("field_emission", "") == "fowler_nordheim")
  {
    _emission = new FowlerNordheim(get_option("work_function", 4.7));
    _emission->set_velocity(get_option("emission_velocity", 1.5e8));
  }

  get_option("area_factor", "");

  /*
  // the default behaviour is to be voltage controlled, but the contact can also
  // be flux controlled
  string eflux = get_option("electron_current", "");
  if (!eflux.empty())
  {
    if (Variable::check_string(eflux))
      get_parameter("electron_current", _eflux);
    else
    {
      istringstream iss(eflux);
      if ((iss >> _eflux).fail())
      {
        _eflux = 0;
        _eflux_sim = SimulationInterface::find_simulation(eflux);
        if (_eflux_sim == NULL)
        {
          string msg("DriftDiffusion Interface: ");
          msg += "no electron current simulation '" + eflux + "' found.";
          throw InitFailedException(msg);
        }

        _flux_predictor = get_option("flux_predictor", false);

        _eflux_id = _eflux_sim->get_solution_id("eCurrentDensity");
        if (_eflux_id == INVALID_ID)
          throw InitFailedException("Module '" +
              eflux + "' does not contain solution variable 'eCurrentDensity'");
      }
      //has_current(true);
    }

    _eflux_controlled = true;

    for (unsigned int i = 0; i <= n_carriers(); i++)
      set_type(i, DIRICHLET);
  }
  */
}




void
DDInterfaceModel::reinit(const Elem* elem, int side)
{
  set_element(elem);
  _side = side;

  // 1. setup the two bulk models
  // 2. decide what should be the conduction and what the valence band
  //    TODO check if there is an explicit band model
  _ddprop_A->reinit(elem);


  map<ID, CarrierProperties*> carriers(_ddprop_A->get_carrier_properties());
  if (_ddprop_B != nullptr)
  {
    _ddprop_B->reinit(elem);

    // for two-sided case, put for each particle the one with the lower
    // reference energy
    const auto carriersB(_ddprop_B->get_carrier_properties());

    for (auto&& it : carriers)
    {
      if (carriersB.count(it.first))
      {
        CarrierProperties* propB = (carriersB.find(it.first))->second;
        CarrierProperties* propA = it.second;

        if (propA->get_band_edge() > propB->get_band_edge())
        {
          if (propA->get_charge() <= 0)
            carriers[it.first] = propB;
        }
        else
        {
          if (propA->get_charge() > 0)
            carriers[it.first] = propB;
        }
      }
    }

  }
  set_carrier_properties(carriers);


  // we can get the lattice temperature from the bulk model
  double kT = _ddprop_A->get_lattice_temperature();
  set_lattice_vt(kT);
  set_carrier_temperatures(kT);
  

}


void
DDInterfaceModel::compute()
{
  for (unsigned int i = 0; i <= n_known_carriers(); i++)
  {
    if (get_type(i) == NEUMANN)
    {
      _coeff_g[i] = 0;
      _jacobian[i].resize(0);
      _jacobian[i].resize(n_known_carriers() + 1, 0);
    }
  }
  calculate_densities();
  calculate_traps();
  do_compute();

  for (auto&& c : _zero_flux)
  {
    if (c == unknown_carrier_id)
      continue;

    set_type(c, NEUMANN);
    _coeff_g[c] = 0;
    _jacobian[c].resize(0);
    _jacobian[c].resize(n_known_carriers() + 1, 0);
  }


  // now add common stuff if needed
  const PointData& pd = get_point_data();

  // for now we cannot put doping on interfaces
  /*
  get_pd().ionized_donor_density = 0;
  get_pd().ionized_acceptor_density = 0;
  get_pd().ionized_donor_density_derivative = 0;
  get_pd().ionized_acceptor_density_derivative = 0;
  */

  /*
  // surface states
  if (!get_bulk_dd_properties()->is_dielectric())
  {
  if (get_type(0) == NEUMANN)
  {
    // NOTE we invert the signs because g = \epsilon \nabla\varphi \hat{n}
    // i.e. refers to the negative charge density
    // 2014-8-19 the above note seems misleading, the signs below are the right ones!
    //calculate_traps();
    double q = (pd.ionized_electron_traps + pd.ionized_hole_traps);
    double dq_dEfn = -pd.ionized_traps_derivative[0];
    double dq_dEfp = -pd.ionized_traps_derivative[1];
    if (is_internal_boundary())
    {
      q /= 2;
      dq_dEfn /= 2;
      dq_dEfp /= 2;
    }
    _coeff_g[0] += q;
    _jacobian[0][0] += dq_dEfn + dq_dEfp;
    _jacobian[0][1] -= dq_dEfn;
    _jacobian[0][2] -= dq_dEfp;
  }
  */

  if (get_number_of_recombination_models() > 0)
  {
    calculate_net_recombination_rates();

    int div = is_internal_boundary() ? 2 : 1;

    for (unsigned int i = 0; i < n_known_carriers(); i++)
    {

      _coeff_g[i] += pd.q_recombination_rate[i] / div;

      for (unsigned int j = 0; j <= n_known_carriers(); j++)
      {
        double dRn_dEf = pd.q_recombination_rate_derivatives[i][j] / div;
        _jacobian[i][j] += dRn_dEf;
      }
    }

  }

  /*
  }

  if (_eflux_controlled)
  {

    double flux = _eflux;
    if (this->is_internal_boundary())
      flux *= 0.5;

    if (_eflux_sim != NULL)
    {
      vector<double> data;
      // we take the flux from the neighbor element
      if (_eflux_sim->get_solution(get_element()->neighbor(_side), _eflux_id,
          data, get_coordinates()))
      {
        flux = data[0] * _normal(0) + data[1] * _normal(1) + data[2] * _normal(2);

        if (_flux_predictor)
        {
          double dV = pd.old_fermi_e - _ref_fermi_e;
          if ((abs(dV) > 1e-9) && (abs(flux) > 0))
          {
            double rho = abs(dV / flux);

            flux = (pd.fermi_e - _ref_fermi_e) / rho;
            //std::cerr << (pd.fermi_e - pd.old_fermi_e) / rho << " " << flux << std::endl;
            _jacobian[1][1] -= 1.0 / rho / Constants::e;
          }
        }

      }
    }


    _coeff_g[1] -= flux / Constants::e;

  }
  */
}

/*

void
DDInterfaceModel::_calculate_traps(double& q, double& dq_dEfn, double& dq_dEfp)
{
  DriftDiffusionProperties* ddprop = get_dd_properties();
  assert(ddprop != NULL);

  q = dq_dEfn = dq_dEfp = 0.0;

  double phi = ddprop->get_electric_potential();
  double Ec = ddprop->get_conduction_band_edge() - phi;
  double Ev = ddprop->get_valence_band_edge() - phi;
  double kT = ddprop->get_lattice_temperature();

  const DriftDiffusionProperties::PointData& pd = ddprop->get_point_data();

  double ionized_electron_traps = 0.0;
  double ionized_electron_traps_derivative = 0.0;
  if (_etraps.size() > 0)
  {
    double nt = 0, dnt = 0;
    set<Trap*>::iterator it(_etraps.begin());
    const set<Trap*>::iterator end(_etraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, phi, -pd.fermi_e, kT);
      nt += (*it)->get_ionized_density();
      dnt += (*it)->get_ionized_density_derivative();
    }

    ionized_electron_traps = nt;
    ionized_electron_traps_derivative = dnt;
  }

  double ionized_hole_traps = 0;
  double ionized_hole_traps_derivative = 0;
  if (_htraps.size() > 0)
  {
    double nt = 0, dnt = 0;
    set<Trap*>::iterator it(_htraps.begin());
    const set<Trap*>::iterator end(_htraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, phi, -pd.fermi_h, kT);
      nt += (*it)->get_ionized_density();
      dnt += (*it)->get_ionized_density_derivative();
    }

    ionized_hole_traps = nt;
    ionized_hole_traps_derivative = dnt;
  }

  q = ionized_electron_traps + ionized_hole_traps;
  dq_dEfn = ionized_electron_traps_derivative;
  dq_dEfp = ionized_hole_traps_derivative;
  if (is_internal_boundary())
  {
    q *= 0.5;
    dq_dEfn *= 0.5;
    dq_dEfp *= 0.5;
  }
}


void
DDInterfaceModel::_calculate_recombination(double rec[6])
{
  rec[0] = rec[1] = rec[2] = rec[3] = rec[4] = rec[5] = 0.0;

  double Re, Rh;
  vector<double> dRe(3), dRh(3);

  set<RecombinationModelInterface*>::iterator it(_recombination_models.begin());
  const set<RecombinationModelInterface*>::iterator end(_recombination_models.end());
  for ( ; it != end; ++it)
  {
    (*it)->get_net_recombination_rates(Re, Rh);
    (*it)->get_net_recombination_rate_derivatives(dRe, dRh);

    rec[0] += Re;
    rec[1] += dRe[0];
    rec[2] += dRe[1];
    rec[3] += Rh;
    rec[4] += dRh[0];
    rec[5] += dRh[1];
  }

  if (is_internal_boundary())
  {
    rec[0] /= 2;
    rec[1] /= 2;
    rec[2] /= 2;
    rec[3] /= 2;
    rec[4] /= 2;
    rec[5] /= 2;
  }
}
*/
