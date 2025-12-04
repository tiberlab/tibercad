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
 * \file DDInterfaceModel.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "DDInterfaceModel.h"
#include "DDBulkModel.h"
#include "tibercad/physics/Material.h"
#include "tibercad/physics/Alloy.h"
#include "tibercad/physics/MaterialBoundary.h"
#include "tibercad/physics/misc/Trap.h"
#include "tibercad/physics/misc/FowlerNordheim.h"
#include "RecombinationModelInterface.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/base/ModelErrorException.h"
#include "tibercad/base/Variable.h"
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
        PhysicalModel::create(_create, _destroy, boundary, options));
  else
    model = dynamic_cast<DDInterfaceModel*>(
        PhysicalModel::create("ddbnd_" + name, boundary, options));

  if (model == NULL)
    throw ModelErrorException("Unknown drift-diffusion "
        "interface model: " + name);

  return model;
}


void
DDInterfaceModel::set_type(unsigned int var, BCType type)
{
  // prepare, if not already done
  if (_coeff_a.size() == 0)
  {
    unsigned int n_carriers = this->n_known_carriers();

    _coeff_a.resize(n_carriers + 1, 0);
    _coeff_b.resize(n_carriers + 1, NEUMANN);
    _coeff_g.resize(n_carriers + 1, 0);

    _jacobian.resize(n_carriers + 1);
    for (unsigned int i = 0; i < n_carriers + 1; i++)
      _jacobian[i].resize(n_carriers + 1, 0);
  }

  _coeff_b[var] = type;
  if (type == DIRICHLET)
  {
    _coeff_a[var] = 1.0;
    _jacobian[var][var] = -1;
  }
  else if (type == NEUMANN)
    _coeff_a[var] = 0.0;
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

  std::map<ID, CarrierProperties*> carriers(_ddprop_A->get_carrier_properties());

  // carriers of other side, if different
  if ((_ddprop_B != nullptr) && (_ddprop_B != _ddprop_A))
  {
    for (auto&& it : _ddprop_B->get_carrier_properties())
    {
      if (carriers.find(it.first) == carriers.end())
        carriers[it.first] = it.second;
    }
  }

  set_carrier_properties(carriers);

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


    // add carriers which are present in B but not in A
    for (auto&& it : carriersB)
    {
      if (carriers.find(it.first) == carriers.end())
        carriers[it.first] = it.second;
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

  // surface states
  if (get_type(n_known_carriers()) == NEUMANN)
  {
    // NOTE we invert the signs because g = \epsilon \nabla\varphi \hat{n}
    // i.e. refers to the negative charge density
    // 2014-8-19 the above note seems misleading, the signs below are the right ones!
    //calculate_traps();
    double q = (pd.ionized_electron_traps + pd.ionized_hole_traps);
    if (is_internal_boundary())
      q /= 2;

    _coeff_g[n_known_carriers()] += q;

    for (unsigned int i = 0; i < n_known_carriers(); i++)
    {
      double dq_dEf = -pd.ionized_traps_derivative[i];
      if (is_internal_boundary())
        dq_dEf /= 2;

      _jacobian[n_known_carriers()][n_known_carriers()] += dq_dEf;
      _jacobian[n_known_carriers()][i] -= dq_dEf;
    }
  }

  if (get_number_of_recombination_models() > 0)
  {
    calculate_net_recombination_rates();

    // 2023-06-06, there is some uncertainty to whether a
    // factor of 1/2 should be used for internal boundaries.
    // We should check if it is actually counted two times in
    // assembly.
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
