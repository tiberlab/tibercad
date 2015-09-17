// $Id$

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
  _coeff_a.resize(3, 0);
  _coeff_b.resize(3, NEUMANN);
  _coeff_g.resize(3, 0);

  _jacobian.resize(3);
  for (unsigned int i = 0; i < 3; i++)
    _jacobian[i].resize(3, 0);

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

  // Set bulk conduction and valence band since submodels may need them
  set_conduction_band(&_ddprop_A->get_conduction_band());
  set_valence_band(&_ddprop_A->get_valence_band());

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
    set_type(1, NEUMANN);
    set_type(2, NEUMANN);
  }



  if (get_option("field_emission", "") == "fowler_nordheim")
  {
    _emission = new FowlerNordheim(get_option("work_function", 4.7));
    _emission->set_velocity(get_option("emission_velocity", 1.5e8));
  }

  get_option("area_factor", "");

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
    set_type(0, NEUMANN);
    set_type(1, NEUMANN);
    set_type(2, NEUMANN);
  }

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
  set_conduction_band(&_ddprop_A->get_conduction_band());
  set_valence_band(&_ddprop_A->get_valence_band());

  // this one may be NULL
  if (_ddprop_B != NULL)
  {
    _ddprop_B->reinit(elem);
    double cb_A = _ddprop_A->get_conduction_band_edge();
    double vb_A = _ddprop_A->get_valence_band_edge();

    if (_ddprop_B->get_conduction_band_edge() < cb_A)
      set_conduction_band(&_ddprop_B->get_conduction_band());

    if (_ddprop_B->get_valence_band_edge() > vb_A)
      set_valence_band(&_ddprop_B->get_valence_band());
  }

  // we can get the lattice temperature from the bulk model
  double kT = _ddprop_A->get_lattice_temperature();
  set_lattice_vt(kT);
  set_carrier_temperatures(kT, kT);

}


void
DDInterfaceModel::compute()
{
  for (unsigned int i = 0; i < 3; i++)
  {
    if (get_type(i) == NEUMANN)
    {
      _coeff_g[i] = 0;
      _jacobian[i][0] = 0;
      _jacobian[i][1] = 0;
      _jacobian[i][2] = 0;
    }
  }
  calculate_densities();
  calculate_traps();

  do_compute();

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

  if (get_number_of_recombination_models() > 0)
  {
    calculate_net_recombination_rates();

    int div = is_internal_boundary() ? 2 : 1;

    _coeff_g[1] += pd.electron_recombination_rate / div;
    double dRn_dEfn = -pd.electron_recombination_rate_derivatives[0] *
        pd.electron_density_derivative / div;
    double dRn_dEfp = -pd.electron_recombination_rate_derivatives[1] *
        pd.hole_density_derivative / div;
    _jacobian[1][0] -= dRn_dEfn + dRn_dEfp;
    _jacobian[1][1] += dRn_dEfn;
    _jacobian[1][2] += dRn_dEfp;

    _coeff_g[2] += pd.hole_recombination_rate / div;
    double dRp_dEfn = -pd.hole_recombination_rate_derivatives[0] *
        pd.electron_density_derivative / div;
    double dRp_dEfp = -pd.hole_recombination_rate_derivatives[1] *
        pd.hole_density_derivative / div;
    _jacobian[2][0] -= dRp_dEfn + dRp_dEfp;
    _jacobian[2][1] += dRp_dEfn;
    _jacobian[2][2] += dRp_dEfp;
  }
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

}

