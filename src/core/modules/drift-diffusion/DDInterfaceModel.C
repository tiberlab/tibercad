// $Id$

#include "DDInterfaceModel.h"
#include "Material.h"
#include "MaterialBoundary.h"
#include "Trap.h"
#include "FowlerNordheim.h"
#include "RecombinationModelInterface.h"
#include "DriftDiffusionProperties.h"
#include "SimulationInterface.h"
#include "ModelErrorException.h"
#include "Variable.h"
#include "elem.h"

using namespace std;


DDInterfaceModel::DDInterfaceModel(const ModelOptions& options) :
  PhysicalModel(options),
  _internal_bd(false),
  _has_current(false),
  _emission(NULL),
  _eflux(0.0),
  _eflux_sim(NULL),
  _eflux_controlled(false)
{
  _coeff_a.resize(3, 0);
  _coeff_b.resize(3, NEUMANN);
  _coeff_g.resize(3, 0);

  _jacobian.resize(3);
  for (unsigned int i = 0; i < 3; i++)
    _jacobian[i].resize(3, 0);

}


DriftDiffusionProperties*
DDInterfaceModel::get_dd_properties(void) const
{
  DriftDiffusionProperties* ddprop = NULL;

  if (get_material() != NULL)
    ddprop = static_cast<DriftDiffusionProperties*>(
        get_material()->get_model(get_simulator_id()));

  return ddprop;
}


void
DDInterfaceModel::prepare_submodels(void)
{
  // for each trap we add an SRH recombination model
  ModelOptions::submodel_iterator it(get_options().submodels_begin("trap"));
  ModelOptions::submodel_iterator end(get_options().submodels_end("trap"));
  for (; it != end; ++it)
  {
    ModelOptions opts(it->second);
    if (opts.get_option("recombination_center", true))
    {
      it->second.delete_option("recombination_center");
      opts.set_option("trap", true);
      opts.set_option("type", "srh");
      opts.set_option("name", "recombination");
      get_options().add_submodel("recombination", opts);
    }
  }

  vector<PhysicalModelInterface*> pd;
  create_submodels(pd, "recombination");

  // traps
  create_submodels(pd, "trap");

}



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

  const Material* mat = bnd->get_material_A();
  SimulationInterface* si = SimulationInterface::get_simulation(get_simulator_id());
  if (!si->includes_region(bnd->get_id_A()))
    mat = bnd->get_material_B();

  assert(mat != NULL);

  // we set a bulk material, just in case a submodel needs it
  set_material(mat);


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

  if (!_recombination_models.empty())
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
        _eflux_sim = SimulationInterface::find_simulation(eflux);
        if (_eflux_sim == NULL)
        {
          string msg("DriftDiffusion Interface: ");
          msg += "no electron current simulation '" + eflux + "' found.";
          throw InitFailedException(msg);
        }

        _eflux_id = _eflux_sim->get_solution_id("eCurrentDensity");
        if (_eflux_id == INVALID_ID)
          throw InitFailedException("Module '" +
              eflux + "' does not contain solution variable 'eCurrentDensity'");
      }
    }

    _eflux_controlled = true;
    set_type(0, NEUMANN);
    set_type(1, NEUMANN);
    set_type(2, NEUMANN);
  }

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

  do_compute();

  // now add common stuff if needed

  // surface states
  if (get_type(0) == NEUMANN)
  {
    double q, dq_dEfp, dq_dEfn;
    _calculate_traps(q, dq_dEfn, dq_dEfp);
    // NOTE we invert the signs because g = \epsilon \nabla\varphi \hat{n}
    // i.e. refers to the negative charge density
    _coeff_g[0] += q;
    _jacobian[0][0] += dq_dEfn + dq_dEfp;
    _jacobian[0][1] -= dq_dEfn;
    _jacobian[0][2] -= dq_dEfp;
  }

  if (!_recombination_models.empty())
  {
    const DriftDiffusionProperties::PointData& pd =
        get_dd_properties()->get_point_data();

    double rec[6];
    _calculate_recombination(rec);

    _coeff_g[1] += rec[0];
    double dRn_dEfn = -rec[1] * pd.electron_density_derivative;
    double dRn_dEfp = -rec[2] * pd.hole_density_derivative;
    _jacobian[1][0] -= dRn_dEfn + dRn_dEfp;
    _jacobian[1][1] += dRn_dEfn;
    _jacobian[1][2] += dRn_dEfp;

    _coeff_g[2] += rec[3];
    double dRp_dEfn = -rec[4] * pd.electron_density_derivative;
    double dRp_dEfp = -rec[5] * pd.hole_density_derivative;
    _jacobian[2][0] -= dRp_dEfn + dRp_dEfp;
    _jacobian[2][1] += dRp_dEfn;
    _jacobian[2][2] += dRp_dEfp;
  }

  if (_eflux_controlled)
  {

    double flux = _eflux;
    if (this->is_internal_boundary())
      flux *= 0.5;

    if (_eflux_sim != NULL)
    {
      DriftDiffusionProperties& dd = *get_dd_properties();
      vector<double> data;
      // we take the flux from the neighbor element
      _eflux_sim->get_solution(dd.get_element()->neighbor(_side), _eflux_id,
          data, dd.get_coordinates());

      flux = data[0] * _normal(0) + data[1] * _normal(1) + data[2] * _normal(2);
    }


    _coeff_g[1] -= flux / Constants::e;

  }

}




void
DDInterfaceModel::_calculate_traps(double& q, double& dq_dEfn, double& dq_dEfp)
{
  DriftDiffusionProperties* ddprop = get_dd_properties();
  assert(ddprop != NULL);

  q = dq_dEfn = dq_dEfp = 0.0;

  double Ec = ddprop->get_conduction_band_edge() - ddprop->get_electric_potential();
  double Ev = ddprop->get_valence_band_edge() - ddprop->get_electric_potential();

  const DriftDiffusionProperties::PointData& pd = ddprop->get_point_data();

  double ionized_electron_traps = 0.0;
  double ionized_electron_traps_derivative = 0.0;
  if (_etraps.size() > 0)
  {
    double nt = 0, dnt = 0;
    double kT = pd.electron_vt;
    set<Trap*>::iterator it(_etraps.begin());
    const set<Trap*>::iterator end(_etraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, -pd.fermi_e, kT);
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
    double kT = pd.hole_vt;
    set<Trap*>::iterator it(_htraps.begin());
    const set<Trap*>::iterator end(_htraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, -pd.fermi_h, kT);
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
