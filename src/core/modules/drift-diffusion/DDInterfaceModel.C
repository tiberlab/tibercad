// $Id$

#include "DDInterfaceModel.h"
#include "Material.h"
#include "MaterialBoundary.h"
#include "Trap.h"
#include "DriftDiffusionProperties.h"
#include "SimulationInterface.h"
#include "ModelErrorException.h"

using namespace std;


DDInterfaceModel::DDInterfaceModel(const ModelOptions& options) :
  PhysicalModel(options),
  _internal_bd(false),
  _has_current(false),
  _ddprop(NULL)
{
  _coeff_a.resize(3, 0);
  _coeff_b.resize(3, NEUMANN);
  _coeff_g.resize(3, 0);

  _jacobian.resize(3);
  for (unsigned int i = 0; i < 3; i++)
    _jacobian[i].resize(3, 0);

}



DDInterfaceModel*
DDInterfaceModel::create(const ModelOptions& options)
{
  DDInterfaceModel* model = NULL;

  string name("interface");
  name = options.get_option("type", name);

  if (name == "interface")
    model = dynamic_cast<DDInterfaceModel*>(
        PhysicalModelInterface::create(_create, _destroy, options));
  else
    model = dynamic_cast<DDInterfaceModel*>(
        PhysicalModelInterface::create("ddbnd_" + name, options));

  if (model == NULL)
    throw ModelErrorException("Unknown drift-diffusion "
        "interface model: " + name);

  return model;
}


void
DDInterfaceModel::do_init(void)
{
  MaterialBoundary* bnd = dynamic_cast<MaterialBoundary*>(get_owner());
  if (bnd == NULL)
    throw ModelErrorException("DriftDiffusion boundary models can "
        "be used only on region boundaries");

  Material* mat = bnd->get_material_A();
  SimulationInterface* si = SimulationInterface::get_simulation(get_simulator_id());
  if (!si->includes_region(bnd->get_id_A()))
    mat = bnd->get_material_B();

  assert(mat != NULL);

  _ddprop = static_cast<DriftDiffusionProperties*>(
      mat->get_model(get_simulator_id()));

  assert(_ddprop != NULL);

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
    cerr << "surface recombination: " << it->second->get_type() << endl;
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
}


void
DDInterfaceModel::_calculate_traps(double& q, double& dq_dEfn, double& dq_dEfp)
{

  double Ec = _ddprop->get_conduction_band_edge() - _ddprop->get_electric_potential();
  double Ev = _ddprop->get_valence_band_edge() - _ddprop->get_electric_potential();

  const DriftDiffusionProperties::PointData& pd = _ddprop->get_point_data();

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

