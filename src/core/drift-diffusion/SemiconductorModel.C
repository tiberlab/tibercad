// $Id$

#include "Alloy.h"

#include "SemiconductorModel.h"
#include "ZbDDsemiconductor.h"
#include "WzDDsemiconductor.h"

#include "point.h"
#include "elem.h"
#include "getpot.h"

#include <iostream>

using namespace DriftDiffusionDefs;
using namespace std;


SemiconductorModel::~SemiconductorModel(void)
{
  delete _bulk_model;
}

SemiconductorModel::SemiconductorModel(void)
  : DriftDiffusionProperties(),
    _recombination(0),
    _coupling(BOTH),
    _bulk_model(NULL)
{
  // spin degeneracy will be included in the degeneracy from
  // DDsemiconductor
  _DOS_factor = std::pow(2 * M_PI * Constants::me /
      (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;
}

SemiconductorModel::SemiconductorModel(
    const SemiconductorModel& model)
  : DriftDiffusionProperties(model),
    _recombination(model._recombination),
    _DOS_factor(model._DOS_factor),
    _conduction_band(model._conduction_band),
    _valence_band(model._valence_band),
    _coupling(model._coupling),
    _electron_recombination_time(model._electron_recombination_time),
    _hole_recombination_time(model._hole_recombination_time),
    _direct_rec_param(model._direct_rec_param),
    _filename(model._filename)
{
}

void
SemiconductorModel::read_database(const Dummy&)
{
  GetPot data(_filename);

  const std::string structure = data("structure", "zb");

  if (structure == "zb")
  {
    // defaults for GaAs
    ZbDDsemiconductor::ZbDDparameters params;

    params.EgGamma = data("Eg_G", 1.519);
    params.EgL = data("Eg_L", 1.815);
    params.EgX = data("Eg_X", 1.981);
    params.Ev = data("E_v", 1.346);

    params.m_G = data("m_G", 0.067);
    params.m_l_L = data("m_L_l", 1.9);
    params.m_t_L = data("m_L_t", 0.0754);
    params.m_l_X = data("m_X_l", 1.3);
    params.m_t_X = data("m_X_t", 0.23);

    params.a_c = data("a_c", -9.36);
    params.a_v = data("a_v", -1.21);
    params.b = data("b", -2.0);
    params.d = data("d", -4.8);

    params.delta = data("delta", 0.341);
    params.gamma1 = data("gamma1", 6.98);
    params.gamma2 = data("gamma2", 2.06);
    params.gamma3 = data("gamma3", 2.93);

    params.def_vol_X = data("abs_def_pot_X", -0.16);
    params.def_uniax_X = data("uniax_def_pot_X", 14.26);
    params.def_vol_L = data("abs_def_pot_L", -4.91);
    params.def_uniax_L = data("uniax_def_pot_L", 6.5);

    ZbDDsemiconductor* zbsc = new ZbDDsemiconductor(params);
    _bulk_model = zbsc;

    permittivity = data("permittivity", 12.93);
    _conduction_band.mobility = data("electron_mobility", 1000.0);
    _valence_band.mobility = data("hole_mobility", 200.0);

  }
  else
  {
    // defaults for GaN
    WzDDsemiconductor::WzDDparameters params;

    params.EgGamma = data("Eg_G", 3.51);
    params.Ev = data("E_v", -0.726);

    params.m_c_zz = data("m_c_zz", 0.20);
    params.m_c_xx = data("m_c_xx", 0.20);

    params.A1 = data("A1", -7.21);
    params.A2 = data("A2", -0.44);
    params.A3 = data("A3", 6.68);
    params.A4 = data("A4", -3.46);
    params.A5 = data("A5", -3.40);
    params.A6 = data("A6", -4.90); 

    params.a_x = data("a_x", -4.9);
    params.a_z = data("a_z", -11.3);

    params.D1 = data("D1", -3.7);
    params.D2 = data("D2", 4.5);
    params.D3 = data("D3", 8.2);
    params.D4 = data("D4", -4.1);
    params.D5 = data("D5", -4.0);
    params.D6 = data("D6", -5.5);
    params.delta_s = data("delta_s", 0.017);
    params.delta_cr = data("delta_cr", 0.010);

    WzDDsemiconductor* wzsc = new WzDDsemiconductor(params);
    _bulk_model = wzsc;
    wzsc->energy_cutoff = 2.0;

    permittivity = data("permittivity", 9.5);
    _conduction_band.mobility = data("electron_mobility", 1000.0);
    _valence_band.mobility = data("hole_mobility", 200.0);
  }
}

void
SemiconductorModel::build_alloy(const std::string& component2,
    const std::string& bowing_params, double content)
{
  GetPot data(component2);
  GetPot bowing(bowing_params);

  const std::string structure = data("structure", "zb");

  double (*alloy)(double, double, double, double) =
    Alloy::calculate_VCA_parameter;

  if (structure == "zb")
  {
    ZbDDsemiconductor* sc = dynamic_cast<ZbDDsemiconductor*>(_bulk_model);
    ZbDDsemiconductor::ZbDDparameters& params = sc->get_parameters();

    params.EgGamma = alloy(data("Eg_G", 1.519), params.EgGamma, content,
       bowing("Eg_G", 0.0));
    //cerr << "EgGamma = " << params.EgGamma << endl;
    params.EgL = alloy(data("Eg_L", 1.815), params.EgL, content,
       bowing("Eg_L", 0.0));
    //cerr << "EgL = " << params.EgL << endl;
    params.EgX = alloy(data("Eg_X", 1.981), params.EgX, content,
       bowing("Eg_X", 0.0));
    //cerr << "EgX = " << params.EgX << endl;
    params.Ev = alloy(data("E_v", 1.346), params.Ev, content,
       bowing("E_v", 0.0));
    //cerr << "Ev = " << params.Ev << endl;

    params.m_G = alloy(data("m_G", 0.067), params.m_G, content,
       bowing("m_G", 0.0));
    //cerr << "m_G = " << params.m_G << endl;
    params.m_t_L = alloy(data("m_L_t", 0.0754), params.m_t_L, content,
       bowing("m_L_t", 0.0));
    //cerr << "m_t_L = " << params.m_t_L << endl;
    params.m_l_L = alloy(data("m_L_l", 1.9), params.m_l_L, content,
       bowing("m_L_l", 0.0));
    //cerr << "m_l_L = " << params.m_l_L << endl;
    params.m_t_X = alloy(data("m_X_t", 1.3), params.m_t_X, content,
       bowing("m_X_t", 0.0));
    //cerr << "m_t_X = " << params.m_t_X << endl;
    params.m_l_X = alloy(data("m_X_l", 0.23), params.m_l_X, content,
       bowing("m_X_l", 0.0));
    //cerr << "m_l_X = " << params.m_l_X << endl;

    params.a_c = alloy(data("a_c", -9.36), params.a_c, content,
       bowing("a_c", 0.0));
    //cerr << "a_c = " << params.a_c << endl;
    params.a_v = alloy(data("a_v", -1.21), params.a_v, content,
       bowing("a_v", 0.0));
    //cerr << "a_v = " << params.a_v << endl;
    params.b = alloy(data("b", -2.0), params.b, content,
       bowing("b", 0.0));
    //cerr << "b = " << params.b << endl;
    params.d = alloy(data("d", -4.8), params.d, content,
       bowing("d", 0.0));
    //cerr << "d = " << params.d << endl;

    params.delta = alloy(data("delta", 0.341), params.delta, content,
       bowing("delta", 0.0));
    //cerr << "delta = " << params.delta << endl;
    params.gamma1 = alloy(data("gamma1", 6.98), params.gamma1, content,
       bowing("gamma1", 0.0));
    //cerr << "gamma1 = " << params.gamma1 << endl;
    params.gamma2 = alloy(data("gamma2", 2.06), params.gamma2, content,
       bowing("gamma2", 0.0));
    //cerr << "gamma2 = " << params.gamma2 << endl;
    params.gamma3 = alloy(data("gamma3", 2.93), params.gamma3, content,
       bowing("gamma3", 0.0));
    //cerr << "gamma3 = " << params.gamma3 << endl;

    params.def_vol_X = alloy(data("abs_def_pot_X", -0.16),
        params.def_vol_X, content, bowing("abs_def_pot_X", 0.0));
    //cerr << "def_vol_X = " << params.def_vol_X << endl;
    params.def_uniax_X = alloy(data("uniax_def_pot_X", 14.26),
        params.def_uniax_X, content, bowing("uniax_def_pot_X", 0.0));
    //cerr << "def_uniax_X = " << params.def_uniax_X << endl;
    params.def_vol_L = alloy(data("abs_def_pot_L", -4.91),
        params.def_vol_L, content, bowing("abs_def_pot_L", 0.0));
    //cerr << "def_vol_L = " << params.def_vol_L << endl;
    params.def_uniax_L = alloy(data("uniax_def_pot_L", 6.5),
        params.def_uniax_L, content, bowing("uniax_def_pot_L", 0.0));
    //cerr << "def_uniax_L = " << params.def_uniax_L << endl;

    permittivity = alloy(data("permittivity", 12.93), permittivity, content,
       bowing("permittivity", 0.0));
    _conduction_band.mobility = alloy(data("electron_mobility", 1000.0),
        _conduction_band.mobility, content, bowing("electron_mobility", 0.0));
    _valence_band.mobility = alloy(data("hole_mobility", 200.0),
        _valence_band.mobility, content, bowing("hole_mobility", 0.0));
  }
  else
  {
    WzDDsemiconductor* sc = dynamic_cast<WzDDsemiconductor*>(_bulk_model);
    WzDDsemiconductor::WzDDparameters& params = sc->get_parameters();

    params.EgGamma = alloy(data("Eg_G", 3.51), params.EgGamma, content,
        bowing("Eg_G", 0.0));
    params.Ev = alloy(data("E_v", -0.726), params.Ev, content,
        bowing("E_v", 0.0));

    params.m_c_zz = alloy(data("m_c_zz", 0.20), params.m_c_zz, content,
        bowing("m_c_zz", 0.0));
    params.m_c_xx = alloy(data("m_c_xx", 0.20), params.m_c_xx, content,
        bowing("m_c_xx", 0.0));

    params.A1 = alloy(data("A1", -7.21), params.A1, content,
        bowing("A1", 0.0));
    params.A2 = alloy(data("A2", -0.44), params.A2, content,
        bowing("A2", 0.0));
    params.A3 = alloy(data("A3", 6.68), params.A3, content,
        bowing("A3", 0.0));
    params.A4 = alloy(data("A4", -3.46), params.A4,content, 
        bowing("A4", 0.0));
    params.A5 = alloy(data("A5", -3.40), params.A5, content,
        bowing("A5", 0.0));
    params.A6 = alloy(data("A6", -4.90), params.A6,content, 
        bowing("A6", 0.0));

    params.a_x = alloy(data("a_x", -4.9), params.a_x, content,
        bowing("a_x", 0.0));
    params.a_z = alloy(data("a_z", -11.3), params.a_z,content, 
        bowing("a_z", 0.0));

    params.D1 = alloy(data("D1", -3.7), params.D1, content,
        bowing("D1", 0.0));
    params.D2 = alloy(data("D2", 4.5), params.D2,content, 
        bowing("D2", 0.0));
    params.D3 = alloy(data("D3", 8.2), params.D3, content,
        bowing("D3", 0.0));
    params.D4 = alloy(data("D4", -4.1), params.D4, content,
        bowing("D4", 0.0));
    params.D5 = alloy(data("D5", -4.0), params.D5, content,
        bowing("D5", 0.0));
    params.D6 = alloy(data("D6", -5.5), params.D6, content,
        bowing("D6", 0.0));
    params.delta_s = alloy(data("delta_s", 0.017), params.delta_s,
        content, bowing("delta_s", 0.0));
    params.delta_cr = alloy(data("delta_cr", 0.010), params.delta_cr,
        content, bowing("delta_cr", 0.0));

    permittivity = alloy(data("permittivity", 9.5), permittivity,
        content, bowing("permittivity", 0.0));
    _conduction_band.mobility = alloy(data("electron_mobility", 1000.0),
      _conduction_band.mobility, content, bowing("electron_mobility", 0.0));
    _valence_band.mobility = alloy(data("hole_mobility", 200.0),
      _valence_band.mobility, content, bowing("hole_mobility", 0.0));
  }
}

void
SemiconductorModel::prepare_element_data(void)
{
  double kT = SimulationOptions::T * Constants::k_B;
  electron_vt = hole_vt = kT;
  
  BandProperties& cb = _conduction_band;
  BandProperties& vb = _valence_band;

  cb.effective_DOS =
    get_DOS_factor() * std::pow(kT * cb.effective_mass, 1.5);

  vb.effective_DOS =
    get_DOS_factor() * std::pow(kT * vb.effective_mass, 1.5);

  conduction_band_edge = cb.band_edge;
  valence_band_edge = vb.band_edge;
}

void
SemiconductorModel::calculate_all(
    double potential, double fermi_e, double fermi_h,
    const Point& coord)
{
  switch (_coupling & BOTH)
  {
    case ELECTRONS:
      calculate_all<ELECTRONS>(potential, fermi_e, fermi_h);
      break;
    case HOLES:
      calculate_all<HOLES>(potential, fermi_e, fermi_h);
      break;
    default:
      calculate_all<BOTH>(potential, fermi_e, fermi_h);
      break;
  }
}

//
// TODO
// Very crude implementation at the moment
// 
void
SemiconductorModel::extract_band_properties(void)
{
  // treat conduction band
  const std::vector<DDsemiconductor::band_extremum>& cbs =
    _bulk_model->get_conduction_band_energy_mass();
  // get minimum
  int id = 0;
  for (int i = 1 ; i < cbs.size(); i++)
  {
    if (cbs[i].energy < cbs[id].energy)
      id = i;
  }
  _conduction_band.band_edge = cbs[id].energy;
  _conduction_band.effective_mass = cbs[id].mass_DOS
      * std::pow(cbs[id].degeneracy, 2.0 / 3.0);

  /*
  cerr << "CB:\n";
  for (int i = 0 ; i < cbs.size(); i++)
  {
    cerr << " Ec = " << cbs[i].energy
      << ", m = " << cbs[i].mass_DOS
      << ", d = " << cbs[i].degeneracy << endl;
  }
  */

  
  // treat valence band
  const std::vector<DDsemiconductor::band_extremum>& vbs =
    _bulk_model->get_valence_band_energy_mass();
  // get maximum
  id = 0;
  int id2 = -1;
  double delta_min = 5e-3;
  for (int i = 1 ; i < vbs.size(); i++)
  {
    // if the bands are (nearly) degenerate we have to take both into
    // account
    if (fabs((vbs[id].energy + delta_min) - vbs[i].energy) < 2 * delta_min)
    {
      id2 = id;
      id = i;
      break;
    }
    if (vbs[i].energy > (vbs[id].energy + delta_min))
      id = i;
  }
  _valence_band.band_edge = vbs[id].energy;
  if (id2 >= 0)
  {
    // two degenerate bands
    double tmp = vbs[id].degeneracy * std::pow(vbs[id].mass_DOS, 1.5)
      + vbs[id2].degeneracy * std::pow(vbs[id2].mass_DOS, 1.5);
    _valence_band.effective_mass = std::pow(tmp, 2.0 / 3.0);
  }
  else
    _valence_band.effective_mass = vbs[id].mass_DOS
      * std::pow(vbs[id].degeneracy, 2.0 / 3.0);
  
  /*
  cerr << "VB:\n";
  for (int i = 0 ; i < vbs.size(); i++)
  {
    cerr << " Ev = " << vbs[i].energy
      << ", m = " << vbs[i].mass_DOS
      << ", d = " << vbs[i].degeneracy << endl;
  }
  */
}


void
SemiconductorModel::calculate_equilibrium_properties(int coupling,
    double temperature)
{

  assert(_bulk_model != NULL);

  // calculate conduction and valence band data
  _bulk_model->calculate_conduction_band_extremum();
  _bulk_model->calculate_valence_band_extremum();

  // get the band properties from _bulk_model
  extract_band_properties();

  
  // call this method to properly set conduction and valence band DOS
  // and energy
  SemiconductorModel::prepare_element_data();


  // remember the coupling
  int coupling_bkp = _coupling;
  _coupling = BOTH;

  SNES           snes;
  KSP            ksp;
  PC             pc;
  Vec            x, r;
  Mat            J;
  PetscErrorCode ierr;
  PetscScalar    guess, *result;

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
  ierr = VecCreateSeq(PETSC_COMM_SELF, 1, &x);
  ierr = VecDuplicate(x, &r);
  ierr = MatCreate(PETSC_COMM_SELF,
      PETSC_DECIDE, PETSC_DECIDE, 1, 1, &J);
  ierr = MatSetFromOptions(J);

  ierr = SNESGetKSP(snes, &ksp);
  ierr = KSPGetPC(ksp, &pc);
  ierr = PCSetType(pc, PCLU);
  ierr = KSPSetTolerances(ksp, 1.e-9, 1e-12, PETSC_DEFAULT, 1500);
  ierr = SNESSetTolerances(snes, PETSC_DEFAULT,
      1.e-12, PETSC_DEFAULT, 50, 5000);

  double kT = SimulationOptions::T * Constants::k_B;

  ierr = SNESSetType(snes, SNESLS);
  // set ls_maxstep to something reasonable
  ierr = SNESSetLineSearchParams(snes, PETSC_DEFAULT, 100 * kT,
      PETSC_DEFAULT);


  ierr = SNESSetFunction(snes, r, function, (void *) this);
  ierr = SNESSetJacobian(snes, J, J, jacobian, (void *) this);

  // calculate first guess according to doping
  // assume complete ionization
  const BandProperties& cb = _conduction_band;
  const BandProperties& vb = _valence_band;
  double Nd = get_donor_density();
  double Na = get_acceptor_density();

  double ni2 = cb.effective_DOS * vb.effective_DOS
    * std::exp(-0.5 * get_band_gap() / kT);
  double ni = std::sqrt(ni2);
  
  // Hmm... Is there a better guess?
  if (Nd > Na)
  {
    guess = cb.band_edge + kT
      * std::log(cb.effective_DOS / (Nd + ni));
  }
  else
  {
    guess = vb.band_edge + kT
      * std::log(vb.effective_DOS / (Na + ni));
  }

  ierr = VecSet(&guess, x);

  ierr = SNESSolve(snes, x);

  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason(snes, &reason);
  int n_iterations = 0;
  ierr = SNESGetIterationNumber(snes, &n_iterations);
  if (reason < 0)
  {
    cerr << "ATTENTION Equilibrium properties calculation:\n";
    cerr << "  # iterations: "  << n_iterations
      << ", converged reason: " << reason << "\n";
  }

  ierr = VecGetArray(x, &result);

  equilibrium_fermi_level =  result[0];

  equilibrium_electron_density = electron_density;
  equilibrium_hole_density = hole_density;

  ierr = VecRestoreArray(x, &result);

  ierr = VecDestroy(x);
  ierr = MatDestroy(J);

  // restore original coupling
  _coupling = coupling_bkp;
}


PetscErrorCode
SemiconductorModel::jacobian(SNES snes, Vec x,
    Mat *jac, Mat *B, MatStructure *flag, void *sc)
{
  PetscScalar    *xx, A[1];
  PetscErrorCode ierr;
  PetscInt       idx[1] = {0};

  ierr = VecGetArray(x, &xx);

  SemiconductorModel* s = static_cast<SemiconductorModel*>(sc);
  s->calculate_all(xx[0], 0.0, 0.0, (s->elem)->centroid());

  A[0] = s->get_charge_density_derivatives()[0];
/*
  std::cerr << "u = " << xx[0] << ", n = " << s->get_charge_density()
    << ", dn = " << A[0] << "\n";
*/
  ierr = MatSetValues(*jac, 1, idx, 1, idx, A, INSERT_VALUES);
  *flag = SAME_NONZERO_PATTERN;

  ierr = VecRestoreArray(x, &xx);

  ierr = MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);

  return ierr;
}

PetscErrorCode
SemiconductorModel::function(SNES snes, Vec x, Vec f, void *sc)
{
  PetscErrorCode ierr;
  PetscScalar    *xx, *ff;

  ierr = VecGetArray(x, &xx);
  ierr = VecGetArray(f, &ff);

  SemiconductorModel* s = static_cast<SemiconductorModel*>(sc);
  s->calculate_all(xx[0], 0.0, 0.0, (s->elem)->centroid());

  ff[0] = s->get_charge_density();

  ierr = VecRestoreArray(x, &xx);
  ierr = VecRestoreArray(f, &ff);

  return ierr;
}


void
SemiconductorModel::print_info(void) const
{
  cout << " - conduction band:\n";
  cout << "    Ec = " << _conduction_band.band_edge
    << ", Nc = " << _conduction_band.effective_DOS
    << ", n0 = " << get_equilibrium_electron_density()
    << ", m_DOS = " << _conduction_band.effective_mass;
  cout << endl;
  cout << " - valence band:\n";
  cout << "    Ev = " << _valence_band.band_edge
    << ", Nv = " << _valence_band.effective_DOS
    << ", p0 = " << get_equilibrium_hole_density()
    << ", m_DOS = " << _valence_band.effective_mass;
  cout << endl;
  cout << " - Ef0 = " << get_equilibrium_fermi_level()
    << ", ni^2 = " << get_intrinsic_density_squared();
  cout << endl;
}
