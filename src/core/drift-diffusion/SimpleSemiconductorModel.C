// $Id$

#include "SimpleSemiconductorModel.h"

#include "point.h"
#include "elem.h"

#include <iostream>

using namespace DriftDiffusionDefs;


SimpleSemiconductorModel::SimpleSemiconductorModel(void)
  : DriftDiffusionProperties(),
    _recombination(0),
    _is_prepared(false),
    _coupling(BOTH)
{
  _DOS_factor = 2 * std::pow(2 * M_PI * Constants::me /
      (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;
}

SimpleSemiconductorModel::SimpleSemiconductorModel(
    const SimpleSemiconductorModel& model)
  : DriftDiffusionProperties(model),
    _recombination(model._recombination),
    _DOS_factor(model._DOS_factor),
    _conduction_band(model._conduction_band),
    _valence_band(model._valence_band),
    _is_prepared(model._is_prepared),
    _coupling(model._coupling),
    _electron_recombination_time(model._electron_recombination_time),
    _hole_recombination_time(model._hole_recombination_time),
    _direct_rec_param(model._direct_rec_param)
{
}

void
SimpleSemiconductorModel::prepare_element_data(void)
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
SimpleSemiconductorModel::calculate_all(
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


void
SimpleSemiconductorModel::calculate_equilibrium_properties(int coupling,
    double temperature)
{

  // if equilibrium properties were already calculated, we
  // just update all properties to equilibrium ones
  if (_is_prepared)
  {
    switch (_coupling & BOTH)
    {
      case ELECTRONS:
        calculate_all<ELECTRONS>(get_equilibrium_fermi_level(), 0, 0);
        break;
      case HOLES:
        calculate_all<HOLES>(get_equilibrium_fermi_level(), 0, 0);
        break;
      default:
        calculate_all<BOTH>(get_equilibrium_fermi_level(), 0, 0);
        break;
    }
    return;
  }

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
  ierr = PCSetType(pc, PCNONE);
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
///*
  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason(snes, &reason);
  int n_iterations = 0;
  ierr = SNESGetIterationNumber(snes, &n_iterations);
  std::cerr << "Equilibrium properties calculation:\n";
  std::cerr << "  # iterations: "  << n_iterations
    << ", converged reason: " << reason << "\n";
//*/

  ierr = VecGetArray(x, &result);

  equilibrium_fermi_level =  result[0];

  equilibrium_electron_density = electron_density;
  equilibrium_hole_density = hole_density;

  ierr = VecRestoreArray(x, &result);

  ierr = VecDestroy(x);
  ierr = MatDestroy(J);

  _is_prepared = true;

  // restore original coupling
  _coupling = coupling_bkp;
}


PetscErrorCode
SimpleSemiconductorModel::jacobian(SNES snes, Vec x,
    Mat *jac, Mat *B, MatStructure *flag, void *sc)
{
  PetscScalar    *xx, A[1];
  PetscErrorCode ierr;
  PetscInt       idx[1] = {0};

  ierr = VecGetArray(x, &xx);

  SimpleSemiconductorModel* s = static_cast<SimpleSemiconductorModel*>(sc);
  s->calculate_all(xx[0], 0.0, 0.0, (s->elem)->centroid());

  A[0] = s->get_charge_density_derivatives()[0];
///*
  std::cerr << "u = " << xx[0] << ", n = " << s->get_charge_density()
    << ", dn = " << A[0] << "\n";
//*/
  ierr = MatSetValues(*jac, 1, idx, 1, idx, A, INSERT_VALUES);
  *flag = SAME_NONZERO_PATTERN;

  ierr = VecRestoreArray(x, &xx);

  ierr = MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);

  return ierr;
}

PetscErrorCode
SimpleSemiconductorModel::function(SNES snes, Vec x, Vec f, void *sc)
{
  PetscErrorCode ierr;
  PetscScalar    *xx, *ff;

  ierr = VecGetArray(x, &xx);
  ierr = VecGetArray(f, &ff);

  SimpleSemiconductorModel* s = static_cast<SimpleSemiconductorModel*>(sc);
  s->calculate_all(xx[0], 0.0, 0.0, (s->elem)->centroid());

  ff[0] = s->get_charge_density();

  ierr = VecRestoreArray(x, &xx);
  ierr = VecRestoreArray(f, &ff);

  return ierr;
}


