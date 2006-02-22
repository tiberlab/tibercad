// $Id$

#include "SemiconductorModel.h"
#include "Constants.h"
#include "DriftDiffusionDefs.h"

#include "point.h"

#include <iostream>

using namespace DriftDiffusionDefs;

SemiconductorModel::SemiconductorModel(void)
  : DriftDiffusionProperties(),
    _recombination(0),
    polarization(3, 0.0)
{
    _DOS_factor = 2 * std::pow(2 * M_PI * Constants::me /
        (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;
}

SemiconductorModel::SemiconductorModel(const SemiconductorModel& model)
  : DriftDiffusionProperties(model),
    _recombination(model._recombination),
    _DOS_factor(model._DOS_factor),
    _conduction_band(model._conduction_band),
    _valence_band(model._valence_band),
    _material(model._material),
    polarization(model.polarization)
{
}


void
SemiconductorModel::calculate_equilibrium_properties(double temperature)
{
  double thermal_voltage = Constants::k_B * temperature;

  BandProperties& cb = conduction_band_properties();
  BandProperties& vb = valence_band_properties();

  cb.band_edge = _material.conduction_band_edge;
  vb.band_edge = _material.valence_band_edge;

  cb.low_field_mobility = _material.electron_mobility;
  vb.low_field_mobility = _material.hole_mobility;

  cb.effective_mass = _material.electron_effective_mass;
  vb.effective_mass = _material.hole_effective_mass;

  cb.effective_DOS =
    get_DOS_factor() * std::pow(thermal_voltage * cb.effective_mass, 1.5);

  vb.effective_DOS =
    get_DOS_factor() * std::pow(thermal_voltage * vb.effective_mass, 1.5);
  
  equilibrium_fermi_level = calculate_equilibrium_potential();

  // they were calculated during 'calculate_equilibrium_potential()'
  equilibrium_electron_density = electron_density;
  equilibrium_hole_density = hole_density;
}
    
void
SemiconductorModel::calculate_all(double potential, double fermi_e,
    double fermi_h,
    const Point& coord, const Elem* elem, int coupling)
{

  double kT = Constants::k_B * SimulationOptions::T;
  
  const BandProperties& cb = get_conduction_band_properties();
  const BandProperties& vb = get_valence_band_properties();

  double Ec = cb.band_edge;
  double Ev = vb.band_edge;
  
  // 1.) electron and hole density
  double n, p, dn, dp, dn2, dp2, dn_over_n, dp_over_p;
  double arg_e = (fermi_e + potential - Ec) / kT;
  double arg_h = -(fermi_h + potential - Ev) / kT;
  if (get_statistics() == TiberCad::FERMIDIRAC)
  {
    density_and_derivatives<TiberCad::FERMIDIRAC>(arg_e,
        n, dn, dn2, dn_over_n);

    density_and_derivatives<TiberCad::FERMIDIRAC>(arg_h,
        p, dp, dp2, dp_over_p);
  }
  else
  {
    density_and_derivatives<TiberCad::BOLTZMANN>(arg_e,
        n, dn, dn2, dn_over_n);

    density_and_derivatives<TiberCad::BOLTZMANN>(arg_h,
        p, dp, dp2, dp_over_p);
  }
  
  double Nc = _conduction_band.effective_DOS;
  n *= Nc;
  dn *= Nc / kT;
  dn2 *= Nc / (kT * kT);
  dn_over_n /= kT;

  double Nv = _valence_band.effective_DOS;
  p *= Nv;
  dp *= -Nv / kT;
  dp2 *= Nv / (kT * kT);
  dp_over_p /= -kT;

  electron_density = n;
  electron_density_derivative = dn;
  hole_density = p;
  hole_density_derivative = dp;

  // 2.) ionized dopant densities
  double Nd, Na, dNd, dNa;
  calculate_ionized_donors(arg_e, kT, Nd, dNd);
  calculate_ionized_acceptors(arg_h, kT, Na, dNa);
  ionized_donor_density = Nd;
  ionized_acceptor_density = Na;

  // 3.) total charge density
  charge_density = p - n + Nd - Na;
  charge_density_derivatives[0] = dp - dn + dNd - dNa;
  charge_density_derivatives[1] =    - dn + dNd;
  charge_density_derivatives[2] = dp            - dNa;
  
  // 4.) mobilities
  // For both statistics:
  // 
  //   mu_n = e * D_n * (1 / n) * (dn / dEf_e)
  //
  // NOTE: kT := kB * T / e includes already e
  double electron_diffusivity = kT * cb.low_field_mobility;
  double hole_diffusivity = kT * vb.low_field_mobility;
  electron_mobility = electron_diffusivity * dn_over_n;
  hole_mobility = -hole_diffusivity * dp_over_p;
  
  // 5.) conductivities
  electron_conductivity = electron_diffusivity * dn;
  hole_conductivity = -hole_diffusivity * dp;
  electron_conductivity_derivatives[0] = electron_diffusivity * dn2;
  electron_conductivity_derivatives[1] = electron_diffusivity * dn2;
  hole_conductivity_derivatives[0] = -hole_diffusivity * dp2;
  hole_conductivity_derivatives[2] = -hole_diffusivity * dp2;

  // 6.) Recombination
  electron_recombination_rate = 0;
  electron_recombination_rate_derivatives[0] = 0;
  electron_recombination_rate_derivatives[1] = 0;
  electron_recombination_rate_derivatives[2] = 0;
  hole_recombination_rate = 0;
  hole_recombination_rate_derivatives[0] = 0;
  hole_recombination_rate_derivatives[1] = 0;
  hole_recombination_rate_derivatives[2] = 0;
  if (_recombination & SRH)
    calculate_SRH_recombination();
  if (_recombination & AUGER)
    calculate_Auger_recombination();
  if (_recombination & DIRECT)
    calculate_direct_recombination();


  DriftDiffusionProperties::polarization(0) = polarization[0];
  DriftDiffusionProperties::polarization(1) = polarization[1];
  DriftDiffusionProperties::polarization(2) = polarization[2];
}


double
SemiconductorModel::calculate_equilibrium_potential(void) const
{
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
  ierr = KSPSetTolerances(ksp, 1.e-4, 1e-12, PETSC_DEFAULT, 500);
  ierr = SNESSetTolerances(snes, PETSC_DEFAULT,
      1.e-12, PETSC_DEFAULT, 50, 500);

  ierr = SNESSetType(snes, SNESLS);
  // set ls_maxstep to something reasonable
  double kT = Constants::k_B * SimulationOptions::T;
  ierr = SNESSetLineSearchParams(snes, PETSC_DEFAULT, 100 * kT,
      PETSC_DEFAULT);


  ierr = SNESSetFunction(snes, r, function, (void *) this);
  ierr = SNESSetJacobian(snes, J, J, jacobian, (void *) this);

  // calculate first guess according to doping
  // assume complete ionization
  const BandProperties& cb = get_conduction_band_properties();
  const BandProperties& vb = get_valence_band_properties();
  double Nd = get_donor_density();
  double Na = get_acceptor_density();

  double ni2 = cb.effective_DOS * vb.effective_DOS
    * std::exp(-0.5 * get_band_gap() / kT);
  double ni = std::sqrt(ni2);
  // Hmm... Is there a better guess?
  if (Nd > Na)
  {
    guess = cb.band_edge + kT * std::log(cb.effective_DOS / (Nd + ni));
  }
  else
  {
    guess = vb.band_edge + kT * std::log(vb.effective_DOS / (Na + ni));
  }
  ierr = VecSet(&guess, x);

  ierr = SNESSolve(snes, x);
///*
  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason(snes, &reason);
  int n_iterations = 0;
  ierr = SNESGetIterationNumber(snes, &n_iterations);
  std::cerr << "iterations: "  << n_iterations
    << "  converged reason: " << reason << "\n";
//*/

  ierr = VecGetArray(x, &result);
  double eq_pot = result[0];
  ierr = VecRestoreArray(x, &result);

  ierr = VecDestroy(x);
  ierr = MatDestroy(J);

  return eq_pot;
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
  Point p;
  s->calculate_all(xx[0], 0.0, 0.0, p, NULL, 0);

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
SemiconductorModel::function(SNES snes, Vec x, Vec f, void *sc)
{
  PetscErrorCode ierr;
  PetscScalar    *xx, *ff;

  ierr = VecGetArray(x, &xx);
  ierr = VecGetArray(f, &ff);

  SemiconductorModel* s = static_cast<SemiconductorModel*>(sc);
  Point p;
  s->calculate_all(xx[0], 0.0, 0.0, p, NULL, 0);

  ff[0] = s->get_charge_density();

  ierr = VecRestoreArray(x, &xx);
  ierr = VecRestoreArray(f, &ff);

  return ierr;
}
