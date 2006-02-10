// $Id: SemiconductorModel.C 14 2006-01-27 11:44:23Z maufder $

#include "SemiconductorModel.h"
#include "Constants.h"
#include "iostream"

SemiconductorModel::SemiconductorModel(void)
  : _statistics(SimulationOptions::BOLTZMANN),
    _recombination(0),
    _thermal_voltage(SimulationOptions::T * Constants::k_B),
    polarization(3, 0.0)
{
    _DOS_factor = 2 * std::pow(2 * M_PI * Constants::me /
        (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;
}

SemiconductorModel::SemiconductorModel(const SemiconductorModel& model)
  : _statistics(model._statistics),
    _recombination(model._recombination),
    _thermal_voltage(model._thermal_voltage),
    _DOS_factor(model._DOS_factor),
    _equilibrium_electron_density(model._equilibrium_electron_density),
    _equilibrium_hole_density(model._equilibrium_hole_density),
    _equilibrium_fermi_level(model._equilibrium_fermi_level),
    _conduction_band(model._conduction_band),
    _valence_band(model._valence_band),
    _material(model._material),
    polarization(model.polarization)
{
}


void
SemiconductorModel::calculate_equilibrium_properties(double temperature)
{
  _thermal_voltage = Constants::k_B * temperature;

  BandProperties& cb = conduction_band_properties();
  BandProperties& vb = valence_band_properties();

  cb.band_edge = _material.conduction_band_edge;
  vb.band_edge = _material.valence_band_edge;

  cb.low_field_mobility = _material.electron_mobility;
  vb.low_field_mobility = _material.hole_mobility;

  cb.effective_mass = _material.electron_effective_mass;
  vb.effective_mass = _material.hole_effective_mass;

  cb.effective_DOS =
    get_DOS_factor() * std::pow(_thermal_voltage * cb.effective_mass, 1.5);

  vb.effective_DOS =
    get_DOS_factor() * std::pow(_thermal_voltage * vb.effective_mass, 1.5);
  
  _equilibrium_fermi_level = calculate_equilibrium_potential();

  CalculatedProperties prop;
  double arg_e = (_equilibrium_fermi_level -
      cb.band_edge) / _thermal_voltage;
  double arg_h = -(_equilibrium_fermi_level -
      vb.band_edge) / _thermal_voltage;
  calculate_e_h_densities(arg_e, arg_h, _thermal_voltage, prop);
  _equilibrium_electron_density = prop.electron_density;
  _equilibrium_hole_density = prop.hole_density;
}
    
void
SemiconductorModel::calculate_all(double potential, double Ef_e, double Ef_h,
    CalculatedProperties& result)
{

  double kT = get_thermal_voltage();
  
  const BandProperties& cb = get_conduction_band_properties();
  const BandProperties& vb = get_valence_band_properties();

  double Ec = cb.band_edge;
  double Ev = vb.band_edge;
  
  // 1.) electron and hole density
  double arg_e = (Ef_e + potential - Ec) / kT;
  double arg_h = -(Ef_h + potential - Ev) / kT;
  calculate_e_h_densities(arg_e, arg_h, kT, result);

  // 2.) ionized dopant density
  calculate_ionized_dopants(arg_e, arg_h, kT, result);

  // 3.) charge density
  calculate_charge_density(result);

  // just for nicer code
  double& n  = result.electron_density;
  double& p  = result.hole_density;
  double& dn  = result.electron_density_derivative;
  double& dp  = result.hole_density_derivative;

  // 4.) mobilities
  // For both statistics:
  // 
  //   mu_n = e * D_n * (1 / n) * (dn / dEf_e)
  //
  // NOTE: kT := kB * T / e includes already e
  result.electron_diffusivity = kT * cb.low_field_mobility;
  result.hole_diffusivity = kT * vb.low_field_mobility;
  result.electron_mobility = result.electron_diffusivity * result.dn_over_n;
  result.hole_mobility = -result.hole_diffusivity * result.dp_over_p;

  // 5.) Recombination
  result.net_recombination_rate = 0;
  result.net_recombination_rate_derivatives[0] = 0;
  result.net_recombination_rate_derivatives[1] = 0;
  result.net_recombination_rate_derivatives[2] = 0;
  if (_recombination && SRH)
    calculate_SRH_recombination(result);
  if (_recombination && AUGER)
    calculate_Auger_recombination(result);
  if (_recombination && DIRECT)
    calculate_direct_recombination(result);

  //TODO this is just for testing
  result.polarization[0] = polarization[0];
  result.polarization[1] = polarization[1];
  result.polarization[2] = polarization[2];
 
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
    ierr = SNESSetLineSearchParams(snes, PETSC_DEFAULT, 100 * _thermal_voltage,
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
      * std::exp(-0.5 * get_band_gap() / get_thermal_voltage());
    double ni = std::sqrt(ni2);
    // Hmm... Is there a better guess?
    if (Nd > Na)
    {
      guess = cb.band_edge + _thermal_voltage
        * std::log(cb.effective_DOS / (Nd + ni));
    }
    else
    {
      guess = vb.band_edge + _thermal_voltage
        * std::log(vb.effective_DOS / (Na + ni));
    }
    ierr = VecSet(&guess, x);

    ierr = SNESSolve(snes, x);
/*
    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason(snes, &reason);
    int n_iterations = 0;
    ierr = SNESGetIterationNumber(snes, &n_iterations);
    std::cerr << "iterations: "  << n_iterations
      << "  converged reason: " << reason << "\n";
*/

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

    SemiconductorModel::CalculatedProperties properties;
    SemiconductorModel* s = static_cast<SemiconductorModel*>(sc);
    s->calculate_all(xx[0], 0.0, 0.0, properties);

    A[0] = properties.charge_density_derivatives[0];
/*
    std::cerr << "u = " << xx[0] << ", n = " << properties.charge_density
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

    SemiconductorModel::CalculatedProperties properties;
    SemiconductorModel* s = static_cast<SemiconductorModel*>(sc);
    s->calculate_all(xx[0], 0.0, 0.0, properties);

    ff[0] = properties.charge_density;

    ierr = VecRestoreArray(x, &xx);
    ierr = VecRestoreArray(f, &ff);

    return ierr;
}
