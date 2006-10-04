#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "Dopant.h"
#include "Constants.h"

#include "elem.h"

#include <cmath>

// we calculate in cm, therefore the factor 1e6
// the electron charge enters because we take k*T in electron volts
const double
DriftDiffusionProperties::_DOS_factor = std::pow(2.0 * M_PI * Constants::me /
      (Constants::h * Constants::h) * Constants::e, 1.5) / 1e6;


DriftDiffusionProperties::DriftDiffusionProperties(void)
  : PhysicalProperties("DriftDiffusionProperties"),
    _elem(NULL),
    charge_density_derivatives(3, 0.0),
    electron_conductivity_derivatives(3, 0.0),
    hole_conductivity_derivatives(3, 0.0),
    electron_recombination_rate_derivatives(3, 0.0),
    hole_recombination_rate_derivatives(3, 0.0),
    _statistics(TiberCad::BOLTZMANN),
    _coupling(DriftDiffusionDefs::BOTH)
{
}

DriftDiffusionProperties::DriftDiffusionProperties(
    const DriftDiffusionProperties& rhs)
  : PhysicalProperties("DriftDiffusionProperties"),
    _elem(NULL),
    charge_density_derivatives(3, 0.0),
    electron_conductivity_derivatives(3, 0.0),
    hole_conductivity_derivatives(3, 0.0),
    electron_recombination_rate_derivatives(3, 0.0),
    hole_recombination_rate_derivatives(3, 0.0),
    _statistics(rhs._statistics),
    _coupling(rhs._coupling)
{
}



DriftDiffusionProperties::~DriftDiffusionProperties(void)
{
}


void
DriftDiffusionProperties::add_dopant(Dopant* dopant)
{
  if (dopant != NULL)
  {
    if (dopant->get_type() == Dopant::N_TYPE)
      _donors.insert(dopant);
    else
      _acceptors.insert(dopant);
  }
    
}

void
DriftDiffusionProperties::add_recombination_model(
    RecombinationModelInterface* recomb_model)
{
  if (recomb_model != NULL)
  {
    ID id = recomb_model->get_id();
    _recombination_models[id] = recomb_model;
    recomb_model->set_driftdiffusionproperties(this);
  }
}


double
DriftDiffusionProperties::get_total_donor_density(void) const
{
  double Nd = 0;
  dopant_iterator it = _donors.begin();
  dopant_iterator end = _donors.end();
  for ( ; it != end; ++it)
    Nd += (*it)->get_doping_density();

  return Nd;
}

double
DriftDiffusionProperties::get_total_acceptor_density(void) const
{
  double Na = 0;
  dopant_iterator it = _acceptors.begin();
  dopant_iterator end = _acceptors.end();
  for ( ; it != end; ++it)
    Na += (*it)->get_doping_density();

  return Na;
}

void
DriftDiffusionProperties::clear_doping(void)
{
  dopant_iterator it = _donors.begin();
  dopant_iterator end = _donors.end();
  for ( ; it != end; ++it) delete (*it);

  it = _acceptors.begin();
  end = _acceptors.end();
  for ( ; it != end; ++it) delete (*it);

  _donors.clear();
  _acceptors.clear();
}

void
DriftDiffusionProperties::clear_recombination(void)
{
  recomb_iterator it = _recombination_models.begin();
  recomb_iterator end = _recombination_models.end();
  for ( ; it != end; ++it) delete (it->second);

  _recombination_models.clear();
}


void
DriftDiffusionProperties::calculate_densities(double potential,
    double fermi_e, double fermi_h)
{

  // for now, all are equal
  double kT = electron_vt;
  double kTe = electron_vt;
  double kTh = hole_vt;
  
  // 1 conduction band
  const BandProperties& cb = conduction_band;
  const BandProperties& vb = valence_band;

  double Ec = get_conduction_band_edge();
  double Ev = get_valence_band_edge();
  
  // 1.) electron and hole density
  double n = 0, dn = 0, dn2 = 0, dn_over_n = 0, arg_e;
  double p = 0, dp = 0, dp2 = 0, dp_over_p = 0, arg_h;
  //if (_coupling & DriftDiffusionDefs::ELECTRONS)
  //{
    arg_e = -fermi_e + potential - Ec;
    if (get_statistics() == TiberCad::FERMIDIRAC)
    {
      density_and_derivatives<TiberCad::FERMIDIRAC>(arg_e / kTe,
          n, dn, dn2, dn_over_n);
    }
    else
    {
      density_and_derivatives<TiberCad::BOLTZMANN>(arg_e / kTe,
          n, dn, dn2, dn_over_n);
    }
  
    double Nc = cb.effective_DOS;
    n *= Nc;
    dn *= Nc / kTe;
    dn2 *= Nc / (kTe * kTe);
    dn_over_n /= kTe;

    electron_density = n;
    electron_density_derivative = dn;
  //}

  //if (_coupling & DriftDiffusionDefs::HOLES)
  //{
    arg_h = fermi_h - potential + Ev;

    if (get_statistics() == TiberCad::FERMIDIRAC)
    {
      density_and_derivatives<TiberCad::FERMIDIRAC>(arg_h / kTh,
          p, dp, dp2, dp_over_p);
    }
    else
    {
      density_and_derivatives<TiberCad::BOLTZMANN>(arg_h / kTh,
          p, dp, dp2, dp_over_p);
    }

    double Nv = vb.effective_DOS;
    p *= Nv;
    dp *= -Nv / kTh;
    dp2 *= Nv / (kTh * kTh);
    dp_over_p /= -kTh;

    hole_density = p;
    hole_density_derivative = dp;
  //}

}


void
DriftDiffusionProperties::calculate_all(double potential,
    double fermi_e, double fermi_h, const Point& coord)
{
  // remember the coordinates we are working on, models could need it
  _coord = &coord;

  // for now, all are equal
  double kT = electron_vt;
  double kTe = electron_vt;
  double kTh = hole_vt;
  
  // 1 conduction band
  const BandProperties& cb = conduction_band;
  const BandProperties& vb = valence_band;

  double Ec = get_conduction_band_edge();
  double Ev = get_valence_band_edge();
  
  // 1.) electron and hole density
  double n = 0, dn = 0, dn2 = 0, dn_over_n = 0, arg_e;
  double p = 0, dp = 0, dp2 = 0, dp_over_p = 0, arg_h;
  //if (_coupling & DriftDiffusionDefs::ELECTRONS)
  //{
    arg_e = -fermi_e + potential - Ec;
    if (get_statistics() == TiberCad::FERMIDIRAC)
    {
      density_and_derivatives<TiberCad::FERMIDIRAC>(arg_e / kTe,
          n, dn, dn2, dn_over_n);
    }
    else
    {
      density_and_derivatives<TiberCad::BOLTZMANN>(arg_e / kTe,
          n, dn, dn2, dn_over_n);
    }
  
    double Nc = cb.effective_DOS;
    n *= Nc;
    dn *= Nc / kTe;
    dn2 *= Nc / (kTe * kTe);
    dn_over_n /= kTe;

    electron_density = n;
    electron_density_derivative = dn;
  //}

  //if (_coupling & DriftDiffusionDefs::HOLES)
  //{
    arg_h = fermi_h - potential + Ev;

    if (get_statistics() == TiberCad::FERMIDIRAC)
    {
      density_and_derivatives<TiberCad::FERMIDIRAC>(arg_h / kTh,
          p, dp, dp2, dp_over_p);
    }
    else
    {
      density_and_derivatives<TiberCad::BOLTZMANN>(arg_h / kTh,
          p, dp, dp2, dp_over_p);
    }

    double Nv = vb.effective_DOS;
    p *= Nv;
    dp *= -Nv / kTh;
    dp2 *= Nv / (kTh * kTh);
    dp_over_p /= -kTh;

    hole_density = p;
    hole_density_derivative = dp;
  //}

  // 2.) ionized dopant densities
  double Nd = 0, dNd = 0;
  double Na = 0, dNa = 0;
  {
    dopant_iterator it = _donors.begin();
    dopant_iterator end = _donors.end();
    for ( ; it != end; ++it)
    {
      Nd += (*it)->get_ionized_dopant_density(arg_e, kT);
      dNd += (*it)->get_ionized_dopant_density_derivative(arg_e, kT);
    }
    ionized_donor_density = Nd;
    ionized_donor_density_derivative = dNd;

    it = _acceptors.begin();
    end = _acceptors.end();
    for ( ; it != end; ++it)
    {
      Na += (*it)->get_ionized_dopant_density(arg_h, kT);
      dNa -= (*it)->get_ionized_dopant_density_derivative(arg_h, kT);
    }
    ionized_acceptor_density = Na;
    ionized_acceptor_density_derivative = dNa;
  }


  // 3.) total charge density
  // NOTE: the sign change comes from the fact, that d/dphi = -d/dFn
  charge_density = p - n + Nd - Na;
  charge_density_derivatives[0] =  dp - dn + dNd - dNa;
  charge_density_derivatives[1] =       dn - dNd;
  charge_density_derivatives[2] = -dp            + dNa;

/*
  // 4.) mobilities / conductivities
  // For both statistics:
  // 
  //   mu_n = e * D_n * (1 / n) * (dn / dEf_e)
  //
  // NOTE: kT := kB * T / e includes already e
  if (_coupling & DriftDiffusionDefs::ELECTRONS)
  {
    double electron_diffusivity = kT * cb.mobility;
    electron_mobility = electron_diffusivity * dn_over_n;
    electron_conductivity = electron_diffusivity * dn;
    electron_conductivity_derivatives[0] = electron_diffusivity * dn2;
    electron_conductivity_derivatives[1] = electron_diffusivity * dn2;
  }
  if (_coupling & DriftDiffusionDefs::HOLES)
  {
    double hole_diffusivity = kT * vb.mobility;
    hole_mobility = -hole_diffusivity * dp_over_p;
    hole_conductivity = -hole_diffusivity * dp;
    hole_conductivity_derivatives[0] = -hole_diffusivity * dp2;
    hole_conductivity_derivatives[2] = -hole_diffusivity * dp2;
  }
*/
  
  {
    electron_recombination_rate = 0;
    electron_recombination_rate_derivatives[0] = 0;
    electron_recombination_rate_derivatives[1] = 0;
    electron_recombination_rate_derivatives[2] = 0;
    hole_recombination_rate = 0;
    hole_recombination_rate_derivatives[0] = 0;
    hole_recombination_rate_derivatives[1] = 0;
    hole_recombination_rate_derivatives[2] = 0;

    double Re, Rh;
    std::vector<double> dRe(3), dRh(3);

    recomb_iterator it = _recombination_models.begin();
    recomb_iterator end = _recombination_models.end();
    for ( ; it != end; ++it)
    {
      (it->second)->get_net_recombination_rates(Re, Rh);
      (it->second)->get_net_recombination_rate_derivatives(dRe, dRh);
      
      electron_recombination_rate += Re;
      electron_recombination_rate_derivatives[0] += dRe[0];
      electron_recombination_rate_derivatives[1] += dRe[1];
      electron_recombination_rate_derivatives[2] += dRe[2];
      hole_recombination_rate += Rh;
      hole_recombination_rate_derivatives[0] += dRh[0];
      hole_recombination_rate_derivatives[1] += dRh[1];
      hole_recombination_rate_derivatives[2] += dRh[2];
    }
  }
}

// TODO
void
DriftDiffusionProperties::get_net_recombination_rates(
    std::vector<double>& rates)
{
}

int
DriftDiffusionProperties::get_net_recombination_rate_IDs(
    std::vector<ID>& ids)
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

double
DriftDiffusionProperties::get_net_recombination_rate(ID id)
{
  double r = 0.0, dummy;
  
  RecombinationModelInterface* rec =
    get_recombination_model(id);
  if (rec != NULL)
    rec->get_net_recombination_rates(r, dummy);

  return r;
}




void
DriftDiffusionProperties::calculate_equilibrium_properties(int coupling,
    double temperature)
{
  
  // call this method to properly set conduction and valence band DOS
  // and energy
  setup_band_edges();

  // FIXME change this implementation, don't use PETSc

  // remember the coupling
  int coupling_bkp = _coupling;
  _coupling = DriftDiffusionDefs::BOTH;

  double kT = SimulationOptions::T * Constants::k_B;

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

  ierr = SNESSetType(snes, SNESLS);
  // set ls_maxstep to something reasonable
  ierr = SNESSetLineSearchParams(snes, PETSC_DEFAULT, 1 * kT,
      PETSC_DEFAULT);


  ierr = SNESSetFunction(snes, r, function, (void *) this);
  ierr = SNESSetJacobian(snes, J, J, jacobian, (void *) this);

  // calculate first guess according to doping
  // assume complete ionization
  const BandProperties& cb = conduction_band;
  const BandProperties& vb = valence_band;
  double Nd = get_total_donor_density();
  double Na = get_total_acceptor_density();

  double ni2 = cb.effective_DOS * vb.effective_DOS
    * std::exp(-get_band_gap() / kT);
  double ni = std::sqrt(ni2);
  equilibrium_electron_density = ni;
  equilibrium_hole_density = ni;
  
  // Hmm... Is there a better guess?
  if (Nd > Na)
  {
    guess = cb.band_edge - kT
      * std::log(cb.effective_DOS / (Nd + ni));
  }
  else
  {
    guess = vb.band_edge + kT
      * std::log(vb.effective_DOS / (Na + ni));
  }
  guess = guess;

  ierr = VecSet(&guess, x);

  ierr = SNESSolve(snes, x);

  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason(snes, &reason);
  int n_iterations = 0;
  ierr = SNESGetIterationNumber(snes, &n_iterations);
  if (reason < 0)
  {
    std::cerr << "ATTENTION Equilibrium properties calculation:\n";
    std::cerr << "  # iterations: "  << n_iterations
      << ", converged reason: " << reason << "\n";
  }

  ierr = VecGetArray(x, &result);

  equilibrium_fermi_level =  result[0];

  equilibrium_electron_density = electron_density;
  equilibrium_hole_density = hole_density;

  ierr = VecRestoreArray(x, &result);

  ierr = VecDestroy(x);
  ierr = MatDestroy(J);

  ierr = SNESDestroy(snes);

  
  // restore original coupling
  _coupling = coupling_bkp;

}


PetscErrorCode
DriftDiffusionProperties::jacobian(SNES snes, Vec x,
    Mat *jac, Mat *B, MatStructure *flag, void *sc)
{
  PetscScalar    *xx, A[1];
  PetscErrorCode ierr;
  PetscInt       idx[1] = {0};

  ierr = VecGetArray(x, &xx);

  DriftDiffusionProperties* s = static_cast<DriftDiffusionProperties*>(sc);
  s->calculate_all(xx[0], 0.0, 0.0, (s->get_element())->centroid());

  A[0] = s->get_charge_density_derivatives()[0];
  //std::cerr << xx[0] << " " << s->get_charge_density() << 
  //  " " << A[0] << "\n";

  ierr = MatSetValues(*jac, 1, idx, 1, idx, A, INSERT_VALUES);
  *flag = SAME_NONZERO_PATTERN;

  ierr = VecRestoreArray(x, &xx);

  ierr = MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);

  return ierr;
}


PetscErrorCode
DriftDiffusionProperties::function(SNES snes, Vec x, Vec f, void *sc)
{
  PetscErrorCode ierr;
  PetscScalar    *xx, *ff;

  ierr = VecGetArray(x, &xx);
  ierr = VecGetArray(f, &ff);

  DriftDiffusionProperties* s = static_cast<DriftDiffusionProperties*>(sc);
  s->calculate_all(xx[0], 0.0, 0.0, (s->get_element())->centroid());

  ff[0] = s->get_charge_density();
  //std::cerr << xx[0] << " " << ff[0] << "\n";

  ierr = VecRestoreArray(x, &xx);
  ierr = VecRestoreArray(f, &ff);

  return ierr;
}


