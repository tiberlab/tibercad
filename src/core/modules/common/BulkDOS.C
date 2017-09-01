// $Id$

#include "BulkDOS.h"
#include "Constants.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include "InitFailedException.h"
#include "SimulationOptions.h"
#include "Messages.h"

#include "TiberModule.h"


using namespace std;

BulkDOS::BulkDOS(const ModelOptions& options) :
  DensityOfStates(options)
{
}



////////////
/// NOTE ///
////////////
//
// Here we apply first VCA and the calculate the temperature corrected
// band gap. I'm not sure if it should be done the other way round
//


void
BulkDOS::read_database(void)
{
  // when reading from the database, we use the same data
  // as for kp
  const Database& db = get_database();

  if (get_particle() == 'e')
  {
    // TODO: should bowing be applied to Eg(T) or Eg(0) ?

    db.set_section("valenceband");
    double Ev = db.get("E_v", 0.0);

    db.set_section("bandgap");
    double bandgap = db.get("Eg_G", 1e3);
    double varshni_a = db.get("varshni_alpha_G", 0.0);
    double varshni_b = db.get("varshni_beta_G", 0.0);

    double Eg_L = db.get("Eg_L", 1e3);
    double Eg_X = db.get("Eg_X", 1e3);

    db.set_section("conductionband");
    double mass = db.get("m_G", 1.0);
    // spin degeneracy
    unsigned short int deg = 2;

    //effective_mass() = mass;

    if ((Eg_L < bandgap) || (get_option("valley", "G") == "L"))
    {
      bandgap = Eg_L;

      db.set_section("bandgap");
      varshni_a = db.get("varshni_alpha_L", 0.0);
      varshni_b = db.get("varshni_beta_L", 0.0);

      db.set_section("conductionband");
      double m_L_t = db.get("m_L_t", 1.0);
      double m_L_l = db.get("m_L_l", 1.0);
      // 2 = spin degeneracy, 8 = valley degeneracy
      deg = 2 * 8;
      mass = std::pow(m_L_t * m_L_t * m_L_l, 1.0/3.0 );

      //effective_mass() = mass * std::pow(64.0, 1.0/3.0);
    }
    if ((Eg_X < bandgap) || (get_option("valley", "G") == "X"))
    {
      bandgap = Eg_X;

      db.set_section("bandgap");
      varshni_a = db.get("varshni_alpha_X", 0.0);
      varshni_b = db.get("varshni_beta_X", 0.0);

      db.set_section("conductionband");
      double m_X_t = db.get("m_X_t", 1.0);
      double m_X_l = db.get("m_X_l", 1.0);
      // 2 = spin degeneracy, 6 = valley degeneracy
      deg = 2 * 6;
      mass = std::pow(m_X_t * m_X_t * m_X_l, 1.0/3.0 );

      //effective_mass() = mass * std::pow(36.0, 1.0/3.0);
    }

    varshni_a = get_option("varshni_a", varshni_a);
    varshni_b = get_option("varshni_b", varshni_b);

    _dos_mass.resize(1, mass);
    _degeneracy.resize(1, deg);
    GapParameters gapparm;
    gapparm.Ev = Ev;
    gapparm.Eg0 = bandgap;
    gapparm.varshni_a = varshni_a;
    gapparm.varshni_b = varshni_b;
    _gap_params.resize(1, gapparm);

    effective_mass()[0] = mass;
  }
  else if (get_particle() == 'h')
  {
    db.set_section("valenceband");
    reference_energy()[0] = db.get("E_v", 0.0);

    // we put the degeneracy to two, since the DOS mass contains the band degeneracy
    unsigned short int deg = db.get("degeneracy", 2);
    double mass = db.get("m_dos", 1.0);

    _dos_mass.resize(1, mass);
    _degeneracy.resize(1, deg);

    effective_mass()[0] = mass;
  }
}


void
BulkDOS::do_init(void)
{

  _ref_energies.resize(0);
  get_parameter("level", reference_energy()[0]);
  _ref_energies.push_back(reference_energy()[0]);

  get_option("levels", _ref_energies);
  if (_ref_energies.size() > 1)
    throw InitFailedException("bulk DOS model can not yet handle multiple bands");

  // for the electrons, reference energy is the CB edge, but the database values
  // are referred to the VB edge _Ev
  if (get_particle() == 'e')
  {
    if (has_parameter("level") || has_option("levels"))
      for (int i = 0; i < _gap_params.size(); ++i)
        _gap_params[i].zero();
    else
    {
      _ref_energies.resize(_gap_params.size(), 0.0);
      for (int i = 0; i < _gap_params.size(); ++i)
        reference_energy()[0] = _gap_params[i].gap(SimulationOptions::T) + _gap_params[i].Ev;
    }
  }
  else
    reference_energy()[0] = _ref_energies[0];


  get_option("dos_mass", _dos_mass);


  get_option("degeneracy", _degeneracy);


  if (_dos_mass.size() != _ref_energies.size())
    throw InitFailedException("Bulk DOS model for material " + get_material()->get_name() +
        ": number of energy levels different from number of DOS masses");

  if (_degeneracy.size() != _ref_energies.size())
    throw InitFailedException("Bulk DOS model for material " + get_material()->get_name() +
        ": number of energy levels different from number of degeneracies");

  _dos_factor = pow(2.0 * M_PI *
        Constants::me / (Constants::h * Constants::h) *
        Constants::e, 1.5) / 1e6;

  effective_dos() = _dos_factor * _degeneracy[0] * std::pow(
      SimulationOptions::T * Constants::k_B * _dos_mass[0], 1.5);

}


void
BulkDOS::do_print_info(void)
{
  ostringstream os;
  os << "band : energy  DOS mass  degeneracy\n";
  for (int i = 0; i < _ref_energies.size(); ++i)
  {
    double erg = _ref_energies[i];

    if (get_particle() == 'e')
      erg += _gap_params[i].Ev + _gap_params[i].gap(SimulationOptions::T);

    os << (i + 1) << " : " << erg << "  " << _dos_mass[i]
       << "  " << _degeneracy[i] << "\n";
  }
  Messages::info(os.str());
}



std::pair<double, double>
BulkDOS::calculate_density_and_derivative(double Ef, double Epot,
    double kT, double kTlattice, const Elem* elem, const Point& p) const
{
  return calculate_density_and_derivative(Ef, Epot, kT, kTlattice);
}

std::pair<double, double>
BulkDOS::calculate_density_and_derivative(double Ef, double Epot, double kT, double kTlattice) const
{
  double density = 0.0;
  double derivative = 0.0;

  const double arg_max = 150;
  const double arg_min = -50;
  const double min_dens = 1e-64;

  _th_el_power = 0;

  for (int i = 0; i < _ref_energies.size(); ++i)
  {
    double dens, der;

    double energy = _ref_energies[i];
   

    if (get_particle() == 'e')
    {
      double T = kTlattice / Constants::k_B;
      // _ref_energy[i] is 0.0 if _gap_params are used, and vice versa
      energy += _gap_params[i].Ev + _gap_params[i].gap(T);
    }
    else if (get_particle() == 'h')
      energy = -energy;

    double arg = (Ef - energy - Epot) / kT;

    double Neff = _dos_factor * _degeneracy[i] * std::pow(kT * _dos_mass[i], 1.5);

    if (arg < arg_min)
    {
      dens = exp(arg);
      der = dens;
    }
    else if (arg < arg_max)
    {
      dens = TiberMath::fermidirac_half(arg);
      der = TiberMath::fermidirac_mhalf(arg);
    }
    else
    {
      dens = 2.0 * M_2_SQRTPI / 3.0 * std::pow(arg, 1.5);
      der = M_2_SQRTPI * std::sqrt(arg);
    }

    dens *= Neff;

    // is this needed???
    if (dens > min_dens)
      der *= Neff / kT;
    else
      der = 0.0;

    density += dens;
    derivative += der;

    // calculate thermoelectric power
    double temp = kT / Constants::k_B;
    double Pth = -dens / der * 1.5 / temp + Constants::k_Boltzmann * (arg - 1);
    if (get_particle() == 'h')
    {
      Pth *= -1;
    }
    _th_el_power += dens * Pth;

  }

  _th_el_power /= density;


  return make_pair(density, derivative);
}
