// $Id$

#ifndef _SEMICONDUCTOR_H_
#define _SEMICONDUCTOR_H_

#include "SimulationOptions.h"
#include "Dopant.h"

// GNU scientific library
#include <gsl/gsl_sf_fermi_dirac.h>

#include <vector>
#include <cmath>

extern "C" {
  #include "petscsnes.h"
}

// forward declarations
//class Point;
class Elem;

class SemiconductorModel
{
  public:

    /**
     * This structure holds the properties that are calculated by
     * \p calculate_all()
     */
    struct CalculatedProperties
    {
      CalculatedProperties(void);

      double ionized_donor_density;
      double ionized_donor_density_derivative;
      double ionized_acceptor_density;
      double ionized_acceptor_density_derivative;

      double electron_density;
      double electron_density_derivative;
      double electron_density_2nd_derivative;
      double electron_mobility;
      double electron_diffusivity;
      double dn_over_n;

      double hole_density;
      double hole_density_derivative;
      double hole_density_2nd_derivative;
      double hole_mobility;
      double hole_diffusivity;
      double dp_over_p;
      
      double electron_conductivity;
      double electron_conductivity_derivative;
      double hole_conductivity;
      double hole_conductivity_derivative;

      double charge_density;
      std::vector<double> charge_density_derivatives;

      double net_electron_recombination_rate;
      std::vector<double> net_electron_recombination_rate_derivatives;
      double net_hole_recombination_rate;
      std::vector<double> net_hole_recombination_rate_derivatives;

      std::vector<double> polarization;
    };

    /**
     * This contains all reference material properties for this model
     */
    struct MaterialDescriptor
    {
      Dopant n_dopant;
      Dopant p_dopant;

      double relative_permittivity;

      double valence_band_edge;
      double conduction_band_edge;

      double electron_effective_mass;
      double hole_effective_mass;

      double electron_mobility;
      double hole_mobility;

      double electron_recombination_time;
      double hole_recombination_time;
    };

    /**
     * This structure holds the basic properties of a band for given
     * conditions (temp etc.)
     */
    struct BandProperties
    {
      double effective_mass;
      double effective_DOS;
      double band_edge;
      double low_field_mobility;
    };


    /**
     * The known recombination models
     */
    enum RecombinationModels { SRH = 1, AUGER = 2, DIRECT = 4 };

    SemiconductorModel(void);
    SemiconductorModel(const SemiconductorModel& model);
    virtual ~SemiconductorModel(void) {};

    /**
     * Set the material for which the properties should be calculated.
     */
    void set_material_descriptor(const MaterialDescriptor& desc);


    /**
     * Set the statistics to be used.
     * 
     * The statistics can be one of \p SimulationOptions::BOLTZMANN or
     * \p SimulationOptions::FERMIDIRAC
     */
    void set_statistics(SimulationOptions::Statistics statistics);

    /**
     * Set conduction band properties
     */
    void set_conduction_band_properties(double band_edge,
        double effective_mass, double mobility);

    /**
     * Set valence band properties
     */
    void set_valence_band_properties(double band_edge,
        double effective_mass, double mobility);
    
    /**
     * Set n-type doping
     */
    void set_n_dopant(const Dopant& dopant);

    /**
     * Set p-type doping
     */
    void set_p_dopant(const Dopant& dopant);

    /**
     * Set the relative permittivity
     */
    void set_relative_permittivity(double epsilon_r);
    


    /**
     * Add recombination model to be used.
     */
    void add_recombination_model(RecombinationModels recomb_model);

    /**
     * Remove a recombination model.
     */
    void remove_recombination_model(RecombinationModels recomb_model);

    /**
     * Set Shockley-Read-Hall recombination parameters
     */
    void set_SRH_parameters(double tau_n, double tau_p);


    /**
     * @returns the statistics used for this model
     */
    SimulationOptions::Statistics get_statistics(void) const;


    /**
     * Calculates the equilibrium properties.
     *
     * This method has to be called before any call to \p calculate_all()
     * but after setting up all material parameters
     */
    virtual void calculate_equilibrium_properties(double temperature =
        SimulationOptions::T);


    /**
     * Assume strained semiconductor.
     *
     * This method has to be called before any call to \p calculate_all()
     * It has to be implemented if strain has to be modeled.
     */
    //virtual void calculate_strained_properties(const Elem* elem) {};
    virtual void calculate_strained_properties(Elem* elem) {};

    /**
     * Calulate all densities, recombination rates, mobilities
     *
     */
    virtual void calculate_all(double potential, double Ef_e, double Ef_h,
        CalculatedProperties& result);
    
    /**
     * @returns the thermal voltage
     */
    double get_thermal_voltage(void) const;

    /**
     * @returns the relative dielectric constant
     */
    double get_relative_permittivity(void) const;

    /**
     * Get the conduction band properties
     */
    const BandProperties& get_conduction_band_properties(void) const;

    /**
     * Get the valence band properties
     */
    const BandProperties& get_valence_band_properties(void) const;

    /**
     * @returns the donor density
     */
    double get_donor_density(void) const;

    /**
     * @returns the acceptor density
     */
    double get_acceptor_density(void) const;

    /**
     * @returns the intrinsic density n_i
     */
    double get_intrinsic_density(void) const;

    /**
     * @returns the square of the intrinsic density n_i
     */
    double get_intrinsic_density_squared(void) const;

    /**
     * @returns the band gap
     */
    double get_band_gap(void) const;

    /**
     * @returns the equilibrium Fermi level
     */
    double get_equilibrium_fermi_level(void) const;

    /**
     * @returns the equilibrium electron density
     */
    double get_equilibrium_electron_density(void) const;

    /**
     * @returns the equilibrium hole density
     */
    double get_equilibrium_hole_density(void) const;


    /**
     * Calculates and returns the equilibrium electric potential with
     * respect to a quasi fermi level of 0 eV
     */
    double calculate_equilibrium_potential(void) const;

    
    /**
     * @returns the material currently associated with this
     * \p SemiconductorModel object.
     */
    const MaterialDescriptor& get_material_descriptor(void) const;


    /**
     * The polarization
     */
    std::vector<double> polarization;

  protected:

    /**
     * @returns a writeable reference to the conduction band properties
     */
    BandProperties& conduction_band_properties(void);

    /**
     * @returns a writeable reference to the valence band properties
     */
    BandProperties& valence_band_properties(void);

    /**
     * @returns the factor 2 * pow(2 * PI / h^2)^1.5
     * for calculating the effective density of states
     */
    double get_DOS_factor(void) const;
    
    /**
     * Calculate electron and hole densities and derivatives
     *
     * The arguments are:
     *   arg_e = (Ef_e - Ec + potential) / kT
     *   arg_h = (Ev - Ef_h - potential) / kT
     */
    void calculate_e_h_densities(double arg_e, double arg_h, double kT,
        CalculatedProperties& result);

    /**
     * Calculate ionized donor and acceptor densities and derivatives
     *
     * The arguments are:
     *   arg_e = (Ef_e - Ec + potential) / kT
     *   arg_h = (Ev - Ef_h - potential) / kT
     */
    void calculate_ionized_dopants(double arg_e, double arg_h, double kT,
        CalculatedProperties& result);

    /**
     * Calculate total charge densitiy and derivatives
     */
    void calculate_charge_density(CalculatedProperties& result);

    /**
     * Calculate Shockley-Read-Hall recombination
     */
    void calculate_SRH_recombination(CalculatedProperties& result);

    /**
     * Calculate Auger recombination
     */
    void calculate_Auger_recombination(CalculatedProperties& result);

    /**
     * Calculate direct recombination
     */
    void calculate_direct_recombination(CalculatedProperties& result);

  private:

    /**
     * The statistics used for this model
     */
    SimulationOptions::Statistics _statistics;

    /**
     * The recombination models used
     */
    int _recombination;

    /**
     * The thermal voltage Vt = kB * T
     */
    double _thermal_voltage;

    /**
     * The factor 2 * pow(2 * PI / h^2)^1.5
     * for calculating the effective density of states
     */
    double _DOS_factor;

    /**
     * The equilibrium densities
     */
    double _equilibrium_electron_density;
    double _equilibrium_hole_density;

    /**
     * The equilibrium Fermi level
     */
    double _equilibrium_fermi_level;

    /**
     * The band properties
     */
    BandProperties _conduction_band;
    BandProperties _valence_band;

    /**
     * The basic (reference) material properties
     */
    MaterialDescriptor _material;


    /**
     * Calculate the density and its derivatives
     * 
     * The method also returns the ratio 'first derivative over density'
     * which is used for the mobility.
     */
    template<SimulationOptions::Statistics S>
    void density_and_derivatives(double arg, double& density,
        double& derivative, double& _2nd_derivative,
        double& derivative_over_density) const;
    

    static PetscErrorCode jacobian(SNES snes, Vec x,
        Mat *jac, Mat *B, MatStructure *flag, void *sc);

    static PetscErrorCode function(SNES snes, Vec x, Vec f, void *sc);
};


//
// inline member functions
//

inline
SemiconductorModel::CalculatedProperties::CalculatedProperties(void)
  : charge_density_derivatives(3),
    net_electron_recombination_rate_derivatives(3, 0.0),
    net_hole_recombination_rate_derivatives(3, 0.0),
    polarization(3, 0.0)
{
}

inline
void
SemiconductorModel::set_material_descriptor(const MaterialDescriptor& desc)
{
  _material = desc;
}

inline
void
SemiconductorModel::set_statistics(SimulationOptions::Statistics statistics)
{
  _statistics = statistics;
}

inline
SimulationOptions::Statistics
SemiconductorModel::get_statistics(void) const
{
  return _statistics;
}

inline
void
SemiconductorModel::set_conduction_band_properties(double band_edge,
    double effective_mass, double mobility)
{
  _material.conduction_band_edge = band_edge;
  _material.electron_effective_mass = effective_mass;
  _material.electron_mobility = mobility;
}

inline
void
SemiconductorModel::set_valence_band_properties(double band_edge,
    double effective_mass, double mobility)
{
  _material.valence_band_edge = band_edge;
  _material.hole_effective_mass = effective_mass;
  _material.hole_mobility = mobility;
}

inline
void
SemiconductorModel::set_n_dopant(const Dopant& dopant)
{
  _material.n_dopant = dopant;
  _material.n_dopant.set_type(Dopant::N_TYPE);
}

inline
void
SemiconductorModel::set_p_dopant(const Dopant& dopant)
{
  _material.p_dopant = dopant;
  _material.p_dopant.set_type(Dopant::P_TYPE);
}

inline
void
SemiconductorModel::set_relative_permittivity(double epsilon_r)
{
  _material.relative_permittivity = epsilon_r;
}

inline
void
SemiconductorModel::add_recombination_model(RecombinationModels recomb_model)
{
  _recombination |= recomb_model;
}

inline
void
SemiconductorModel::remove_recombination_model(RecombinationModels recomb_model)
{
  _recombination &= !recomb_model;
}

inline
void
SemiconductorModel::set_SRH_parameters(double tau_n, double tau_p)
{
  _material.electron_recombination_time = tau_n;
  _material.hole_recombination_time = tau_p;
}

inline
double
SemiconductorModel::get_thermal_voltage(void) const
{
  return _thermal_voltage;
}

inline
double
SemiconductorModel::get_relative_permittivity(void) const
{
  return _material.relative_permittivity;
}

inline
double
SemiconductorModel::get_donor_density(void) const
{
  return _material.n_dopant.get_doping_density();
}

inline
double
SemiconductorModel::get_acceptor_density(void) const
{
  return _material.p_dopant.get_doping_density();
}


inline
double
SemiconductorModel::get_intrinsic_density(void) const
{
  return std::sqrt(get_intrinsic_density_squared());
}

inline
double
SemiconductorModel::get_intrinsic_density_squared(void) const
{
  return _equilibrium_electron_density * _equilibrium_hole_density;
}

inline
double
SemiconductorModel::get_band_gap(void) const
{
  return _conduction_band.band_edge - _valence_band.band_edge;
}

inline
double
SemiconductorModel::get_equilibrium_fermi_level(void) const
{
  return _equilibrium_fermi_level;
}

inline
double
SemiconductorModel::get_equilibrium_electron_density(void) const
{
  return _equilibrium_electron_density;
}

inline
double
SemiconductorModel::get_equilibrium_hole_density(void) const
{
  return _equilibrium_hole_density;
}

inline
const SemiconductorModel::BandProperties&
SemiconductorModel::get_conduction_band_properties(void) const
{
  return _conduction_band;
}

inline
const SemiconductorModel::BandProperties&
SemiconductorModel::get_valence_band_properties(void) const
{
  return _valence_band;
}

inline
SemiconductorModel::BandProperties&
SemiconductorModel::conduction_band_properties(void)
{
  return _conduction_band;
}


inline
SemiconductorModel::BandProperties&
SemiconductorModel::valence_band_properties(void)
{
  return _valence_band;
}

inline
const SemiconductorModel::MaterialDescriptor&
SemiconductorModel::get_material_descriptor(void) const
{
  return _material;
}

inline
double
SemiconductorModel::get_DOS_factor(void) const
{
  return _DOS_factor;
}

inline
void
SemiconductorModel::calculate_e_h_densities(double arg_e, double arg_h,
        double kT, CalculatedProperties& result)
{

  double Nc = _conduction_band.effective_DOS;
  double Nv = _valence_band.effective_DOS;

  if (_statistics == SimulationOptions::FERMIDIRAC)
  {
    density_and_derivatives<SimulationOptions::FERMIDIRAC>(arg_e,
        result.electron_density, result.electron_density_derivative,
        result.electron_density_2nd_derivative, result.dn_over_n);

    density_and_derivatives<SimulationOptions::FERMIDIRAC>(arg_h,
        result.hole_density, result.hole_density_derivative,
        result.hole_density_2nd_derivative, result.dp_over_p);
  }
  else
  {
    density_and_derivatives<SimulationOptions::BOLTZMANN>(arg_e,
        result.electron_density, result.electron_density_derivative,
        result.electron_density_2nd_derivative, result.dn_over_n);

    density_and_derivatives<SimulationOptions::BOLTZMANN>(arg_h,
        result.hole_density, result.hole_density_derivative,
        result.hole_density_2nd_derivative, result.dp_over_p);
  }
  
  result.electron_density *= Nc;
  result.electron_density_derivative *= Nc / kT;
  result.electron_density_2nd_derivative *= Nc / (kT * kT);
  result.dn_over_n /= kT;

  result.hole_density *= Nv;
  result.hole_density_derivative *= -Nv / kT;
  result.hole_density_2nd_derivative *= Nv / (kT * kT);
  result.dp_over_p /= -kT;
}

inline
void
SemiconductorModel::calculate_ionized_dopants(double arg_e, double arg_h,
        double kT, CalculatedProperties& result)
{
  
  const double arg_max = 150;
  
  // ionized donor density
  double Ed = _material.n_dopant.get_ionisation_energy();
  double arg = arg_e + Ed / kT;
  if (arg > arg_max)
  {
    result.ionized_donor_density = 0.0;
    result.ionized_donor_density_derivative = 0.0;
  }
  else
  {
    double tmp = std::exp(arg);
    double g  = _material.n_dopant.get_g_factor();
    double denom = 1 + g * tmp;
    double Nd = _material.n_dopant.get_doping_density();

    result.ionized_donor_density = Nd / denom;
    result.ionized_donor_density_derivative =
          -g * tmp * result.ionized_donor_density / (kT * denom);
  }

  // ionized acceptor density
  Ed = _material.p_dopant.get_ionisation_energy();
  arg = arg_h + Ed / kT;
  if (arg > arg_max)
  {
    result.ionized_acceptor_density = 0.0;
    result.ionized_acceptor_density_derivative = 0.0;
  }
  else
  {
    double tmp = std::exp(arg);
    double g  = _material.p_dopant.get_g_factor();
    double denom = 1 + g * tmp;
    double Na = _material.p_dopant.get_doping_density();

    result.ionized_acceptor_density = Na / denom;
    result.ionized_acceptor_density_derivative =
          g * tmp * result.ionized_acceptor_density / (kT * denom);
  }
}

inline
void
SemiconductorModel::calculate_charge_density(CalculatedProperties& result)
{
  double& n  = result.electron_density;
  double& p  = result.hole_density;
  double& Nd = result.ionized_donor_density;
  double& Na = result.ionized_acceptor_density;
  double& dn  = result.electron_density_derivative;
  double& dp  = result.hole_density_derivative;
  double& dNd = result.ionized_donor_density_derivative;
  double& dNa = result.ionized_acceptor_density_derivative;
  result.charge_density = p - n - Na + Nd;
  result.charge_density_derivatives[0] = dp - dn - dNa + dNd; // d / dpotential
  result.charge_density_derivatives[1] =    - dn       + dNd; // d / dfermi_e
  result.charge_density_derivatives[2] = dp      - dNa;       // d / dfermi_h
}

inline
void
SemiconductorModel::calculate_SRH_recombination(CalculatedProperties& result)
{
  double& n  = result.electron_density;
  double& p  = result.hole_density;
  double& dn  = result.electron_density_derivative;
  double& dp  = result.hole_density_derivative;
  double tn  = _material.electron_recombination_time;
  double tp  = _material.hole_recombination_time;
  double ni2 = get_intrinsic_density_squared();
  double ni  = std::sqrt(ni2);
  double denom = tp * (n + ni) + tn * (p + ni);
  double SRH = (n * p - ni2) / denom;
  double a = (p - tp * SRH) * dn / denom;
  double b = (n - tn * SRH) * dp / denom; 
  result.net_electron_recombination_rate += SRH;
  result.net_electron_recombination_rate_derivatives[1] += a;
  result.net_electron_recombination_rate_derivatives[2] += b;
  result.net_electron_recombination_rate_derivatives[0] += a + b;
  result.net_hole_recombination_rate += SRH;
  result.net_hole_recombination_rate_derivatives[1] += a;
  result.net_hole_recombination_rate_derivatives[2] += b;
  result.net_hole_recombination_rate_derivatives[0] += a + b;
}

inline
void
SemiconductorModel::calculate_Auger_recombination(CalculatedProperties& result)
{
}

inline
void
SemiconductorModel::calculate_direct_recombination(CalculatedProperties& result)
{
}

template<SimulationOptions::Statistics S>
inline
void
SemiconductorModel::density_and_derivatives(double arg, double& density,
    double& derivative, double& _2nd_derivative,
    double& derivative_over_density) const
{
  
  const double arg_max = 50;
  const double arg_min = -10;

  switch (S)
  {
    case SimulationOptions::FERMIDIRAC:
      if (arg < arg_max)
      {
        if (arg < arg_min)
        {
          density = std::exp(arg);
          derivative = density;
          _2nd_derivative = density;
          derivative_over_density = 1;
        }
        else
        {
          density = gsl_sf_fermi_dirac_half(arg);
          derivative = gsl_sf_fermi_dirac_mhalf(arg);
          // TODO implement d/dx F_-1/2(x)
          double step = 1e-2 * std::fabs(arg) + 1e-2;
          const double eps = 1e-6;
          double error = 1;
          _2nd_derivative = 0.5 * (gsl_sf_fermi_dirac_mhalf(arg + step)
              - gsl_sf_fermi_dirac_mhalf(arg - step)) / step;
          double old_d;
          while (std::fabs(error) > eps)
          {
            old_d = _2nd_derivative;
            step *= 0.5;
            _2nd_derivative = 0.5 * (gsl_sf_fermi_dirac_mhalf(arg + step)
                - gsl_sf_fermi_dirac_mhalf(arg - step)) / step;
            error = _2nd_derivative - old_d;
            if (std::fabs(_2nd_derivative) > eps)
              error /= _2nd_derivative;
          }
          derivative_over_density = derivative / density;
        }
      }
      else
      {
        density = 2 * M_2_SQRTPI / 3 * std::pow(arg, 1.5);
        derivative = M_2_SQRTPI * std::sqrt(arg);
        _2nd_derivative = 0.5 * M_2_SQRTPI / std::sqrt(arg);
        derivative_over_density = derivative / density;
      }
      break;

    default:
      if (arg < arg_max)
      {
        density = std::exp(arg);
        derivative = density;
        _2nd_derivative = density;
        derivative_over_density = 1;
      }
      else
      {
        density = std::exp(arg_max);
        derivative = density;
        _2nd_derivative = density;
        derivative_over_density = 1;
        //density = 2 * M_2_SQRTPI / 3 * std::pow(arg, 1.5);
        //derivative = M_2_SQRTPI * std::sqrt(arg);
        //_2nd_derivative = 0.5 * M_2_SQRTPI / std::sqrt(arg);
        //derivative_over_density = derivative / density;
      }
      break;
  }
}

#endif //_SEMICONDUCTOR_H_
