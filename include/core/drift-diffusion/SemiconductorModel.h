// $Id$

#ifndef _SEMICONDUCTOR_H_
#define _SEMICONDUCTOR_H_

#include "SimulationOptions.h"
#include "Dopant.h"
#include "DriftDiffusionProperties.h"
#include "TiberCad.h"

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

class SemiconductorModel : public DriftDiffusionProperties
{
  public:

    /**
     * This contains all reference material properties for this model
     */
    struct MaterialDescriptor
    {

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


    SemiconductorModel(void);
    SemiconductorModel(const SemiconductorModel& model);
    virtual ~SemiconductorModel(void) {};

    /**
     * Set the material for which the properties should be calculated.
     */
    void set_material_descriptor(const MaterialDescriptor& desc);


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
    void add_recombination_model(
        DriftDiffusionDefs::RecombinationModels recomb_model);

    /**
     * Remove a recombination model.
     */
    void remove_recombination_model(
        DriftDiffusionDefs::RecombinationModels recomb_model);

    /**
     * Set Shockley-Read-Hall recombination parameters
     */
    void set_SRH_parameters(double tau_n, double tau_p);


    /**
     * Calculates the equilibrium properties.
     *
     * This method has to be called before any call to \p calculate_all()
     * but after setting up all material parameters
     */
    virtual void calculate_equilibrium_properties(double temperature =
        SimulationOptions::T);

    virtual void read_database(const Dummy&) {};

    /**
     * Assume strained semiconductor.
     *
     * This method has to be called before any call to \p calculate_all()
     * It has to be implemented if strain has to be modeled.
     */
    //virtual void calculate_strained_properties(const Elem* elem) {};
    virtual void calculate_strained_properties(Elem* elem) {};
        
    /*! \copydoc DDSemiconductor::calculate_all()
     * 
     * This implementation models the most simple semiconductor equations
     */
    virtual void calculate_all(double potential,
      double fermi_e, double fermi_h,
      const Point& p, const Elem* elem, int coupling);
    
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
    //double get_intrinsic_density(void) const;

    /**
     * @returns the square of the intrinsic density n_i
     */
    //double get_intrinsic_density_squared(void) const;

    /**
     * @returns the band gap
     */
    double get_band_gap(void) const;

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
    
    void SemiconductorModel::calculate_ionized_donors(double arg_e, double kT,
        double& Nd, double& dNd);

    void SemiconductorModel::calculate_ionized_acceptors(double arg_h, double kT,
        double& Na, double& dNa);

    /**
     * Calculate Shockley-Read-Hall recombination
     */
    void calculate_SRH_recombination(void);

    /**
     * Calculate Auger recombination
     */
    void calculate_Auger_recombination(void);

    /**
     * Calculate direct recombination
     */
    void calculate_direct_recombination(void);

  private:

    /**
     * The recombination models used
     */
    int _recombination;

    /**
     * The factor 2 * pow(2 * PI / h^2)^1.5
     * for calculating the effective density of states
     */
    double _DOS_factor;

    /**
     * The band properties
     */
    BandProperties _conduction_band;
    BandProperties _valence_band;

    /**
     * The basic (reference) material properties
     */
    MaterialDescriptor _material;

    //! The n-type doping
    Dopant _n_dopant;
    
    //! The p-type doping
    Dopant _p_dopant;


    /**
     * Calculate the density and its derivatives
     * 
     * The method also returns the ratio 'first derivative over density'
     * which is used for the mobility.
     */
    template<TiberCad::Statistics S>
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
void
SemiconductorModel::set_material_descriptor(const MaterialDescriptor& desc)
{
  _material = desc;
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
  _n_dopant = dopant;
  _n_dopant.set_type(Dopant::N_TYPE);
}

inline
void
SemiconductorModel::set_p_dopant(const Dopant& dopant)
{
  _p_dopant = dopant;
  _p_dopant.set_type(Dopant::P_TYPE);
}

inline
void
SemiconductorModel::set_relative_permittivity(double epsilon_r)
{
  permittivity = epsilon_r;
}

inline
void
SemiconductorModel::add_recombination_model(
        DriftDiffusionDefs::RecombinationModels recomb_model)
{
  _recombination |= recomb_model;
}

inline
void
SemiconductorModel::remove_recombination_model(
        DriftDiffusionDefs::RecombinationModels recomb_model)
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
SemiconductorModel::get_donor_density(void) const
{
  return _n_dopant.get_doping_density();
}

inline
double
SemiconductorModel::get_acceptor_density(void) const
{
  return _p_dopant.get_doping_density();
}

inline
double
SemiconductorModel::get_band_gap(void) const
{
  return _conduction_band.band_edge - _valence_band.band_edge;
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
SemiconductorModel::calculate_ionized_donors(double arg_e, double kT,
    double& Nd, double& dNd)
{
  const double arg_max = 150;

  double Ed = _n_dopant.get_ionisation_energy();
  double arg = arg_e + Ed / kT;
  if (arg > arg_max)
  {
    Nd = 0.0;
    dNd = 0.0;
  }
  else
  {
    double tmp = std::exp(arg);
    double g  = _n_dopant.get_g_factor();
    double denom = 1 + g * tmp;

    Nd = _n_dopant.get_doping_density() / denom;
    dNd = -g * tmp * Nd / (kT * denom);
  }
}

inline
void
SemiconductorModel::calculate_ionized_acceptors(double arg_h, double kT,
    double& Na, double& dNa)
{
  const double arg_max = 150;

  double Ed = _p_dopant.get_ionisation_energy();
  double arg = arg_h + Ed / kT;
  if (arg > arg_max)
  {
    Na = 0.0;
    dNa = 0.0;
  }
  else
  {
    double tmp = std::exp(arg);
    double g  = _p_dopant.get_g_factor();
    double denom = 1 + g * tmp;

    Na = _p_dopant.get_doping_density() / denom;
    dNa = g * tmp * Na / (kT * denom);
  }

}


inline
void
SemiconductorModel::calculate_SRH_recombination(void)
{
  double& n  = electron_density;
  double& p  = hole_density;
  double& dn  = electron_density_derivative;
  double& dp  = hole_density_derivative;
  double tn  = _material.electron_recombination_time;
  double tp  = _material.hole_recombination_time;
  double ni2 = get_intrinsic_density_squared();
  double ni  = std::sqrt(ni2);
  double denom = tp * (n + ni) + tn * (p + ni);
  double SRH = (n * p - ni2) / denom;
  double a = (p - tp * SRH) * dn / denom;
  double b = (n - tn * SRH) * dp / denom; 
  electron_recombination_rate += SRH;
  electron_recombination_rate_derivatives[1] += a;
  electron_recombination_rate_derivatives[2] += b;
  electron_recombination_rate_derivatives[0] += a + b;
  hole_recombination_rate += SRH;
  hole_recombination_rate_derivatives[1] += a;
  hole_recombination_rate_derivatives[2] += b;
  hole_recombination_rate_derivatives[0] += a + b;
}

inline
void
SemiconductorModel::calculate_Auger_recombination(void)
{
}

inline
void
SemiconductorModel::calculate_direct_recombination(void)
{
}


template<TiberCad::Statistics S>
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
    case TiberCad::FERMIDIRAC:
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
