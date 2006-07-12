// $Id$

#ifndef _SIMPLESEMICONDUCTORMODEL_H_
#define _SIMPLESEMICONDUCTORMODEL_H_

#include "SimulationOptions.h"
#include "Dopant.h"
#include "DriftDiffusionProperties.h"
#include "Constants.h"

#include <vector>
#include <cmath>

// forward declarations
//class Point;
class Elem;

//! A simple semiconductor model
/*!
 * This simple model implements Boltzmann or Fermi-Dirac statistics,
 * constant mobility and SRH and direct recombination
 */
class SimpleSemiconductorModel : public DriftDiffusionProperties
{
  public:

    SimpleSemiconductorModel(void);
    SimpleSemiconductorModel(const SimpleSemiconductorModel& model);
    virtual ~SimpleSemiconductorModel(void) {};

    //! Set conduction band properties
    /*!
     * The DOS effective mass has to include the spin degeneracy.
     */
    void set_conduction_band_properties(double band_edge,
        double effective_mass, double mobility);

    //! Set valence band properties
    /*!
     * The DOS effective mass has to include the spin degeneracy.
     */
    void set_valence_band_properties(double band_edge,
        double effective_mass, double mobility);
    
    //! Set the relative permittivity
    void set_relative_permittivity(double epsilon_r);
    
    //! Add recombination model to be used.
    void add_recombination_model(DriftDiffusionDefs::RecombinationModel
        recomb_model);

    //! Remove a recombination model.
    void remove_recombination_model(DriftDiffusionDefs::RecombinationModel
        recomb_model);

    //! Set Shockley-Read-Hall recombination parameters
    void set_SRH_parameters(double tau_n, double tau_p);

    //! Set Direct Recombination parameters
    void set_direct_rec_parameters(double C)
      { _direct_rec_param = C; };


    virtual void read_database(const Dummy&) {};

    /*! \copydoc DriftDiffusionProperties::calculate_all()
     * 
     * This implementation models the most simple semiconductor equations
     */
    virtual void calculate_all(double potential,
      double fermi_e, double fermi_h, const Point& coord);
    
  protected:

    //! \copydoc DriftDiffusionProperties::prepare_element_data()
    /*!
     * This implementation calculates the effective density of states
     */
    virtual void prepare_element_data(void);

    //! Calculate Shockley-Read-Hall recombination
    void calculate_SRH_recombination(void);

    //! Calculate Auger recombination
    void calculate_Auger_recombination(void);

    //! Calculate direct recombination
    void calculate_direct_recombination(void);

  private:

    typedef DriftDiffusionProperties Parent;
    
    /**
     * The recombination models used
     */
    int _recombination;

    //! \c true if equilibrium properties are calculated
    bool _is_prepared;

    double _electron_recombination_time;
    double _hole_recombination_time;

    double _direct_rec_param;

};


//
// inline member functions
//


inline
void
SimpleSemiconductorModel::set_conduction_band_properties(double band_edge,
    double effective_mass, double mobility)
{
  get_conduction_band().band_edge = band_edge;
  get_conduction_band().effective_mass = effective_mass;
  electron_mobility = mobility;
}

inline
void
SimpleSemiconductorModel::set_valence_band_properties(double band_edge,
    double effective_mass, double mobility)
{
  get_valence_band().band_edge = band_edge;
  get_valence_band().effective_mass = effective_mass;
  hole_mobility = mobility;
}

inline
void
SimpleSemiconductorModel::set_relative_permittivity(double epsilon_r)
{
  permittivity = epsilon_r;
}

inline
void
SimpleSemiconductorModel::add_recombination_model(
    DriftDiffusionDefs::RecombinationModel recomb_model)
{
  _recombination |= recomb_model;
}

inline
void
SimpleSemiconductorModel::remove_recombination_model(
    DriftDiffusionDefs::RecombinationModel recomb_model)
{
  _recombination &= !recomb_model;
}

inline
void
SimpleSemiconductorModel::set_SRH_parameters(double tau_n, double tau_p)
{
  _electron_recombination_time = tau_n;
  _hole_recombination_time = tau_p;
}

inline
void
SimpleSemiconductorModel::calculate_SRH_recombination(void)
{
  double n  = electron_density;
  double p  = hole_density;
  double dn  = electron_density_derivative;
  double dp  = hole_density_derivative;
  double tn  = _electron_recombination_time;
  double tp  = _hole_recombination_time;
  double ni2 = get_intrinsic_density_squared();
  double ni  = std::sqrt(ni2);
  double denom = tp * (n + ni) + tn * (p + ni);
  double SRH = ni2 / denom;
  double nn = n / ni;
  double pp = p / ni;
  SRH *= (nn * pp - 1);
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
SimpleSemiconductorModel::calculate_direct_recombination(void)
{
  double n  = electron_density;
  double p  = hole_density;
  double dn  = electron_density_derivative;
  double dp  = hole_density_derivative;
  double C  = _direct_rec_param;
  double ni2 = get_intrinsic_density_squared();

  double rec = C * (n * p - ni2);
  double a = C * (dn * p);
  double b = C * (n * dp);
  electron_recombination_rate += rec;
  electron_recombination_rate_derivatives[1] += a;
  electron_recombination_rate_derivatives[2] += b;
  electron_recombination_rate_derivatives[0] += a + b;
  hole_recombination_rate += rec;
  hole_recombination_rate_derivatives[1] += a;
  hole_recombination_rate_derivatives[2] += b;
  hole_recombination_rate_derivatives[0] += a + b;
}
  

inline
void
SimpleSemiconductorModel::calculate_Auger_recombination(void)
{
}


#endif //_SIMPLESEMICONDUCTORMODEL_H_
