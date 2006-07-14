// $Id$

#ifndef _DRIFTDIFFUSIONPROPERTIES_H_
#define _DRIFTDIFFUSIONPROPERTIES_H_

//#include "tensor_value.h"
#include "vector_value.h"

#include "SimulationOptions.h"
#include "PhysicalProperties.h"
#include "DriftDiffusionDefs.h"
#include "TiberCad.h"
#include "Constants.h"
#include "TypeDefs.h"

// GNU scientific library
#include <gsl/gsl_sf_fermi_dirac.h>

#include <vector>
#include <set>
#include <map>


extern "C" {
  #include "petscsnes.h"
}

// forward declarations
class Point;
class Elem;
class Dopant;
class RecombinationModelInterface;

class DriftDiffusionProperties : public PhysicalProperties
{
    
  public:
       
    //! A default (empty) destructor.
    virtual ~DriftDiffusionProperties(void);

    //! Set the statistics to be used
    /*!
     * \param statistics the statistics
     */
    void set_statistics(TiberCad::Statistics statistics);

    //! Add a dopant
    void add_dopant(Dopant* dopant);

    //! Add a recombination model
    /*!
     * \param id the unique ID of this model
     */
    void add_recombination_model(RecombinationModelInterface* recomb_model);

    //! Get the statistics to be used
    /*!
     * \return the statistics
     */
    TiberCad::Statistics get_statistics(void) const;

    //! (Re-)Initialize for the given element
    /*!
     * \c reinit() calls \c prepare_element_data() which needs to be
     * implemented in derived classes
     */
    void reinit(const Elem* elem);
    
    //! Set the coupling type
    void set_coupling_type(DriftDiffusionDefs::Coupling coupling)
      { _coupling = (int) coupling; };
    
    //! Set the coupling type
    void set_coupling_type(int coupling)
      { _coupling = coupling; };
    
    //! Get the coupling type
    DriftDiffusionDefs::Coupling get_coupling_type(void) const
      { return (DriftDiffusionDefs::Coupling) _coupling; };

    //! Get the element we are currently working on
    const Elem* get_element(void) const;

    //! Get the coordinates of the point we are currently working on
    const Point* get_coordinates(void) const;

    //! Setup the band edge data
    /*!
     * This implementation calculates the effective density of states
     * and sets the band edges.
     */
    void setup_band_edges(void);

    //! Get the total n-doping
    double get_total_donor_density(void) const;

    //! Get the total p-doping
    double get_total_acceptor_density(void) const;

    //! Get the total net doping density
    /*!
     * The return value is \f$N_d - N_a\f$
     */
    double get_net_doping_density(void) const;


    //! Calculates the equilibrium properties.
    /*!
     *
     * This method has to be called before any call to \p calculate_all()
     * but after setting up all material parameters.
     * Call this method a derived one.
     *
     * \pre { \c reinit() has to be called before }
     * \post { all equilibrium properties are accessible without
     *  explicitly calling \c calculate_all() }
     */
    virtual void calculate_equilibrium_properties(
        int coupling = DriftDiffusionDefs::BOTH,
        double temperature = SimulationOptions::T);

    //! The method that will calculate all needed properties
    /*!
     * This method can be reimplemented in derived classes if necessary
     * 
     * \li the dielectric tensor
     * \li the total electric polarization
     * \li the thermal voltage \f$v_T=k_BT/e\f$ for the electrons
     * \li the thermal voltage \f$v_T=k_BT/e\f$ for the holes
     * \li the electron density and its derivative
     * \li the hole density and its derivative
     * \li the ionized donor density and its derivative
     * \li the ionized acceptor density and its derivative
     * \li the total charge density and its derivatives
     * \li the electron conductivity and its derivatives
     * \li the electron mobility
     * \li the hole mobility
     * \li the hole conductivity and its derivatives
     * \li the net electron recombination rate and its derivatives
     * \li the net hole recombination rate and its derivatives
     * 
     * \param potential the electric potential
     * \param fermi_e the electron electro-chemical potential
     * \param fermi_h the hole electro-chemical potential
     * \param p the coordinates in real space
     *
     */
    virtual void calculate_all(double potential,
      double fermi_e, double fermi_h, const Point& coord);
      

    //! Get the electron density
    /*!
     * Get the electron density as calculated by \c calculate_all(...)
     * 
     * \return the electron density
     */
    double get_electron_density(void) const
      { return electron_density; };
     
    //! Get the electron density derivative
    /*!
     * \return the electron density derivative with respect to the potential
     */
    double get_electron_density_derivative(void) const
      { return electron_density_derivative; };
    
    //! Get the hole density
    /*!
     * Get the hole density as calculated by \c calculate_all(...)
     * 
     * \return the hole density
     */
    double get_hole_density(void) const
      { return hole_density; };
     
    //! Get the ehole density derivative
    /*!
     * \return the hole density derivative with respect to the potential
     */
    double get_hole_density_derivative(void) const
      { return hole_density_derivative; };
    
    //! Get the ionized donor density
    /*!
     * Get the ionized donor density as calculated by \c calculate_all(...)
     * 
     * \return the ionized donor density
     */
    double get_ionized_donor_density(void) const
      { return ionized_donor_density; };
     
    //! Get the ionized donor density derivative
    /*!
     * \return the ionized donor density derivative with respect to the potential
     */
    double get_ionized_donor_density_derivative(void) const
      { return ionized_donor_density_derivative; };
        
    //! Get the ionized acceptor density
    /*!
     * Get the ionized acceptor density as calculated by \c calculate_all(...)
     * 
     * \return the ionized acceptor density
     */
    double get_ionized_acceptor_density(void) const
      { return ionized_acceptor_density; };
     
    //! Get the ionized acceptor density derivative
    /*!
     * \return the ionized acceptor density derivative with respect to the
     * potential
     */
    double get_ionized_acceptor_density_derivative(void) const
      { return ionized_acceptor_density_derivative; };
    
    
    //! Get the total charge density
    /*!
     * Get the total charge density as calculated by \c calculate_all(...)
     * \f$ \rho = p - n + N_D^+ - N_A^- \f$
     * 
     * \return the total charge density \f$\rho\f$
     *
     * \note
     * The charge density is returned in units of the elementary charge
     * \f$e\f$, \em not in Coulomb (= As)!
     */
    double get_charge_density(void) const;
    
    //! Get the derivatives of the charge density
    /*!
     * Get the derivatives of the total charge density with respect to the 
     * electric potential and the two electro-chemical potentials:
     * \f[
     *   \frac{\partial\rho}{\partial\varphi},\;
     *   \frac{\partial\rho}{\partial\phi_n},\;
     *   \frac{\partial\rho}{\partial\phi_p}
     * \f]
     * 
     * \return the derivatives as a const vector reference
     */
    const std::vector<double>& get_charge_density_derivatives(void) const;
    
    //! Get the derivatives of the charge density
    /*!
     * Get the derivatives of the total charge density with respect to the 
     * electric potential one of the two electro-chemical potentials:
     * \f{eqnarray*}
     *   \frac{\partial\rho}{\partial\varphi} & \mathsf{if}\, i=0, \\
     *   \frac{\partial\rho}{\partial\phi_n} & \mathsf{if}\, i=1, \\
     *   \frac{\partial\rho}{\partial\phi_p} & \mathsf{if}\, i=2
     * \f}
     */
    double get_charge_density_derivative(int i) const;
    
    //! Get the net electron recombination rate
    /*!
     * Get \f$R_{net} = R - G\f$ as
     * calculated by \c calculate_all(...)
     */
    double get_net_electron_recombination_rate(void) const
      { return electron_recombination_rate; };
      
    //! Get the net electron recombination rate derivative
    /*!
     * Get \f$\frac{\partial R_{net}}{\partial\varphi}\f$
     *
     * \return the derivatives as a const vector reference
     */
    const std::vector<double>&
      get_net_electron_recombination_rate_derivatives(void) const
        { return electron_recombination_rate_derivatives; };
    
    //! Get the net hole recombination rate
    /*!
     * Get \f$R_{net} = R - G\f$ as
     * calculated by \c calculate_all(...)
     *
     */
    double get_net_hole_recombination_rate(void) const
      { return hole_recombination_rate; };
      
    //! Get the net hole recombination rate derivative
    /*!
     * Get \f$\frac{\partial R_{net}}{\partial\varphi}\f$
     *
     * \return the derivatives as a const vector reference
     */
    const std::vector<double>&
      get_net_hole_recombination_rate_derivatives(void) const
        { return hole_recombination_rate_derivatives; };

    //! Get the total electric polarization
    /*!
     * The total electric polarization \b P is the sum of the
     * pyroelectric and piezoelectric polarization
     */
    const RealVectorValue& get_total_polarization(void) const
      { return polarization; };

    //! Get the relative permittivity tensor
    //const RealTensorValue& get_relative_permittivity(void) const
    double get_relative_permittivity(void) const
      { return permittivity; };
      
    //! Get the electron conductivity
    /*!
     * \return the electron conductivity \f$\sigma_n = \mu_n n\f$
     */
    double get_electron_conductivity(void) const
      { return electron_conductivity; };
      
    //! Get the hole conductivity
    /*!
     * \return the hole conductivity \f$\sigma_p = \mu_p p\f$
     */
    double get_hole_conductivity(void) const
      { return hole_conductivity; };
      
    //! Get the electron mobility
    /*!
     * \return the electron mobility
     */
    double get_electron_mobility(void) const
      { return electron_mobility; };
      
    //! Get the hole mobility
    /*!
     * \return the hole mobility
     */
    double get_hole_mobility(void) const
      { return hole_mobility; };


    //! Get the electron conductivity derivatives
    const std::vector<double>& get_electron_conductivity_derivatives(void) const
      { return electron_conductivity_derivatives; };
      
    //! Get the hole conductivity derivatives
    const std::vector<double>& get_hole_conductivity_derivatives(void) const
      { return hole_conductivity_derivatives; };

    //! Get the equilibrium electron density
    double get_equilibrium_electron_density(void) const
      { return equilibrium_electron_density; };

    //! Get the equilibrium hole density
    double get_equilibrium_hole_density(void) const
      { return equilibrium_hole_density; };

    //! Get the square of the intrinsic density
    double get_intrinsic_density_squared(void) const
      { return equilibrium_electron_density * equilibrium_hole_density; };

    //! Get the intrinsic density
    double get_intrinsic_density(void) const
      { return std::sqrt(get_intrinsic_density_squared()); };

    //! Get equilibrium fermi level
    double get_equilibrium_fermi_level(void) const
      { return equilibrium_fermi_level; };

    //! Get the conduction band edge
    double get_conduction_band_edge(void) const
      { return conduction_band.band_edge; };

    //! Get the valence band edge
    double get_valence_band_edge(void) const
      { return valence_band.band_edge; };

    //! Get the band gap
    double get_band_gap(void) const
      { return conduction_band.band_edge - valence_band.band_edge; };

    void get_net_recombination_rates(std::vector<double>& rates) {};

    RecombinationModelInterface* get_recombination_model(ID id);

    //! clear all doping
    void clear_doping(void);

    //! clear all recombination rates
    void clear_recombination(void);

  protected:
      
    /*!
     * This structure holds the basic properties of a band for given
     * conditions (temp etc.)
     */
    struct BandProperties
    {
      //! The effective mass for the DOS
      /*!
       * It includes any degeneration, i.e. also spin
       */
      double effective_mass;
      
      //! The effective density of states
      double effective_DOS;

      //! The band edge energy
      double band_edge;
    };

    //! The empty constructor.
    DriftDiffusionProperties(void);


    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     * This method can be used overiden by derived classes.
     */
    virtual void prepare_element_data(void) {};
    

    //! The thermal voltage for the electrons
    double electron_vt;

    //! The thermal voltage for the holes
    double hole_vt;

    //! The electron density
    double electron_density;

    //! The electron density derivative
    double electron_density_derivative;

    //! The hole density
    double hole_density;

    //! The hole density derivative
    double hole_density_derivative;
  
    //! The ionized donor density
    double ionized_donor_density;

    //! The ionized donor density derivative
    double ionized_donor_density_derivative;
        
    //! The ionized acceptor density
    double ionized_acceptor_density;

    //! The ionized acceptor density derivative
    double ionized_acceptor_density_derivative;

    //! The total charge density
    double charge_density;

    //! The total charge density derivatives
    std::vector<double> charge_density_derivatives;

    //! The electron mobility
    double electron_mobility;

    //! The hole mobility
    double hole_mobility;
    
    //! The electron conductivity
    double electron_conductivity;
    
    //! The derivatives of the electron conductivity
    std::vector<double> electron_conductivity_derivatives;
    
    //! The hole conductivity
    double hole_conductivity;
    
    //! The derivatives of the hole conductivity
    std::vector<double> hole_conductivity_derivatives;
    
    //! The net electron recombination rate
    double electron_recombination_rate;
    
    //! The derivatives of the net electron recombination rate
    std::vector<double> electron_recombination_rate_derivatives;
    
    //! The net hole recombination rate
    double hole_recombination_rate;
    
    //! The derivatives of the net hole recombination rate
    std::vector<double> hole_recombination_rate_derivatives;
    
    //! The total electric polarization
    RealVectorValue polarization;
    
    //! The relative permittivity tensor
    //RealTensorValue permittivity;
    double permittivity;

    //! The equilibrium fermi level
    /*!
     * The fermi level such that \f$n=n_0,\,p=p_0\f$
     */
    double equilibrium_fermi_level;

    //! The equilibrium electron density
    /*!
     * \f$n_i = n_0 p_0\f$
     */
    double equilibrium_electron_density;

    //! The equilibrium hole density
    /*!
     * \f$n_i = n_0 p_0\f$
     */
    double equilibrium_hole_density;
        
    //! Calculate the density and its derivatives
    /*! 
     * The method also returns the ratio 'first derivative over density'
     * which is used for the mobility.
     */
    template<TiberCad::Statistics S>
    void density_and_derivatives(double arg, double& density,
        double& derivative, double& _2nd_derivative,
        double& derivative_over_density) const;
    
    //! Calculate the density for a given argument
    /*!
     * Basically an approximation for the Fermi integral with index 0.5 is
     * returned if statistics is Fermi-Dirac
     */
    template<TiberCad::Statistics S> double density(double arg) const;

    //! Get the conduction band properties
    const BandProperties& get_conduction_band(void) const
      { return conduction_band; };

    //! Get the valence band properties
    const BandProperties& get_valence_band(void) const
      { return valence_band; };

    //! Get the conduction band properties
    BandProperties& get_conduction_band(void)
      { return conduction_band; };

    //! Get the valence band properties
    BandProperties& get_valence_band(void)
      { return valence_band; };

    //! Get the constant factor to calculate the effective density of states
    /*!
     * \return the factor pow(2 * PI / h^2)^1.5
     *
     * The spin degeneracy has to be included in the effective mass.
     */
    static double get_DOS_factor(void)
      { return _DOS_factor; }

    static PetscErrorCode jacobian(SNES snes, Vec x,
        Mat *jac, Mat *B, MatStructure *flag, void *sc);

    static PetscErrorCode function(SNES snes, Vec x, Vec f, void *sc);


  private:

    typedef std::set<Dopant*>::iterator dopant_iterator;
    typedef
      std::map<ID, RecombinationModelInterface*>::iterator recomb_iterator;

    //! The copy constructor
    DriftDiffusionProperties(const DriftDiffusionProperties& rhs);

    //! The element we are currently working on
    const Elem* _elem;

    //! The coordinates of the point we are working on
    const Point* _coord;

    //! The statistics used 
    TiberCad::Statistics _statistics;

    //! Type of coupling (particles) we want to study
    /*!
     * This can be one of \c ELECTRONS, \c HOLES or \c BOTH
     */
    int _coupling;

    //! The conduction band properties
    /*!
     * Band properties are assumed to be elemental data, \em not nodal data
     */
    BandProperties conduction_band;

    //! The conduction band properties
    /*!
     * Band properties are assumed to be elemental data, \em not nodal data
     */
    BandProperties valence_band;

    //! The donors
    std::set<Dopant*> _donors;

    //! The acceptors
    std::set<Dopant*> _acceptors;

    //! The recombination models
    std::map<ID, RecombinationModelInterface*> _recombination_models;

    //MobilityModel* _mobility_model;

    //! The constant factor to calculate the effective density of states
    /*!
     * The spin degeneracy has to be included in the effective mass
     */
    static const double _DOS_factor;


};


//
// inline members
//

inline
void
DriftDiffusionProperties::reinit(const Elem* elem)
{
  _elem = elem;
  this->prepare_element_data();
}
    
inline
const Elem*
DriftDiffusionProperties::get_element(void) const
{
  return _elem;
}
    
inline
const Point*
DriftDiffusionProperties::get_coordinates(void) const
{
  return _coord;
}


inline
double
DriftDiffusionProperties::get_charge_density(void) const
{
  return charge_density;
}

inline
const std::vector<double>&
DriftDiffusionProperties::get_charge_density_derivatives(void) const
{
  return charge_density_derivatives;
}

inline
RecombinationModelInterface*
DriftDiffusionProperties::get_recombination_model(ID id)
{
  RecombinationModelInterface* rec = NULL;
  recomb_iterator it = _recombination_models.find(id);
  if (it != _recombination_models.end())
    rec = it->second;

  return rec;
}

inline
double
DriftDiffusionProperties::get_charge_density_derivative(int i) const
{
  double drho;
  
  switch (i)
  {
    case 1:
      drho = charge_density_derivatives[1];
      break;
    case 2:
      drho = charge_density_derivatives[2];
      break;
    default: // i = 0
      drho = charge_density_derivatives[0];
      break;
  }
  
  return drho;
}

inline
void
DriftDiffusionProperties::set_statistics(TiberCad::Statistics statistics)
{
  _statistics = statistics;
}

inline
TiberCad::Statistics
DriftDiffusionProperties::get_statistics(void) const
{
  return _statistics;
}

inline
double
DriftDiffusionProperties::get_net_doping_density(void) const
{
  return (get_total_donor_density() - get_total_acceptor_density());
}

template<TiberCad::Statistics S>
inline
double
DriftDiffusionProperties::density(double arg) const
{
  
  const double arg_max = 150;
  const double arg_min = -50;

  double dens;
  switch (S)
  {
    case TiberCad::FERMIDIRAC:
      if (arg < arg_max)
      {
        if (arg < arg_min)
          dens = std::exp(arg);
        else
          dens = gsl_sf_fermi_dirac_half(arg);
      }
      else
        dens = 2 * M_2_SQRTPI / 3 * std::pow(arg, 1.5);

      break;

    default:
      if (arg < arg_max)
        dens = std::exp(arg);
      else
        dens = std::exp(arg_max);
  }
  return dens;
}


template<TiberCad::Statistics S>
inline
void
DriftDiffusionProperties::density_and_derivatives(double arg, double& density,
    double& derivative, double& _2nd_derivative,
    double& derivative_over_density) const
{
  
  const double arg_max = 150;
  const double arg_min = -50;

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
      arg = (arg > arg_max) ? arg_max : arg;

      density = std::exp(arg);
      derivative = density;
      _2nd_derivative = density;
      derivative_over_density = 1;
      break;
  }
}

inline
void
DriftDiffusionProperties::setup_band_edges(void)
{
  double kT = SimulationOptions::T * Constants::k_B;
  electron_vt = hole_vt = kT;
  
  BandProperties& cb = conduction_band;
  BandProperties& vb = valence_band;

  cb.effective_DOS =
    get_DOS_factor() * std::pow(kT * cb.effective_mass, 1.5);

  vb.effective_DOS =
    get_DOS_factor() * std::pow(kT * vb.effective_mass, 1.5);
}






#endif /* _DRIFTDIFFUSIONPROPERTIES_H_*/
