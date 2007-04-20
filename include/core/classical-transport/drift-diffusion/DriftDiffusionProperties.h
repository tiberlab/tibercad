// $Id$

#ifndef _DRIFTDIFFUSIONPROPERTIES_H_
#define _DRIFTDIFFUSIONPROPERTIES_H_


#ifndef TIBER_MODULE_NAME
# define TIBER_MODULE_NAME dd
#endif


#include "PhysicalModel.h"

#include "SimulationOptions.h"
#include "DriftDiffusionDefs.h"
#include "TiberCad.h"
#include "Constants.h"
#include "TypeDefs.h"

#include "tensor.h"

//#include "tensor_value.h"
#include "vector_value.h"

// GNU scientific library
#include <gsl/gsl_sf_fermi_dirac.h>
#include <values.h>

#include <vector>
#include <set>
#include <map>


// forward declarations
class Point;
class Elem;
class Dopant;
class RecombinationModelInterface;
class MobilityModelInterface;

//! The base class for all drift-diffusion related semiconductor models
/*!
 * \note { it is not yet possible to have twice the same recombination model.
 * Trying to add to identical models will result in a memory leak. This will
 * be corrected in future. }
 */
class DriftDiffusionProperties : public PhysicalModel
{
    
  public:
       
    //! A default (empty) destructor.
    virtual ~DriftDiffusionProperties(void);


    //! Create a named drift-diffusion model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     * 
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static DriftDiffusionProperties* create(const std::string& name,
        const ModelOptions& options = ModelOptions());


    //! Set the statistics to be used
    /*!
     * \param statistics the statistics
     */
    void set_statistics(TiberCad::Statistics statistics);


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

    
    //! Setup the band edge data
    /*!
     * This implementation calculates the effective density of states
     * and sets the band edges.
     */
    void setup_band_edges(void);


    //! Set the coordinates
    void set_coordinates(const Point& p);

    //! Set the carrier temperatures
    /*!
     * The lattice, electron and hole temperatures have to be set before
     * calling any of the \c calculate_xxx methods
     *
     * \param T_lat the lattice temperature
     * \param T_e the electron temperature
     * \param T_h the hole temperature
     */
    void set_carrier_temperatures(double T_e, double T_h);

    //! Set the potentials
    /*!
     * The potentials have to be set before calling any of the
     * \c calculate_xxx methods. If the densities are used as variables,
     * the electro-chemical potentials can be left out.
     *
     * \param potential the electric potential
     * \param Ef_e the electron electro-chemical potential
     * \param Ef_h the hole electro-chemical potential
     */
    void set_potentials(double potential, double Ef_e = 0.0, double Ef_h = 0.0);

    //! Set the carrier densities
    /*!
     * This method can be used if the densities are the independent variables.
     * The densities have to be set before calling any of the
     * \c calculate_xxx methods.
     *
     * \param n the electron density
     * \param p the hole density
     */
    void set_densities(double n, double p);

    //! Set the electric field
    /*!
     * The electric field has to be set before calling any of the
     * \c calculate_xxx methods.
     *
     * \param E the electric field
     */
    void set_electric_field(RealGradient E);


    //! Get the element we are currently working on
    const Elem* get_element(void) const;

    //! Get the coordinates of the point we are currently working on
    const Point& get_coordinates(void) const;

    //! Set the strain
    void set_strain(const Tensor2Sym& strain);

    //! Get the strain
    const Tensor2Sym& get_strain(void) const;


    //! Set the lattice temperature (in K)
    void set_lattice_temperature(double T);
    
    //! Get the lattice temperature (in units of eV)
    double get_lattice_temperature(void) const;


    //! Calculate the electro-chemical potentials for given densities
    void calculate_electro_chemical_potentials(void);
    
    
    //! Calculate electron and hole densities and derivatives
    /*!
     * This method calculates electron and hole densities and their
     * derivatives with respect to the electric and the
     * electro-chemical potentials.
     */
    void calculate_densities(void);

    
    //! Calculate the ionized dopant densities and derivatives
    /*!
     * The total density of ionized donors and acceptors is calculated,
     * respectively. Their derivative with respect to the electric
     * potential are also computed (the derivative with respect to the
     * electro-chemical potentials have the same value but opposite sign).
     */
    void calculate_ionized_dopants(void);

    
    //! Calculate net recombination rates and derivatives
    /*!
     * The recombination models calculate the derivatives with respect to
     * the densities.
     */
    void calculate_net_recombination_rates(void);

    
    //! Calculate the mobilities
    void calculate_mobilities(void);

    
    //! Calculates the equilibrium properties.
    /*!
     *
     * This method has to be called before calculating anything
     * Call this method from a derived one.
     *
     * \pre { \c reinit() has to be called before }
     * \post { all equilibrium properties are accessible without
     *  explicitly calling \c calculate_all() }
     */
    virtual void calculate_equilibrium_properties(void);
    

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
     * \return the electron density derivative with respect to the
     * electric potential
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
     * \return the hole density derivative with respect to the
     * electric potential
     */
    double get_hole_density_derivative(void) const
      { return hole_density_derivative; };
    

    //! Get the ionized donor density
    /*!
     * \return the ionized donor density
     */
    double get_ionized_donor_density(void) const
      { return ionized_donor_density; };
     
    
    //! Get the ionized donor density derivative
    /*!
     * \return the ionized donor density derivative with respect to the
     * electric potential
     */
    double get_ionized_donor_density_derivative(void) const
      { return ionized_donor_density_derivative; };
        
    
    //! Get the ionized acceptor density
    /*!
     * \return the ionized acceptor density
     */
    double get_ionized_acceptor_density(void) const
      { return ionized_acceptor_density; };
     
    
    //! Get the ionized acceptor density derivative
    /*!
     * \return the ionized acceptor density derivative with respect to the
     * electric potential
     */
    double get_ionized_acceptor_density_derivative(void) const
      { return ionized_acceptor_density_derivative; };
    
    
    //! Get the total charge density
    /*!
     *
     * \f$ \rho = p - n + N_D^+ - N_A^- \f$
     * 
     * \return the total charge density \f$\rho\f$
     *
     * \note
     * The charge density is returned in units of the elementary charge
     * \f$e\f$, \em not in Coulomb (= As)!
     */
    double get_charge_density(void) const;
  
    
    //! Get the net electron recombination rate
    /*!
     * Get \f$R_{net} = R - G\f$ as
     */
    double get_net_electron_recombination_rate(void) const
      { return electron_recombination_rate; };
      
    
    //! Get the net electron recombination rate derivatives
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
      
    
    //! Get the net hole recombination rate derivatives
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
    //double get_electron_conductivity(void) const
    //  { return electron_conductivity; };
      
    
    //! Get the hole conductivity
    /*!
     * \return the hole conductivity \f$\sigma_p = \mu_p p\f$
     */
    //double get_hole_conductivity(void) const
    //  { return hole_conductivity; };
      
    
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
    //const std::vector<double>& get_electron_conductivity_derivatives(void) const
    //  { return electron_conductivity_derivatives; };
      
    
    //! Get the hole conductivity derivatives
    //const std::vector<double>& get_hole_conductivity_derivatives(void) const
    //  { return hole_conductivity_derivatives; };


    //! Get the square of the intrinsic density
    double get_intrinsic_density_squared(void) const
      { return intrinsic_density * intrinsic_density; };

    
    //! Get the intrinsic density
    double get_intrinsic_density(void) const
      { return intrinsic_density; };

    
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

    void get_net_recombination_rates(std::vector<double>& rates);
    
    
    //! Get the IDs of the registered recombination models
    /*!
     * \return number of registered models
     */
    int get_net_recombination_rate_IDs(std::vector<ID>& ids);

    
    //! Get the recombination model with ID \c id
    RecombinationModelInterface* get_recombination_model(ID id);

    
    //! Get the recombination model with ID \c id
    const RecombinationModelInterface* get_recombination_model(ID id) const;

    
    //! get the net recombination rate of model \c id
    double get_net_recombination_rate(ID id);

    
    //! Clear all recombination rates
    void clear_recombination(void);

    
    //! Returns the number of recombination models
    int get_number_of_recombination_models(void) const;


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

    //! Initialize this model
    /*!
     * This reads the database and calls init for all submodels
     * A derived class which reimplements this method has to call 
     * explicitly the one of this class!
     */
    virtual void do_init(void);

    
    /*! \copydoc PhysicalModel::calculate_VCA() */
    virtual void calculate_VCA(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);


    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     * This method can be used overiden by derived classes.
     */
    virtual void prepare_element_data(void) {};

    
    //! \copydoc PhysicalModel::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    
    //! \copydoc PhysicalModel::copy_from()
    virtual void copy_from(const PhysicalModelInterface* rhs);

    
    //! The lattice temperature in eV (\f$= k_B T_{lat} / e\f$)
    double lattice_vt;

    //! The electron temperature in eV (\f$= k_B T_e / e\f$)
    double electron_vt;

    //! The hole temperature in eV (\f$= k_B T_h / e\f$)
    double hole_vt;


    //! The electric poential
    double electric_potential;

    //! The electron electro-chemical potential
    double fermi_e;

    //! The hole electro-chemical potential
    double fermi_h;


    //! The electric field
    RealGradient electric_field;


    //! The electron density
    /*!
     * The electron density is given by
     * \f[n = N_c e^{\frac{e\varphi - e\phi_n - E_c}{k_B T_e}}\f]
     */
    double electron_density;

    //! The electron density derivative with respect to the electric potential
    /*!
     * The derivative with respect to the electro-chemical potential
     * has the same value but opposite sign.
     */
    double electron_density_derivative;

    //! The hole density
    /*!
     * The hole density is given by
     * \f[p = N_v e^{-\frac{e\varphi - e\phi_p - E_v}{k_B T_h}}\f]
     */
    double hole_density;

    //! The hole density derivative
    /*!
     * The derivative with respect to the electro-chemical potential
     * has the same value but opposite sign.
     */
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

    //! The intrinsic density
    double intrinsic_density;
        

    
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



  private:

    //! An iterator for the recombination models
    typedef std::map<ID, RecombinationModelInterface*>::iterator
      recomb_iterator;
    
    //! A const iterator for the recombination models
    typedef std::map<ID, RecombinationModelInterface*>::const_iterator
      const_recomb_iterator;

    
    //! The copy constructor is disabled
    DriftDiffusionProperties(const DriftDiffusionProperties& rhs);

    
    //! The assignment operator is disabled
    DriftDiffusionProperties& operator=(const DriftDiffusionProperties& rhs);

    
    //! Add a recombination model
    /*!
     * Creates and adds a new recombination model from the given name and
     * options.
     *
     * \param model_name the name of the model
     * \param options the options for the model
     */
     void add_recombination_model(const std::string& model_name,
        const ModelOptions& options = ModelOptions());

    
    
    //! Create a mobility model
    /*!
     * Creates a mobility model from the given model name
     * and options.
     */
    MobilityModelInterface* create_mobility_model(
        const ModelOptions& options = ModelOptions());



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
    
    //! The strain
    Tensor2Sym _strain;

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

    
    //! The recombination models
    std::map<ID, RecombinationModelInterface*> _recombination_models;

    //! The electron mobility
    MobilityModelInterface* _electron_mobility;

    //! The hole mobility
    MobilityModelInterface* _hole_mobility;

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
DriftDiffusionProperties*
DriftDiffusionProperties::create(const std::string& name,
    const ModelOptions& options)
{
  return dynamic_cast<DriftDiffusionProperties*>(
      PhysicalModelInterface::create("dd_" + name, options));
}

inline
void
DriftDiffusionProperties::reinit(const Elem* elem)
{
  if (_elem != elem)
  {
    _elem = elem;
    this->prepare_element_data();
  }
}
 
inline
void
DriftDiffusionProperties::set_coordinates(const Point& p)
{
  _coord = &p;
}

inline
void
DriftDiffusionProperties::set_carrier_temperatures(double T_e, double T_h)
{
  electron_vt = T_e;
  hole_vt = T_h;
}

inline
void
DriftDiffusionProperties::set_potentials(double potential, double Ef_e,
    double Ef_h)
{
  electric_potential = potential;
  fermi_e = Ef_e;
  fermi_h = Ef_h;
}

inline
void
DriftDiffusionProperties::set_electric_field(RealGradient E)
{
  electric_field = E;
}

inline
void
DriftDiffusionProperties::set_densities(double n, double p)
{
  electron_density = n;
  hole_density = p;
}



inline
void
DriftDiffusionProperties::set_lattice_temperature(double T)
{
  lattice_vt = T * Constants::k_B;
}



inline
double
DriftDiffusionProperties::get_lattice_temperature(void) const
{
  return lattice_vt;
}



inline
const Elem*
DriftDiffusionProperties::get_element(void) const
{
  return _elem;
}
    
inline
const Point&
DriftDiffusionProperties::get_coordinates(void) const
{
  return *_coord;
}
 

inline
double
DriftDiffusionProperties::get_charge_density(void) const
{
  return hole_density - electron_density + ionized_donor_density -
    ionized_acceptor_density;
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
const RecombinationModelInterface*
DriftDiffusionProperties::get_recombination_model(ID id) const
{
  RecombinationModelInterface* rec = NULL;
  const_recomb_iterator it = _recombination_models.find(id);
  if (it != _recombination_models.end())
    rec = it->second;

  return rec;
}


inline
void
DriftDiffusionProperties::set_strain(const Tensor2Sym& strain)
{
  _strain = strain;
}


inline
const Tensor2Sym&
DriftDiffusionProperties::get_strain(void) const
{
  return _strain;
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
        dens = 2.0 * M_2_SQRTPI / 3.0 * std::pow(arg, 1.5);

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
  
  //const double arg_max = 150;
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
        density = 2.0 * M_2_SQRTPI / 3.0 * std::pow(arg, 1.5);
        derivative = M_2_SQRTPI * std::sqrt(arg);
        _2nd_derivative = 0.5 * M_2_SQRTPI / std::sqrt(arg);
        derivative_over_density = derivative / density;
      }
      break;

    default:
      arg = (arg > arg_max) ? arg_max : arg;

      density = std::exp(arg) + 2 * MINDOUBLE;
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
  double kT = lattice_vt;
  electron_vt = hole_vt = kT;
  
  BandProperties& cb = conduction_band;
  BandProperties& vb = valence_band;

  cb.effective_DOS =
    get_DOS_factor() * std::pow(kT * cb.effective_mass, 1.5);

  vb.effective_DOS =
    get_DOS_factor() * std::pow(kT * vb.effective_mass, 1.5);
}


inline
int
DriftDiffusionProperties::get_number_of_recombination_models(void) const
{
  return _recombination_models.size();
}


inline
PhysicalModelInterface*
DriftDiffusionProperties::create_new(void) const
{
  return new DriftDiffusionProperties();
}



#endif /* _DRIFTDIFFUSIONPROPERTIES_H_ */
