// $Id$

#ifndef _DRIFTDIFFUSIONPROPERTIES_H_
#define _DRIFTDIFFUSIONPROPERTIES_H_


#ifndef TIBER_MODULE_NAME
# define TIBER_MODULE_NAME dd
#endif


#include "PhysicalModel.h"
#include "Variable.h"

#include "ParticleDensity.h"
#include "TemperatureInterface.h"
#include "StrainInterface.h"
#include "SimulationOptions.h"
#include "DriftDiffusion.h"
#include "DriftDiffusionDefs.h"
#include "TiberCad.h"
#include "Constants.h"
#include "TypeDefs.h"

#include "tensor.h"

//#include "tensor_value.h"
#include "vector_value.h"
#include "point.h"


#include <float.h>
#include <vector>
#include <set>
#include <map>
#include <stack>


// forward declarations
class Elem;
class Dopant;
class RecombinationModelInterface;
class MobilityModelInterface;
class ThermoelectricPower;
class PyroPolarization;



//! The base class for all drift-diffusion related semiconductor models
/*!
 * \note { it is not yet possible to have twice the same recombination model.
 * Trying to add to identical models will result in a memory leak. This will
 * be corrected in future. }
 */
class DriftDiffusionProperties : public PhysicalModel, public Variable
{
    
  public:

    //! A nested class that holds all point data
    class PointData
    {
      public:

        //! Constructor
        PointData(void);

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


        double old_electric_potential;
        double old_fermi_e;
        double old_fermi_h;


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
        //std::vector<double> electron_conductivity_derivatives;

        //! The hole conductivity
        double hole_conductivity;

        //! The derivatives of the hole conductivity
        //std::vector<double> hole_conductivity_derivatives;

        //! The net electron recombination rate
        double electron_recombination_rate;

        //! The derivatives of the net electron recombination rate
        std::vector<double> electron_recombination_rate_derivatives;

        //! The net hole recombination rate
        double hole_recombination_rate;

        //! The derivatives of the net hole recombination rate
        std::vector<double> hole_recombination_rate_derivatives;



      private:

        //! We let DriftDiffusionProperties talk to it directly
        friend class DriftDiffusionProperties;
    };


       
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


    //! Lock the parameters
    /*!
     * Call this method before reinit()
     *
     * You can provide a pointer to a PointData object created outside
     * or let it create one internally
     */
    void lock(PointData* pd = NULL);


    //! Unlock the parameters
    /*!
     * Call this after all calculations on the element hav been done
     */
    void unlock(void);
    
    
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


    //! Tells if this model is for a dielectric
    bool is_dielectric(void) const;

    
    //! Setup the band edge data
    /*!
     * This implementation calculates the effective density of states
     * and sets the band edges.
     */
    void setup_band_edges(void);


    //! Set the polarization vector
    void set_polarization(const RealVectorValue& polarization);


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

    // Set the old potentials
    void set_old_potentials(double potential, double Ef_e = 0.0, double Ef_h = 0.0);

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
    void set_electric_field(const RealGradient& E);
    
    //! Set the gradient of the electron electr-chemical potential
    void set_grad_fermi_e(const RealGradient& grad_Fe);
    
    //! Set the gradient of the hole electr-chemical potential
    void set_grad_fermi_h(const RealGradient& grad_Fh);

    //! Get the electric field
    const RealGradient& get_electric_field(void) const;

    //! Get the gradient of the electron electr-chemical potential
    const RealGradient& get_grad_fermi_e(void) const;

    //! Get the gradient of the hole electr-chemical potential
    const RealGradient& get_grad_fermi_h(void) const;


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
      { return _pd->electron_density; };
     
    

    //! Get the electron density derivative
    /*!
     * \return the electron density derivative with respect to the
     * electric potential
     */
    double get_electron_density_derivative(void) const
      { return _pd->electron_density_derivative; };
    
    
    //! Get the hole density
    /*!
     * Get the hole density as calculated by \c calculate_all(...)
     * 
     * \return the hole density
     */
    double get_hole_density(void) const
      { return _pd->hole_density; };
     
    
    //! Get the ehole density derivative
    /*!
     * \return the hole density derivative with respect to the
     * electric potential
     */
    double get_hole_density_derivative(void) const
      { return _pd->hole_density_derivative; };
    

    //! Get the ionized donor density
    /*!
     * \return the ionized donor density
     */
    double get_ionized_donor_density(void) const
      { return _pd->ionized_donor_density; };
     
    
    //! Get the ionized donor density derivative
    /*!
     * \return the ionized donor density derivative with respect to the
     * electric potential
     */
    double get_ionized_donor_density_derivative(void) const
      { return _pd->ionized_donor_density_derivative; };
        
    
    //! Get the ionized acceptor density
    /*!
     * \return the ionized acceptor density
     */
    double get_ionized_acceptor_density(void) const
      { return _pd->ionized_acceptor_density; };
     
    
    //! Get the ionized acceptor density derivative
    /*!
     * \return the ionized acceptor density derivative with respect to the
     * electric potential
     */
    double get_ionized_acceptor_density_derivative(void) const
      { return _pd->ionized_acceptor_density_derivative; };
    
    
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
      { return _pd->electron_recombination_rate; };
      
    
    //! Get the net electron recombination rate derivatives
    /*!
     * Get \f$\frac{\partial R_{net}}{\partial\varphi}\f$
     *
     * \return the derivatives as a const vector reference
     */
    const std::vector<double>&
      get_net_electron_recombination_rate_derivatives(void) const
        { return _pd->electron_recombination_rate_derivatives; };
    
    
    //! Get the net hole recombination rate
    /*!
     * Get \f$R_{net} = R - G\f$ as
     * calculated by \c calculate_all(...)
     *
     */
    double get_net_hole_recombination_rate(void) const
      { return _pd->hole_recombination_rate; };
      
    
    //! Get the net hole recombination rate derivatives
    /*!
     * Get \f$\frac{\partial R_{net}}{\partial\varphi}\f$
     *
     * \return the derivatives as a const vector reference
     */
    const std::vector<double>&
      get_net_hole_recombination_rate_derivatives(void) const
        { return _pd->hole_recombination_rate_derivatives; };

    
    //! Get the total electric polarization
    /*!
     * The total electric polarization \b P is the sum of the
     * pyroelectric and piezoelectric polarization
     */
    const RealVectorValue& get_total_polarization(void) const
      { return _polarization; };

    
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
      { return _pd->electron_mobility; };

    //! Get the electron mobility model
    MobilityModelInterface* get_electron_mobility_model(void)
      { return _electron_mobility; };
      
    
    //! Get the hole mobility
    /*!
     * \return the hole mobility
     */
    double get_hole_mobility(void) const
      { return _pd->hole_mobility; };

    //! Get the hole mobility model
    MobilityModelInterface* get_hole_mobility_model(void)
      { return _hole_mobility; };
     

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

    
    //! Get the lowest conduction band edge
    double get_conduction_band_edge(void) const
      { return conduction_band.band_edge; };

    //! Get all conduction bands
    const std::vector<double>& get_conduction_bands(void)
      { return conduction_band.band_edges; };

    
    //! Get the highest valence band edge
    double get_valence_band_edge(void) const
      { return valence_band.band_edge; };

    //! Get all valence bands
    const std::vector<double>& get_valence_bands(void)
      { return valence_band.band_edges; };

    
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


    //! Returns the thermoelectric power for electrons
    double get_electron_thermoelectric_power() const;

    //! Returns the thermoelectric power for holes 
    double get_hole_thermoelectric_power() const;

    //! Computes the electron and hole thermoelectric powers
    void compute_thermoelectric_powers(void);

    //! Computes the electron and hole thermoelectric power derivatives
    void compute_thermoelectric_power_gradient(void);

    //!provides holes thermoelectric power [V/K]
    RealGradient get_electron_thermoelectric_power_gradient(void) const;

    //!provides holes thermoelectric power [V/K]
    RealGradient get_hole_thermoelectric_power_gradient(void) const;
  
    //! Get the all nodal temperatures for a given element
    std::vector<double>& get_temperature_at_nodes(void);

    //! Get the electric potential
    double get_electric_potential(void) const;


    //! Get the electron electro-chemical potential
    double get_electron_electro_chemical_potential(void) const
      { return _pd->fermi_e; };


    //! Get the hole electro-chemical potential
    double get_hole_electro_chemical_potential(void) const
      { return _pd->fermi_h; };

    double get_old_phi(void) const
      { return _pd->old_electric_potential; };

    double get_old_fermie(void) const
      { return _pd->old_fermi_e; };

    double get_old_fermih(void) const
      { return _pd->old_fermi_h; };


    //! Set the flag for equilibrium calculation
    void set_driftdiffusion(DriftDiffusion* dd);


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

      //! All the band energies
      std::vector<double> band_edges;
    };


    //! The possible variables
    enum Variables
    {
      UNKNOWN = 0,
      RELAXPOLARIZ //!< Relaxation factor for the polarization
    };

    //! The empty constructor.
    DriftDiffusionProperties(void);


    /*! \copydoc PhysicalModel::read_database() */
    virtual void read_database(void);

    
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


    //! Set the value of a variable
    virtual void set_variable_value(double value, ID id = 0);

    //! Get the value of a variable
    virtual double get_variable_value(ID id = 0);
    

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


    //! \copydoc PhysicalModel::do_print_info(void)
    virtual void do_print_info(void);
    

    //! Get the point data structure
    PointData& get_pd(void);


    //! Get the strain as writable reference
    Tensor2Sym& get_strain(void);
    
    
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


    //! Get the temperature interface
    TemperatureInterface& get_temperature_interface(void);

    //! Get the strain interface
    StrainInterface& get_strain_interface(void);


    //! Tells if we should assume inhomogeneous band parameters
    bool is_inhomogeneous(void) const;


    //! Tells if we are doing equilibrium calculation
    bool has_solution(void) const;


    //! Tells if we should use a predictor for quantum densities
    bool use_predictor(void) const;



  private:

    //! The interface to the lattice temperature simulation
    TemperatureInterface _lattice_temp;


    //! The interface to a strain simulation
    StrainInterface _strain_if;


    //! \c true if we should assume inhomogeneous band parameters
    bool _is_inhomogeneous;

    //! \c true if we should use a predictor for quantum densities
    bool _use_predictor;

    //! The simulation this model is used for
    DriftDiffusion* _driftdiffusion;


    //! The point-wise data
    PointData* _pd;

       
    //! Electron thermoelectric power gradient
    RealGradient _eTEpowerGrad;

    //! Hole thermoelectric power gradient
    RealGradient _hTEpowerGrad; 

    //! Electron thermoelectric power
    double _eTEpower;

    //! Hole thermoelectric power
    double _hTEpower;

    //! The electric field
    RealGradient _electric_field;

    //! The gradient of the electron chemical-potential
    RealGradient _grad_fermi_e;

    //! The gradient of the hole chemical-potential
    RealGradient _grad_fermi_h;


    //! The pyropolarization
    PyroPolarization* _pyropolarization;
    
    //! The total electric polarization
    RealVectorValue _polarization;


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
    Point _coord;

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


    //! The thermoelectric power
    ThermoelectricPower* _thermoelectric_power;


    //! The lattice temperature in eV (\f$= k_B T_{lat} / e\f$)
    double _lattice_vt;


    //! The nodal lattice temperature
    std::vector<double> _nodal_lattice_vt;


    //! The constant factor to calculate the effective density of states
    /*!
     * The spin degeneracy has to be included in the effective mass
     */
    static const double _DOS_factor;


    //! True if this is a dielectric
    bool _is_dielectric;


    //! The electrons 
    ParticleDensity _electrons;


    //! The holes
    ParticleDensity _holes;


    //! The relaxation factor for the polarization
    double _relax_polariz;


    //! Whether to use generalized Einstein relation
    bool _generalized_einstein_relation;


    //! The stack for point data locking
    /*!
     * The second value tells whether to delete the data after use
     */
    std::stack<std::pair<PointData*, bool> > _pd_stack;

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
bool
DriftDiffusionProperties::is_dielectric(void) const
{
  return _is_dielectric;
}


 
inline
void
DriftDiffusionProperties::set_coordinates(const Point& p)
{
  _coord = p;
}

inline
void
DriftDiffusionProperties::set_carrier_temperatures(double T_e, double T_h)
{
  _pd->electron_vt = T_e;
  _pd->hole_vt = T_h;
}



inline
double
DriftDiffusionProperties::get_electric_potential(void) const
{
  return _pd->electric_potential;
}



inline
void
DriftDiffusionProperties::set_potentials(double potential, double Ef_e,
    double Ef_h)
{
  _pd->electric_potential = potential;
  _pd->fermi_e = Ef_e;
  _pd->fermi_h = Ef_h;
}


inline
void
DriftDiffusionProperties::set_old_potentials(double potential, double Ef_e,
    double Ef_h)
{
  _pd->old_electric_potential = potential;
  _pd->old_fermi_e = Ef_e;
  _pd->old_fermi_h = Ef_h;
}


inline
void
DriftDiffusionProperties::set_electric_field(const RealGradient& E)
{
  _electric_field = E;
}

inline
void
DriftDiffusionProperties::set_grad_fermi_e(const RealGradient& grad_Fe)
{
  _grad_fermi_e = grad_Fe;
}


inline
void
DriftDiffusionProperties::set_grad_fermi_h(const RealGradient& grad_Fh)
{
  _grad_fermi_h = grad_Fh;
}



inline
const RealGradient&
DriftDiffusionProperties::get_electric_field(void) const
{
  return _electric_field;
}



inline
const RealGradient&
DriftDiffusionProperties::get_grad_fermi_e(void) const
{
  return _grad_fermi_e;
}



inline
const RealGradient&
DriftDiffusionProperties::get_grad_fermi_h(void) const
{
  return _grad_fermi_h;
}




inline
void
DriftDiffusionProperties::set_densities(double n, double p)
{
  _pd->electron_density = n;
  _pd->hole_density = p;
  _pd->electron_density_derivative = n / _pd->electron_vt;
  _pd->hole_density = p;
  _pd->hole_density_derivative = -p / _pd->hole_vt;
}



inline
void
DriftDiffusionProperties::set_lattice_temperature(double T)
{
  _lattice_vt = T * Constants::k_B;
}



inline
double
DriftDiffusionProperties::get_lattice_temperature(void) const
{
  return _lattice_vt;
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
  return _coord;
}
 

inline
double
DriftDiffusionProperties::get_charge_density(void) const
{
  return _pd->hole_density - _pd->electron_density +
    _pd->ionized_donor_density - _pd->ionized_acceptor_density;
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
Tensor2Sym&
DriftDiffusionProperties::get_strain(void) 
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




inline
void
DriftDiffusionProperties::setup_band_edges(void)
{
  double kT = _lattice_vt;
  _pd->electron_vt = _pd->hole_vt = kT;
  
  BandProperties& cb = conduction_band;
  BandProperties& vb = valence_band;

  cb.effective_DOS =
    get_DOS_factor() * std::pow(kT * cb.effective_mass, 1.5);

  vb.effective_DOS =
    get_DOS_factor() * std::pow(kT * vb.effective_mass, 1.5);
}

inline
void
DriftDiffusionProperties::set_polarization(const RealVectorValue& polarization)
{
  _polarization += _relax_polariz * polarization;
}


inline
int
DriftDiffusionProperties::get_number_of_recombination_models(void) const
{
  return _recombination_models.size();
}


inline
double
DriftDiffusionProperties::get_electron_thermoelectric_power(void) const
{
  return _eTEpower;
}


inline
double
DriftDiffusionProperties::get_hole_thermoelectric_power(void) const
{
  return _hTEpower;
}

inline
RealGradient 
DriftDiffusionProperties::get_electron_thermoelectric_power_gradient(void) const
{

  return  _eTEpowerGrad;

} 

inline
RealGradient 
DriftDiffusionProperties::get_hole_thermoelectric_power_gradient(void) const
{

  return  _hTEpowerGrad;

}
  
  

inline
PhysicalModelInterface*
DriftDiffusionProperties::create_new(void) const
{
  return new DriftDiffusionProperties();
}



inline
DriftDiffusionProperties::PointData&
DriftDiffusionProperties::get_pd(void)
{
  return *_pd;
}


inline
TemperatureInterface&
DriftDiffusionProperties::get_temperature_interface(void)
{
  return _lattice_temp;
}

inline
StrainInterface&
DriftDiffusionProperties::get_strain_interface(void)
{
  return _strain_if;
}


inline
bool
DriftDiffusionProperties::is_inhomogeneous(void) const
{
  return _is_inhomogeneous;
}


inline
bool
DriftDiffusionProperties::has_solution(void) const
{
  return _driftdiffusion->is_solved();
}

inline
bool
DriftDiffusionProperties::use_predictor(void) const
{
  return _use_predictor;
}


inline
void
DriftDiffusionProperties::set_driftdiffusion(DriftDiffusion* dd)
{
  _driftdiffusion = dd;
}


#endif /* _DRIFTDIFFUSIONPROPERTIES_H_ */
