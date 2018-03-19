// $Id$

#ifndef _DRIFTDIFFUSIONPROPERTIES_H_
#define _DRIFTDIFFUSIONPROPERTIES_H_



#include "PhysicalModel.h"

#include "TemperatureInterface.h"
#include "BandProperties.h"
#include "SimulationOptions.h"
#include "DriftDiffusionDefs.h"
#include "TiberCad.h"
#include "Constants.h"
#include "TypeDefs.h"

#include "tensor.h"

#include "tensor_value.h"
#include "vector_value.h"
#include "point.h"


#include <vector>
#include <set>
#include <map>


// forward declarations
//class Elem;
class Dopant;
class Trap;
class SimulationInterface;
class DriftDiffusion;
class RecombinationModelInterface;
class ParticleDensity;


//! The base class for all drift-diffusion related semiconductor models
class DriftDiffusionProperties : public PhysicalModel
{

  private:

    //! An iterator for the recombination models
    typedef std::multimap<ID, RecombinationModelInterface*>::iterator
      recomb_iterator;

    //! A const iterator for the recombination models
    typedef std::multimap<ID, RecombinationModelInterface*>::const_iterator
      const_recomb_iterator;

  public:

    class RecombinationModelIterator : public recomb_iterator
    {
      public:
        RecombinationModelIterator(const recomb_iterator& it) :
          recomb_iterator(it) {}

        RecombinationModelIterator(const RecombinationModelIterator& it) :
          recomb_iterator(it) {}

        ID id(void) { return((*this)->first); }

        RecombinationModelInterface* operator*(void) { return((*this)->second); }

    };

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

        //! \f$\gamma_n\f$
        double gamma_n;

        //! \f$\gamma_n\f$
        double gamma_p;

        //! The ionized donor density
        double ionized_donor_density;

        //! The ionized donor density derivative
        double ionized_donor_density_derivative;

        //! The ionized acceptor density
        double ionized_acceptor_density;

        //! The ionized acceptor density derivative
        double ionized_acceptor_density_derivative;


        //! The trapped electron density
        double ionized_electron_traps;

        //! The trapped hole density
        double ionized_hole_traps;

        //! The trapped density derivatives d/dfermi_e, d/dfermi_h
        std::vector<double> ionized_traps_derivative;

        //! The total charge density
        double charge_density;

        //! The electron mobility
        double electron_mobility;

        double electron_mobility_derivative_potential;

        libMesh::RealGradient electron_mobility_derivative_grad_potential;

        //! The electron mobility derivative w.r.t gradient of Fermi level
        libMesh::RealGradient electron_mobility_derivatives;


        //! The hole mobility
        double hole_mobility;

        double hole_mobility_derivative_potential;

        libMesh::RealGradient hole_mobility_derivative_grad_potential;

        //! The hole mobility derivative w.r.t gradient of Fermi level
        libMesh::RealGradient hole_mobility_derivatives;


        //! The electron conductivity
        double electron_conductivity;

        // ! The derivatives of the electron conductivity
        //std::vector<double> electron_conductivity_derivatives;

        //! The hole conductivity
        double hole_conductivity;

        // ! The derivatives of the hole conductivity
        //std::vector<double> hole_conductivity_derivatives;

        //! The net electron recombination rate
        double electron_recombination_rate;

        //! The derivatives of the net electron recombination rate
        /*!
         * order is: d/dn, d/dp, d/dfermi_e, d/dfermi_h
         */
        std::vector<double> electron_recombination_rate_derivatives;

        //! The net hole recombination rate
        double hole_recombination_rate;

        //! The derivatives of the net hole recombination rate
        /*!
         * order is: d/dn, d/dp, d/dfermi_e, d/dfermi_h
         */
        std::vector<double> hole_recombination_rate_derivatives;

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
        const Material* mat,
        const ModelOptions& options = ModelOptions());




    // ! Lock the parameters
    /* !
     * Call this method before reinit()
     *
     * You can provide a pointer to a PointData object created outside
     * or let it create one internally
     */
    //void lock(PointData* pd = NULL);


    // ! Unlock the parameters
    /* !
     * Call this after all calculations on the element hav been done
     */
    //void unlock(void);


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

    //! Get or set dielectric flag
    bool& is_dielectric(void);

    //! Set the coordinates
    void set_coordinates(const libMesh::Point& p);

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
    void set_electric_field(const libMesh::RealGradient& E);
    void set_old_electric_field(const libMesh::RealGradient& E);


    //! Set the gradient of the electron electr-chemical potential
    void set_grad_fermi_e(const libMesh::RealGradient& grad_Fe);

    //! Set the gradient of the hole electr-chemical potential
    void set_grad_fermi_h(const libMesh::RealGradient& grad_Fh);

    //! Get the electric field
    const libMesh::RealGradient& get_electric_field(void) const;
    const libMesh::RealGradient& get_old_electric_field(void) const;


    //! Get the gradient of the electron electr-chemical potential
    const libMesh::RealGradient& get_grad_fermi_e(void) const;

    //! Get the gradient of the hole electr-chemical potential
    const libMesh::RealGradient& get_grad_fermi_h(void) const;


    //! Get the element we are currently working on
    const libMesh::Elem* get_element(void) const;

    //! Get the coordinates of the point we are currently working on
    const libMesh::Point& get_coordinates(void) const;

    //! Set the strain
    void set_strain(const Tensor2Sym& strain);

    //! Get the strain
    const Tensor2Sym& get_strain(void) const;


    //! Set the lattice temperature (in K)
    void set_lattice_temperature(double T);

    //! Set the lattice temperature (in eV)
    void set_lattice_vt(double vt);

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


    //! Calculated trapped particle densities
    void calculate_traps(void);


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


    //! Get the electron gamma \f$\gamma_n\f$
    /*!
     * \f$\gamma_n = n/n_{cl}\f$
     */
    double get_electron_gamma(void) const
    { return _pd->gamma_n; };


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


    //! Get the hole gamma \f$\gamma_p\f$
    /*!
     * \f$\gamma_p = p/p_{cl}\f$
     */
    double get_hole_gamma(void) const
    { return _pd->gamma_p; };


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


    //! Get the ionized electron trap density
    /*!
     * \return the ionized electron trap density
     */
    double get_ionized_electron_traps(void) const
      { return _pd->ionized_electron_traps; };


    //! Get the ionized hole trap density
    /*!
     * \return the ionized hole trap density
     */
    double get_ionized_hole_traps(void) const
      { return _pd->ionized_hole_traps; };


    //! Get the total charge density
    /*!
     *
     * \f$ \rho = p - n + N_D^+ - N_A^- \f$
     *
     * \return the total charge density \f$\rho\f$
     *
     * \note
     * The charge density is returned in units of the elementary charge
     * \f$e\f$, \em not in Coulomb!
     */
    double get_charge_density(void) const;


    //! Get the derivatives of the total charge density
    /*!
     * The derivatives are given w.r.t. the two quasi Fermi levels
     */
    void get_charge_density_derivatives(double derivatives[2]) const;


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


    // ! Get the total electric polarization
    /* !
     * The total electric polarization \b P is the sum of the
     * pyroelectric and piezoelectric polarization
     */
    //const libMesh::RealVectorValue& get_total_polarization(void) const
    //  { return _polarization; };


    //! Get the relative permittivity tensor
    const libMesh::RealTensor& get_relative_permittivity(void) const
      { return _permittivity; };


    //! Get the electron conductivity
    /*!
     * \return the electron conductivity \f$\sigma_n = \mu_n n\f$
     */
    double get_electron_conductivity(void) const
      { return _pd->electron_conductivity; };


    //! Get the hole conductivity
    /*!
     * \return the hole conductivity \f$\sigma_p = \mu_p p\f$
     */
    double get_hole_conductivity(void) const
      { return _pd->hole_conductivity; };


    //! Get the electron mobility
    /*!
     * \return the electron mobility
     */
    double get_electron_mobility(void) const
      { return _pd->electron_mobility; };


    //! Get the hole mobility
    /*!
     * \return the hole mobility
     */
    double get_hole_mobility(void) const
      { return _pd->hole_mobility; };


    //! Get the electron mobility derivative w.r.t. electric potential
    double get_electron_mobility_derivative_potential(void) const
      { return _pd->electron_mobility_derivative_potential; };

    //! Get the hole mobility derivative w.r.t. electric potential
    double get_hole_mobility_derivative_potential(void) const
      { return _pd->hole_mobility_derivative_potential; };

    //! Get the electron mobility derivative w.r.t. gradient of the electric potential
    void get_electron_mobility_derivative_grad_potential(libMesh::RealGradient& dmu) const
      { dmu = _pd->electron_mobility_derivative_grad_potential; };

    //! Get the hole mobility derivative w.r.t. gradient of the electric potential
    void get_hole_mobility_derivative_grad_potential(libMesh::RealGradient& dmu) const
      { dmu = _pd->hole_mobility_derivative_grad_potential; };

    //! Get the electron mobility derivative w.r.t. gradient of the fermi potential
    void get_electron_mobility_derivative_grad_fermi(libMesh::RealGradient& dmu) const
      { dmu = _pd->electron_mobility_derivatives; };

    //! Get the hole mobility derivative w.r.t. gradient of the fermi potential
    void get_hole_mobility_derivative_grad_fermi(libMesh::RealGradient& dmu) const
      { dmu = _pd->hole_mobility_derivatives; };

    //! Get the square of the intrinsic density
    double get_intrinsic_density_squared(void) const
      { return _intrinsic_density * _intrinsic_density; };


    //! Get the intrinsic density
    double get_intrinsic_density(void) const
      { return _intrinsic_density; };


    //! Get the equilibrium electron density
    double get_equilibrium_electron_density(void) const
      { return _equilibrium_n; };


    //! Get the equilibrium hole density
    double get_equilibrium_hole_density(void) const
      { return _equilibrium_p; };


    //! Get equilibrium fermi level
    double get_equilibrium_fermi_level(void) const
      { return _equilibrium_fermi_level; };


    //! Get the lowest conduction band edge
    double get_conduction_band_edge(void) const
      { return _conduction_band->get_band_edge(); };

    //! Get all conduction bands
    void get_conduction_bands(std::vector<double>& bands)
      { _conduction_band->get_bands(bands); };


    //! Get the highest valence band edge
    double get_valence_band_edge(void) const
      { return _valence_band->get_band_edge(); };

    //! Get all valence bands
    void get_valence_bands(std::vector<double>& bands)
      { _valence_band->get_bands(bands); };


    //! Get the band gap
    double get_band_gap(void) const
      { return _conduction_band->get_band_edge() - _valence_band->get_band_edge(); };

    //void get_net_recombination_rates(std::vector<double>& rates);
    //! Get the iterator to the first recombination model
    RecombinationModelIterator recombination_models_begin(void);


    //! Get the past-the-end iterator for the recombination models
    RecombinationModelIterator recombination_models_end(void);

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
    std::pair<double, double> get_net_recombination_rate(ID id);


    //! Clear all recombination rates
    void clear_recombination(void);


    //! Returns the number of recombination models
    int get_number_of_recombination_models(void) const;


    //! Returns the thermoelectric power for electrons
    //double get_electron_thermoelectric_power() const;

    //! Returns the thermoelectric power for holes
    //double get_hole_thermoelectric_power() const;

    //! Computes the electron and hole thermoelectric powers
    //void compute_thermoelectric_powers(void);

    //! Computes the electron and hole thermoelectric power derivatives
    //void compute_thermoelectric_power_gradient(void);


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


    //! Get the point data
    const PointData& get_point_data(void) const
        { return get_pd(); }


    //! Create the recombination models
    void create_recombination_models(void);


    //! Get the conduction band properties
    const BandProperties& get_conduction_band(void) const
      { return *_conduction_band; };


    //! Get the valence band properties
    const BandProperties& get_valence_band(void) const
      { return *_valence_band; };


    //! Get the conduction band properties
    BandProperties& get_conduction_band(void)
      { return *_conduction_band; };


    //! Get the valence band properties
    BandProperties& get_valence_band(void)
      { return *_valence_band; };



    //! Tells if we are doing equilibrium calculation
    bool has_solution(void) const;

    //! Get the temperature interface
    TemperatureInterface& get_temperature_interface(void);



  protected:



    //! The empty constructor.
    DriftDiffusionProperties(const ModelOptions& options);


    //! Initialize this model
    /*!
     * This reads the database and calls init for all submodels
     * A derived class which reimplements this method has to call
     * explicitly the one of this class!
     */
    virtual void do_init(void);


    //! Create some of the submodels
    virtual void prepare_submodels(void);


    //! Set the element we are currently working on
    void set_element(const Elem* elem);


    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     * This method can be overriden by derived classes.
     */
    virtual void prepare_element_data(void) {};



    //! Get the point data structure
    PointData& get_pd(void) const;


    //! Get the strain as writable reference
    Tensor2Sym& get_strain(void);


    //! Get the strain interface
    //StrainInterface& get_strain_interface(void);


    //! Tells if we should assume inhomogeneous band parameters
    bool is_inhomogeneous(void) const;






    //! Set the relative permittivity tensor
    void set_relative_permittivity(const libMesh::RealTensor& eps)
      { _permittivity = eps; };





    //! Set the equilibrium Fermi level
    void set_equilibrium_fermi_level(double Ef)
      { _equilibrium_fermi_level = Ef; }


    //! Set the intrinsic Fermi level
    /*!
     * Calculates the associated equilibrium densities
     */
    void set_equilibrium_properties(double Ef);


    //! Get the conduction band properties as pointer
    void set_conduction_band(BandProperties* cb);


    //! Get the valence band properties as pointer
    void set_valence_band(BandProperties* vb);


    //! Set the intrinsic density
    void set_intrinsic_density(double ni)
      { _intrinsic_density = ni; };


    //! Set the equilibrium electron density
    void set_equilibrium_n(double n)
      { _equilibrium_n = n; };


    //! Set the equilibrium hole density
    void set_equilibrium_p(double p)
      { _equilibrium_p = p; };



  private:

    //! The interface to the lattice temperature simulation
    TemperatureInterface _lattice_temp;


    //! The interface to a strain simulation
    //StrainInterface _strain_if;


    //! \c true if we should assume inhomogeneous band parameters
    bool _is_inhomogeneous;





    //! The equilibrium fermi level
    /*!
     * The fermi level such that \f$n=n_0,\,p=p_0\f$
     */
    double _equilibrium_fermi_level;

    //! The intrinsic density
    double _intrinsic_density;

    //! The equilibrium electron density
    double _equilibrium_n;

    //! The equilibrium electron density
    double _equilibrium_p;

    //! The point-wise data
    PointData* _pd;


    //! The electric field
    libMesh::RealGradient _electric_field;
    libMesh::RealGradient _old_electric_field;


    //! The gradient of the electron chemical-potential
    libMesh::RealGradient _grad_fermi_e;

    //! The gradient of the hole chemical-potential
    libMesh::RealGradient _grad_fermi_h;

    //libMesh::RealVectorValue _polarization;

    //! The relative permittivity tensor
    libMesh::RealTensor _permittivity;


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



    //! The element we are currently working on
    const libMesh::Elem* _elem;

    //! The coordinates of the point we are working on
    libMesh::Point _coord;


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
    BandProperties* _conduction_band;

    //! The conduction band properties
    /*!
     * Band properties are assumed to be elemental data, \em not nodal data
     */
    BandProperties* _valence_band;


    //! The electron traps
    std::set<Trap*> _etraps;

    //! The hole traps
    std::set<Trap*> _htraps;


    //! The recombination models
    std::multimap<ID, RecombinationModelInterface*> _recombination_models;


    //! The thermoelectric power
    //ThermoelectricPower* _thermoelectric_power;


    //! The polarization models
    //std::vector<PolarizationModel*> _pm;


    //! The lattice temperature in eV (\f$= k_B T_{lat} / e\f$)
    double _lattice_vt;


    //! True if this is a dielectric
    bool _is_dielectric;










    //! The relaxation factor for the polarization
    //double _relax_polariz;


};


//
// inline members
//





inline
bool
DriftDiffusionProperties::is_dielectric(void) const
{
  return _is_dielectric;
}



inline
bool&
DriftDiffusionProperties::is_dielectric(void)
{
  return _is_dielectric;
}


inline
void
DriftDiffusionProperties::set_coordinates(const libMesh::Point& p)
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
DriftDiffusionProperties::set_electric_field(const libMesh::RealGradient& E)
{
  _electric_field = E;
}








inline
void
DriftDiffusionProperties::set_old_electric_field(const libMesh::RealGradient& E)
{
  _old_electric_field = E;
}

inline
void
DriftDiffusionProperties::set_grad_fermi_e(const libMesh::RealGradient& grad_Fe)
{
  _grad_fermi_e = grad_Fe;
}


inline
void
DriftDiffusionProperties::set_grad_fermi_h(const libMesh::RealGradient& grad_Fh)
{
  _grad_fermi_h = grad_Fh;
}



inline
const libMesh::RealGradient&
DriftDiffusionProperties::get_electric_field(void) const
{
  return _electric_field;
}

inline
const libMesh::RealGradient&
DriftDiffusionProperties::get_old_electric_field(void) const
{
  return _old_electric_field;
}








inline
const libMesh::RealGradient&
DriftDiffusionProperties::get_grad_fermi_e(void) const
{
  return _grad_fermi_e;
}



inline
const libMesh::RealGradient&
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
void
DriftDiffusionProperties::set_lattice_vt(double vt)
{
  _lattice_vt = vt;
}


inline
double
DriftDiffusionProperties::get_lattice_temperature(void) const
{
  return _lattice_vt;
}

inline
void
DriftDiffusionProperties::set_element(const Elem* elem)
{
  _elem = elem;
}

inline
const libMesh::Elem*
DriftDiffusionProperties::get_element(void) const
{
  return _elem;
}

inline
const libMesh::Point&
DriftDiffusionProperties::get_coordinates(void) const
{
  return _coord;
}


inline
double
DriftDiffusionProperties::get_charge_density(void) const
{
  long double dens = static_cast<long double>(_pd->hole_density) -
      static_cast<long double>(_pd->electron_density) +
      static_cast<long double>(_pd->ionized_donor_density) -
      static_cast<long double>(_pd->ionized_acceptor_density) +
      static_cast<long double>(_pd->ionized_electron_traps) +
      static_cast<long double>(_pd->ionized_hole_traps);

  return static_cast<double>(dens);
  //return _pd->hole_density - _pd->electron_density +
  //  _pd->ionized_donor_density - _pd->ionized_acceptor_density
  //  + _pd->ionized_electron_traps + _pd->ionized_hole_traps;
}


inline
void
DriftDiffusionProperties::get_charge_density_derivatives(
    double derivatives[2]) const
{
  long double der0 = static_cast<long double>(get_electron_density_derivative())
      - static_cast<long double>(get_ionized_donor_density_derivative())
      + static_cast<long double>(_pd->ionized_traps_derivative[0]);

  long double der1 = -static_cast<long double>(get_hole_density_derivative())
      + static_cast<long double>(get_ionized_acceptor_density_derivative())
      + static_cast<long double>(_pd->ionized_traps_derivative[1]);

  derivatives[0] = static_cast<double>(der0);
  derivatives[1] = static_cast<double>(der1);

  //derivatives[0] = get_electron_density_derivative()
  //  - get_ionized_donor_density_derivative()
  //  - _pd->ionized_electron_traps_derivative;
  //derivatives[1] = -get_hole_density_derivative()
  //  + get_ionized_acceptor_density_derivative()
  //  - _pd->ionized_hole_traps_derivative;
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
int
DriftDiffusionProperties::get_number_of_recombination_models(void) const
{
  return _recombination_models.size();
}



inline
DriftDiffusionProperties::PointData&
DriftDiffusionProperties::get_pd(void) const
{
  return *_pd;
}


inline
TemperatureInterface&
DriftDiffusionProperties::get_temperature_interface(void)
{
  return _lattice_temp;
}

//inline
//StrainInterface&
//DriftDiffusionProperties::get_strain_interface(void)
//{
//  return _strain_if;
//}


inline
bool
DriftDiffusionProperties::is_inhomogeneous(void) const
{
  return _is_inhomogeneous;
}


inline
DriftDiffusionProperties::RecombinationModelIterator
DriftDiffusionProperties::recombination_models_begin(void)
{
  return(RecombinationModelIterator(_recombination_models.begin()));
}


inline
DriftDiffusionProperties::RecombinationModelIterator
DriftDiffusionProperties::recombination_models_end(void)
{
  return(RecombinationModelIterator(_recombination_models.end()));
}


#endif /* _DRIFTDIFFUSIONPROPERTIES_H_ */
