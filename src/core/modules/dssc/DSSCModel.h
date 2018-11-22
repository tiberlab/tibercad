// $Id$

#ifndef _DSSCMODEL_H_
#define _DSSCMODEL_H_


#ifndef TIBER_MODULE_PREFIX
# define TIBER_MODULE_PREFIX dssc
#endif


#include "PhysicalModel.h"

#include "ParticleDensity.h"
#include "TemperatureInterface.h"
#include "SimulationOptions.h"
#include "TiberCad.h"
#include "Constants.h"
#include "TypeDefs.h"
#include "DSSC.h"

#include "vector_value.h"
#include "point.h"


// forward declarations
class Elem;
class Trap;


//! The base class for DSSC models
class DSSCModel : public PhysicalModel
{

  public:

    //! The thermodynamic equilibrium concentrations
    struct EquilibriumConcentrations
    {
      double n;
      double I;
      double I3;
      double C;
    };


    //! The mobilities
    struct Mobilities
    {
      //! The electron mobility
      double n;

      //! The iodide mobility
      double I;

      //! The triiodide mobility
      double I3;

      //! The cation mobility
      double C;
    };


    //! The empty constructor.
    DSSCModel(const ModelOptions& options);


    //! A default (empty) destructor.
    virtual ~DSSCModel(void) { };


    //! Create a named drift-diffusion model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static DSSCModel* create(const std::string& name,
        const ModelOptions& options = ModelOptions());


    //! Create an instance of this class
    static DSSCModel* create(const ModelOptions& options);


    //! (Re-)Initialize for the given element
    /*!
     * \c reinit() calls \c prepare_element_data() which needs to be
     * implemented in derived classes
     */
    void reinit(const Elem* elem);


    //! Tells if this model is for electrolyte
    bool is_electrolyte(void) const;

    //! Tells if this model is for TiO2
    bool is_TiO2(void) const;


    //! Get the TiO2 porosity
    /*!
     * In principle this would make sense only if
     * is_electrolyte() and is_TiO2() bth return true
     */
    double porosity(void) const;


    //! Set the coordinates
    void set_coordinates(const Point& p);


    //! Set the potentials
    /*!
     * The potentials have to be set before calling any of the
     * \c calculate_xxx methods. If the densities are used as variables,
     * the electro-chemical potentials can be left out.
     *
     * \param potential the electric potential
     * \param Ef_n the electron electro-chemical potential
     * \param Ef_I the iodide electro-chemical potential
     * \param Ef_I3 the triiodide electro-chemical potential
     * \param Ef_C the cation electro-chemical potential
     */
    void set_potentials(double electric_potential, double Ef_n = 0.0,
        double Ef_I = 0.0, double Ef_I3 = 0.0, double Ef_C = 0.0);


    //! Set coordinates of the contact under illumination
    void set_x0(double x0);


    //! Set the electric field
    void set_electric_field(const libMesh::RealGradient& E);


    //! Set the gradient of the electron electro-chemical potential
    void set_grad_fermi_n(const libMesh::RealGradient& grad_F);


    //! Set the gradient of the iodide electro-chemical potential
    void set_grad_fermi_I(const libMesh::RealGradient& grad_F);


    //! Set the gradient of the triiodide electro-chemical potential
    void set_grad_fermi_I3(const libMesh::RealGradient& grad_F);


    //! Set the gradient of the cation electro-chemical potential
    void set_grad_fermi_C(const libMesh::RealGradient& grad_F);


    //! Get coordinates of the contact under illumination
    const double get_x0() const;

    
    //! Get the electric field
    const libMesh::RealGradient& get_electric_field(void) const;


    //! Get the gradient of the electron electro-chemical potential
    const libMesh::RealGradient& get_grad_fermi_n(void) const;


    //! Get the gradient of the iodide electro-chemical potential
    const libMesh::RealGradient& get_grad_fermi_I(void) const;


    //! Get the gradient of the triiodide electro-chemical potential
    const libMesh::RealGradient& get_grad_fermi_I3(void) const;


    //! Get the gradient of the cation electro-chemical potential
    const libMesh::RealGradient& get_grad_fermi_C(void) const;


    //! Get the electric potential
    double get_electric_potential(void) const;


    //! Get the element we are currently working on
    const Elem* get_element(void) const;


    //! Get the coordinates of the point we are currently working on
    const Point& get_coordinates(void) const;


    //! Set the lattice temperature (in K)
    void set_lattice_temperature(double T);


    //! Get the lattice temperature (in units of eV)
    double get_lattice_temperature(void) const;


    //! Calculate electron and hole densities and derivatives
    /*!
     * This method calculates electron and hole densities and their
     * derivatives with respect to the electric and the
     * electro-chemical potentials.
     */
    void calculate_densities(void);


    //! Calculate net recombination rate and derivative
    void calculate_net_recombination_rate(void);


    //! Calculate trap density
    void calculate_traps(void);

   
    //! Calculate dark background for traps (only real traps) 
    void calculate_equilibrium_traps(void);


    //! Get the electron density
    double get_density_n(void) const;
    

    //! Get the electron density exponential DOS
    double get_density_n_exp_DOS(void) const;
    

    //! Get the electron density derivative
    double get_density_derivative_n(void) const;


    //! Get the electron density derivative exponential DOS
    double get_density_derivative_n_exp_DOS(void) const;


    //! Get the iodide density
    double get_density_I(void) const;


    //! Get the iodide density derivative
    double get_density_derivative_I(void) const;


    //! Get the triiodide density
    double get_density_I3(void) const;


    //! Get the triiodide density derivative
    double get_density_derivative_I3(void) const;


    //! Get the cation density
    double get_density_C(void) const;


    //! Get the cation density derivative
    double get_density_derivative_C(void) const;


    //! Get the ionized dye density
    double get_ionized_dye_density(void) const;


    //! Get the total charge density
    double get_charge_density(void) const;


    //! Get the net recombination rate
    double get_net_recombination_rate(void) const;


    //! Get the recombination rate
    double get_recombination_rate(void) const;


    //! Get the net recombination rate derivatives
    const std::vector<double>&
      get_net_recombination_rate_derivatives(void) const;


    //! Compute the electron density mobility for electrons
    double get_dens_elec_mobility(void) const;


    //! Get electron trapped density
    double get_ionized_electron_traps(void) const; 


    //! Get electron trapped density derivative
    double  get_ionized_electron_traps_derivative(void) const;


    //! Get generation rate
    double get_generation_rate(void) const;


    //! Get the relative permittivity
    /*!
     * \note what should this be in porous materials?
     */
    double get_relative_permittivity(void) const
      { return _permittivity; };


    //! Get the electron mobility
    double get_mobility_n(void) const
      { return _mobility.n; };


    //! Get the iodide mobility
    double get_mobility_I(void) const
      { return _mobility.I; };


    //! Get the triiodide mobility
    double get_mobility_I3(void) const
      { return _mobility.I3; };


    //! Get the cation mobility
    double get_mobility_C(void) const
      { return _mobility.C; };


    //! Get all the nodal temperatures for a given element
    //std::vector<double>& get_temperature_at_nodes(void);


    //! Get the equilibrium concentrations
    const EquilibriumConcentrations& get_equilibrium_concentrations(void) const;


  protected:

    struct PointData
    {
      //! The coordinates of the point we are working on
      Point coordinates;

      double electric_potential;
      double fermi_n;
      double fermi_I;
      double fermi_I3;
      double fermi_C;

      libMesh::RealGradient electric_field;

      //! The gradient of the electron chemical-potential
      libMesh::RealGradient grad_fermi_n;

      //! The gradient of the iodide chemical-potential
      libMesh::RealGradient grad_fermi_I;

      //! The gradient of the triiodide chemical-potential
      libMesh::RealGradient grad_fermi_I3;

      //! The gradient of the cation chemical-potential
      libMesh::RealGradient grad_fermi_C;


      //! The temperature in eV (\f$= k_B T_{lat} / e\f$)
      double kT;

      double density_n;
      double density_I;
      double density_I3;
      double density_C;

      double ionized_dye;


      //! The recombination
      double recombination_rate;


      //! The recombination rate derivatives
      std::vector<double> recombination_rate_derivatives;


      //! The generation
      double generation_rate;

      double ionized_electron_traps;
      double ionized_electron_traps_derivative;
      double ionized_equilibrium_electron_traps;

    };


    /*! \copydoc PhysicalModel::read_database() */
    //virtual void read_database(void) { };
    /*! \copydoc PhysicalModel::read_database() */
    virtual void read_database(void);


    //! Initialize this model
    /*!
     * This reads the database and calls init for all submodels
     * A derived class which reimplements this method has to call
     * explicitly the one of this class!
     */
    virtual void do_init(void);


    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     * This method can be used overiden by derived classes.
     */
    virtual void prepare_element_data(void) {};

    
    // Initialize submodels
    void prepare_submodels(void); 


    //! \copydoc PhysicalModel::create_new()
    virtual PhysicalModel* create_new(void) const;


    //! \copydoc PhysicalModel::copy_from()
    virtual void copy_from(const PhysicalModel* rhs);


    //! \copydoc PhysicalModel::do_print_info(void)
    virtual void do_print_info(void);


    //! Get the temperature interface
    TemperatureInterface& get_temperature_interface(void);


    //! Tells if we are doing equilibrium calculation
    bool has_solution(void) const;


    //! Get the pointwise data
    PointData& get_pd(void);


  private:

    //! The interface to the lattice temperature simulation
    TemperatureInterface _lattice_temp;


    //! The copy constructor is disabled
    DSSCModel(const DSSCModel& rhs);


    //! The assignment operator is disabled
    DSSCModel& operator=(const DSSCModel& rhs);


    //! The element we are currently working on
    const Elem* _elem;


    //! The porosity
    double _porosity;


    //! True if we are in an electrolyte
    bool _is_electrolyte;


    //! True if we are in TiO2
    bool _is_TiO2;


    //! The relative permittivity tensor
    double _permittivity;

    //! The permittivity oxide
    double _perm_ox;

    //! The permittivity electrolyte
    double _perm_elec;
    
    //! The nodal lattice temperature
    //std::vector<double> _nodal_lattice_vt;
 
    //! The electron traps
    std::vector<Trap*> _etraps;

    //! The electrons
    ParticleDensity _electrons;


    //! The iodide
    ParticleDensity _iodide;


    //! The triiodide
    ParticleDensity _triiodide;


    //! The cation
    ParticleDensity _cation;


    //! The thermodynamic equilibrium densities
    EquilibriumConcentrations _eq_conc;


    //! The mobilities
    Mobilities _mobility;


    //! The pointwise data
    PointData _pd;


    //! Electron relaxation rate (for my friends, Err)
    double _ke;


    //! Rate of sensitizer regeneration (for my friends, Rosr)
    double _k3;


    //! exponential trap density
    double _trap_DOS;


    //! trap exponent
    double _exp_trap;


    //! CB density
    double _CB_DOS;


    //! The generation rate
    double _generation;


    //! The generation model
    std::vector<SimulationInterface*> _generation_model;

    //! The solution ID of the generation models variable
    std::vector<ID> _gen_id;


    //! coordinate of the contact under illumination
    double _x0;


    //! absorption coefficient
    double _alpha;


    //! absorption coefficient
    double _alpha2;


    //! absorption coefficient
    double _deltaG;


    //! beta factor
    double _beta;


};


//
// inline members
//


inline
DSSCModel*
DSSCModel::create(const std::string& name,
    const ModelOptions& options)
{
  //return dynamic_cast<DSSCModel*>(
  //    PhysicalModel::create("dscbulk_" + name, options));
  return new DSSCModel(options);
}


inline
DSSCModel*
DSSCModel::create(const ModelOptions& options)
{
  return new DSSCModel(options);
}


inline
bool
DSSCModel::is_electrolyte(void) const
{
  return _is_electrolyte;
}


inline
bool
DSSCModel::is_TiO2(void) const
{
  return _is_TiO2;
}


inline
double
DSSCModel::porosity(void) const
{
  return _porosity;
}


inline
void
DSSCModel::set_coordinates(const Point& p)
{
  _pd.coordinates = p;
}


inline
double
DSSCModel::get_electric_potential(void) const
{
  return _pd.electric_potential;
}


inline
void
DSSCModel::set_potentials(double electric_potential, double Ef_n,
    double Ef_I, double Ef_I3, double Ef_C)
{
  _pd.electric_potential = electric_potential;
  _pd.fermi_n = Ef_n;
  _pd.fermi_I = Ef_I;
  _pd.fermi_I3 = Ef_I3;
  _pd.fermi_C = Ef_C;
}


inline
void
DSSCModel::set_x0(double x0)
{
  _x0 = x0;
}


inline
void
DSSCModel::set_electric_field(const libMesh::RealGradient& E)
{
  _pd.electric_field = E;
}


inline
void
DSSCModel::set_grad_fermi_n(const libMesh::RealGradient& grad_F)
{
  _pd.grad_fermi_n = grad_F;
}


inline
void
DSSCModel::set_grad_fermi_I(const libMesh::RealGradient& grad_F)
{
  _pd.grad_fermi_I = grad_F;
}


inline
void
DSSCModel::set_grad_fermi_I3(const libMesh::RealGradient& grad_F)
{
  _pd.grad_fermi_I3 = grad_F;
}


inline
void
DSSCModel::set_grad_fermi_C(const libMesh::RealGradient& grad_F)
{
  _pd.grad_fermi_C = grad_F;
}


inline
const double 
DSSCModel::get_x0(void) const
{
  return _x0;
}


inline
const libMesh::RealGradient&
DSSCModel::get_electric_field(void) const
{
  return _pd.electric_field;
}


inline
const libMesh::RealGradient&
DSSCModel::get_grad_fermi_n(void) const
{
  return _pd.grad_fermi_n;
}


inline
const libMesh::RealGradient&
DSSCModel::get_grad_fermi_I(void) const
{
  return _pd.grad_fermi_I;
}


inline
const libMesh::RealGradient&
DSSCModel::get_grad_fermi_I3(void) const
{
  return _pd.grad_fermi_I3;
}


inline
const libMesh::RealGradient&
DSSCModel::get_grad_fermi_C(void) const
{
  return _pd.grad_fermi_C;
}


inline
void
DSSCModel::set_lattice_temperature(double T)
{
  _pd.kT = T * Constants::k_B;
}


inline
double
DSSCModel::get_lattice_temperature(void) const
{
  return _pd.kT;
}


inline
const Elem*
DSSCModel::get_element(void) const
{
  return _elem;
}


inline
const Point&
DSSCModel::get_coordinates(void) const
{
  return _pd.coordinates;
}


inline
double
DSSCModel::get_charge_density(void) const
{
 return _pd.density_C - get_density_n_exp_DOS() + _trap_DOS * pow(_eq_conc.n/_CB_DOS,_exp_trap) +
    _pd.ionized_dye - _pd.density_I - _pd.density_I3 + (_pd.ionized_electron_traps - _pd.ionized_equilibrium_electron_traps);
      
// return _pd.density_C - _pd.density_n + _eq_conc.n +
//    _pd.ionized_dye - _pd.density_I - _pd.density_I3 + _pd.ionized_electron_traps;
}


inline
double
DSSCModel::get_ionized_electron_traps(void) const
{ 
  return _pd.ionized_electron_traps; 
}


inline
double
DSSCModel::get_ionized_electron_traps_derivative(void) const
{ 
  return _pd.ionized_electron_traps_derivative; 
}


inline
PhysicalModel*
DSSCModel::create_new(void) const
{
  return new DSSCModel(get_options());
}


inline
DSSCModel::PointData&
DSSCModel::get_pd(void)
{
  return _pd;
}


inline
TemperatureInterface&
DSSCModel::get_temperature_interface(void)
{
  return _lattice_temp;
}


inline
double
DSSCModel::get_generation_rate(void) const
{
  return _pd.generation_rate;
}


inline
double
DSSCModel::get_recombination_rate(void) const
{
  return _pd.recombination_rate;
}


inline
double
DSSCModel::get_net_recombination_rate(void) const
{
  return (_pd.recombination_rate - _pd.generation_rate);
}


inline
const std::vector<double>&
DSSCModel::get_net_recombination_rate_derivatives(void) const
{
  return _pd.recombination_rate_derivatives;
}


inline
double
DSSCModel::get_ionized_dye_density(void) const
{
  return _pd.ionized_dye;
}


inline
double
DSSCModel::get_density_C(void) const
{
  return _pd.density_C;
}


inline
double
DSSCModel::get_density_derivative_C(void) const
{
  return -get_density_C() / _pd.kT;
}


inline
double
DSSCModel::get_density_n_exp_DOS(void) const
{
  return _trap_DOS * pow(_pd.density_n/_CB_DOS, _exp_trap);
}


inline
double
DSSCModel::get_density_n(void) const
{
  return _pd.density_n;
}


inline
double
DSSCModel::get_density_derivative_n_exp_DOS(void) const
{
  return _exp_trap * pow(get_density_n(), _exp_trap) * pow( 1/_CB_DOS, _exp_trap) * _trap_DOS / _pd.kT;  
    
//  pow(_trap_DOS, 1 - _exp_trap) * _exp_trap * pow(get_density_n(), _exp_trap) / _pd.kT;
}


inline
double
DSSCModel::get_dens_elec_mobility(void) const
{
  return  ( _mobility.n / _exp_trap ) * pow(_CB_DOS / _trap_DOS, _exp_trap) * pow(get_density_n() / _trap_DOS, 1 - _exp_trap); 
}


inline
double
DSSCModel::get_density_derivative_n(void) const
{
  return get_density_n() / _pd.kT;
}


inline
double
DSSCModel::get_density_I(void) const
{
  return _pd.density_I;
}


inline
double
DSSCModel::get_density_derivative_I(void) const
{
  return get_density_I() / _pd.kT;
}


inline
double
DSSCModel::get_density_I3(void) const
{
  return _pd.density_I3;
}


inline
double
DSSCModel::get_density_derivative_I3(void) const
{
  return get_density_I3() / _pd.kT;
}


inline
const DSSCModel::EquilibriumConcentrations&
DSSCModel::get_equilibrium_concentrations(void) const
{
  return _eq_conc;
}

#endif /* _DSSCMODEL_H_ */
