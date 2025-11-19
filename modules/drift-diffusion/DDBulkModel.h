// $Id$

#ifndef _DDBULKMODEL_H_
#define _DDBULKMODEL_H_



#include "DriftDiffusionProperties.h"

//#include "tibercad/physics/TemperatureInterface.h"
#include "tibercad/physics/StrainInterface.h"
#include "tibercad/base/SimulationOptions.h"
#include "DriftDiffusionDefs.h"
#include "tibercad/base/TiberCad.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/base/TypeDefs.h"


#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/point.h"


#include <vector>
#include <set>
#include <map>


// forward declarations
class Elem;
class Dopant;
class Trap;
class SimulationInterface;
class DriftDiffusion;
class MobilityModelInterface;
class ThermoelectricPower;
class PolarizationModel;


//! The base class for all drift-diffusion related semiconductor models
/*!
 * \note { it is not yet possible to have twice the same recombination model.
 * Trying to add to identical models will result in a memory leak. This will
 * be corrected in future. }
 */
class DDBulkModel : public DriftDiffusionProperties
{

  public:

    //! A default (empty) destructor.
    virtual ~DDBulkModel(void);


    //! Create a named drift-diffusion model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static DDBulkModel* create(const std::string& name,
        const Material* mat,
        const ModelOptions& options = ModelOptions());



    //! (Re-)Initialize for the given element
    /*!
     * \c reinit() calls \c prepare_element_data() which needs to be
     * implemented in derived classes
     */
    virtual void do_reinit(const Elem* elem);



    //! Set the polarization vector
    void set_polarization(const libMesh::RealVectorValue& polarization);


    // ! Get the element we are currently working on
    //const Elem* get_element(void) const;


    //! Set the lattice temperature (in K)
    //void set_lattice_temperature(double T);

    //! Get the lattice temperature (in units of eV)
    //double get_lattice_temperature(void) const;



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


    // ! Get the total charge density
    /* !
     *
     * \f$ \rho = p - n + N_D^+ - N_A^- \f$
     *
     * \return the total charge density \f$\rho\f$
     *
     * \note
     * The charge density is returned in units of the elementary charge
     * \f$e\f$, \em not in Coulomb!
     */
    //double get_charge_density(void) const;


    // ! Get the derivatives of the total charge density
    /* !
     * The derivatives are given w.r.t. the two quasi Fermi levels
     */
    //void get_charge_density_derivatives(double derivatives[2]) const;



    //! Get the total electric polarization
    /*!
     * The total electric polarization \b P is the sum of the
     * pyroelectric and piezoelectric polarization
     */
    const libMesh::RealVectorValue& get_total_polarization(void) const
      { return _polarization; };


    //! Get the electron mobility model
    MobilityModelInterface* get_electron_mobility_model(void)
      { return _electron_mobility; };


    //! Get the hole mobility model
    MobilityModelInterface* get_hole_mobility_model(void)
      { return _hole_mobility; };


    //! Returns the thermoelectric power for electrons
    double get_electron_thermoelectric_power() const;

    //! Returns the thermoelectric power for holes
    double get_hole_thermoelectric_power() const;

    //! Computes the electron and hole thermoelectric powers
    void compute_thermoelectric_powers(void);

    //! Computes the electron and hole thermoelectric power derivatives
    void compute_thermoelectric_power_gradient(void);

    //!provides holes thermoelectric power [V/K]
    libMesh::RealGradient get_electron_thermoelectric_power_gradient(void) const;

    //!provides holes thermoelectric power [V/K]
    libMesh::RealGradient get_hole_thermoelectric_power_gradient(void) const;

    //! Get the all nodal temperatures for a given element
    std::vector<double>& get_temperature_at_nodes(void);




  protected:



    //! The empty constructor.
    DDBulkModel(const ModelOptions& options);


    /*! \copydoc PhysicalModel::read_database() */
    virtual void read_database(void);


    //! Initialize this model
    /*!
     * This reads the database and calls init for all submodels
     * A derived class which reimplements this method has to call
     * explicitly the one of this class!
     */
    virtual void do_init(void);


    //! Create some of the submodels
    virtual void prepare_submodels(void);


    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     * This method can be used overiden by derived classes.
     */
    virtual void prepare_element_data(void) {};


    //! \copydoc PhysicalModel::do_print_info(void)
    virtual void do_print_info(void);


    //! Get the temperature interface
    //TemperatureInterface& get_temperature_interface(void);

    //! Get the strain interface
    StrainInterface& get_strain_interface(void);


    //! Tells if we should assume inhomogeneous band parameters
    bool is_inhomogeneous(void) const;


    //! Tells if we are doing equilibrium calculation
    bool has_solution(void) const;


    //! Tells if we should use a predictor for quantum densities
    bool use_predictor(void) const;


    //! Set the intrinsic Fermi level
    /*!
     * Calculates the associated equilibrium densities
     */
    void set_equilibrium_properties(double Ef);


  private:

    //! The interface to the lattice temperature simulation
    //TemperatureInterface _lattice_temp;


    //! The interface to a strain simulation
    StrainInterface _strain_if;


    //! \c true if we should assume inhomogeneous band parameters
    bool _is_inhomogeneous;

    //! \c true if we should use a predictor for quantum densities
    bool _use_predictor;



    //! Electron thermoelectric power gradient
    libMesh::RealGradient _eTEpowerGrad;

    //! Hole thermoelectric power gradient
    libMesh::RealGradient _hTEpowerGrad;

    //! Electron thermoelectric power
    double _eTEpower;

    //! Hole thermoelectric power
    double _hTEpower;


    //! The total electric polarization
    libMesh::RealVectorValue _polarization;


    //! A background conductivity for electrons
    /*!
     * Mainly used for stability reasons in pathologic cases
     * Units are S/cm
     */
    double _background_conductivity_e;

    //! A background conductivity for holes
    /*!
     * Mainly used for stability reasons in pathologic cases
     * Units are S/cm
     */
    double _background_conductivity_h;



    //! The copy constructor is disabled
    DDBulkModel(const DDBulkModel& rhs);


    //! The assignment operator is disabled
    DDBulkModel& operator=(const DDBulkModel& rhs);


    //! Parse the model options
    void parse_options(void);



    //! Create a mobility model
    /*!
     * Creates a mobility model from the given model name
     * and options.
     */
    MobilityModelInterface* create_mobility_model(
        const ModelOptions& options = ModelOptions());



    //! The element we are currently working on
    //const Elem* _elem;


    //! The electron traps
    std::set<Trap*> _etraps;

    //! The hole traps
    std::set<Trap*> _htraps;


    //! The electron mobility
    MobilityModelInterface* _electron_mobility;

    //! The hole mobility
    MobilityModelInterface* _hole_mobility;


    //! The thermoelectric power
    ThermoelectricPower* _thermoelectric_power;


    //! The polarization models
    std::vector<PolarizationModel*> _pm;


    //! The lattice temperature in eV (\f$= k_B T_{lat} / e\f$)
    //double _lattice_vt;


    //! The nodal lattice temperature
    std::vector<double> _nodal_lattice_vt;


    //! The relaxation factor for the polarization
    double _relax_polariz;


};


//
// inline members
//










//inline
//void
//DDBulkModel::set_lattice_temperature(double T)
//{
//  _lattice_vt = T * Constants::k_B;
//}



//inline
//double
//DDBulkModel::get_lattice_temperature(void) const
//{
//  return _lattice_vt;
//}


/*
inline
const Elem*
DDBulkModel::get_element(void) const
{
  return _elem;
}
*/


/*
inline
double
DDBulkModel::get_charge_density(void) const
{
  long double dens = static_cast<long double>(get_pd().hole_density) -
      static_cast<long double>(get_pd().electron_density) +
      static_cast<long double>(get_pd().ionized_donor_density) -
      static_cast<long double>(get_pd().ionized_acceptor_density) +
      static_cast<long double>(get_pd().ionized_electron_traps) +
      static_cast<long double>(get_pd().ionized_hole_traps);

  return static_cast<double>(dens);
  //return get_hole_density() - get_electron_density() +
  //  get_ionized_donor_density() - get_ionized_acceptor_density()
  //  + get_ionized_electron_traps() + get_ionized_hole_traps();
}


inline
void
DDBulkModel::get_charge_density_derivatives(
    double derivatives[2]) const
{
  long double der0 = static_cast<long double>(get_electron_density_derivative())
      - static_cast<long double>(get_ionized_donor_density_derivative())
      - static_cast<long double>(get_pd().ionized_traps_derivative[0]);

  long double der1 = -static_cast<long double>(get_hole_density_derivative())
      + static_cast<long double>(get_ionized_acceptor_density_derivative())
      - static_cast<long double>(get_pd().ionized_traps_derivative[1]);

  derivatives[0] = static_cast<double>(der0);
  derivatives[1] = static_cast<double>(der1);

  //derivatives[0] = get_electron_density_derivative()
  //  - get_ionized_donor_density_derivative()
  //  - _pd->ionized_electron_traps_derivative;
  //derivatives[1] = -get_hole_density_derivative()
  //  + get_ionized_acceptor_density_derivative()
  //  - _pd->ionized_hole_traps_derivative;
}
*/




inline
void
DDBulkModel::set_polarization(const libMesh::RealVectorValue& polarization)
{
  _polarization = _relax_polariz * polarization;
}




inline
double
DDBulkModel::get_electron_thermoelectric_power(void) const
{
  return _eTEpower;
}


inline
double
DDBulkModel::get_hole_thermoelectric_power(void) const
{
  return _hTEpower;
}

inline
libMesh::RealGradient
DDBulkModel::get_electron_thermoelectric_power_gradient(void) const
{

  return  _eTEpowerGrad;

}

inline
libMesh::RealGradient
DDBulkModel::get_hole_thermoelectric_power_gradient(void) const
{

  return  _hTEpowerGrad;

}



/*
inline
TemperatureInterface&
DDBulkModel::get_temperature_interface(void)
{
  return _lattice_temp;
}
*/

inline
StrainInterface&
DDBulkModel::get_strain_interface(void)
{
  return _strain_if;
}


inline
bool
DDBulkModel::is_inhomogeneous(void) const
{
  return _is_inhomogeneous;
}


inline
bool
DDBulkModel::use_predictor(void) const
{
  return _use_predictor;
}








#endif /* _DDBULKMODEL_H_ */
