// $Id$

#ifndef _CARRIERPROPERTIES_H_
#define _CARRIERPROPERTIES_H_

#include "DriftDiffusionModelInterface.h"
#include "TemperatureInterface.h"
#include "DensityOfStates.h"
#include "Constants.h"

#include "libmesh/vector_value.h"

class ModelOptions;
class MobilityModelInterface;


//! Base class for band parameter models
class CarrierProperties : public DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~CarrierProperties(void);

    //! Creator method
    static CarrierProperties* create(const ModelOptions& options);

    //! Get the carrier ID
    ID get_carrier_id(void) const;

    //! Set the temperature interface
    void set_temperature_interface(const TemperatureInterface& temp);

    //! Calculate for a given temperature
    void calculate(double temperature);

    //! Set the temperature
    void set_temperature(double temperature) { _temperature = temperature; };

    //! Get the band edge
    double get_band_edge(void) const;

    //! Get all band energies
    void get_bands(std::vector<double>& bands) const
      { bands = _dos_model->get_reference_energy(); }

    //! Get the effective mass
    double get_effective_mass(void) const;

    //! Get the effective DOS
    double get_effective_DOS(void) const;

    //! Get the maximum density from the DOS model
    double get_maximum_density(void) const;


    //! Get thermal velocity in cm/s
    /*!
     * \param temp the temperature in eV
     */
    double get_thermal_velocity(double temp) const;


    //! Get the particle density and its derivative
    /*!
     * \param Ef the quasi Fermi level
     * \param Epot the electrostatic potential
     *
     * Derivative is given w.r.t to the quasi Fermi level
     */
    std::pair<double, double> get_density_and_derivative(double Ef, double Epot);


    //! Get the conductivity
    double get_conductivity(void) const;

    //! Get the electric conductivity
    /*!
     * the signed value (negative for negatively charged carriers)
     * is returned, without multiplication with the elementary charge
     */
    double get_electric_conductivity(void) const;

    //! Get the mobility
    double get_mobility(void) const;

    /*!
     * \brief Get the conductivity and its derivatives
     */
    double get_conductivity_and_derivatives(double& _deriv_qFermi,
        libMesh::RealGradient& deriv_grad_qFermi,
        libMesh::RealGradient& deriv_grad_potential) const;


    //! Get the thermoelectric power
    double get_thermoelectric_power(void) const;


    //! Set the band edge
    void set_band_edge(double band_edge);

    //! Set the effective mass
    void set_effective_mass(double effective_mass) { _effective_mass = effective_mass; }

    //! Get the carrier type
    const char get_carrier_type(void) const
      { return _particle; };

    //! Get logical name of the particle
    std::string get_particle_name(void) const;

    //! Get the particle charge in units of \c e
    double get_charge(void) const;

    //! Get the particle spin  in units of \c h/2pi
    double get_spin(void) const;

    //! Tell if the carrier is a dopant
    const bool is_dopant(void) const;

    //! Get the sign of the charge
    /*!
     * 0 charge is returned as negative sign, because formulation
     * of carrier flux is such that it is consistent with the
     * formulas for electrons.
     */
    int get_charge_sign(void) const;

    //! Do we have quantum density?
    /*!
     * \note this is only for the time being
     */
    bool has_quantum(void) const;

    //! Should we use quantum?
    void use_quantum(bool use_quantum = true);


  protected:

    //! Constructor
    CarrierProperties(const ModelOptions& options);

    //! Get the temperature in eV
    double get_temperature(void) const { return _temperature; }

    //! Get the lattice temperature in eV
    double get_lattice_temperature(void) const;

    //! Prepare submodels
    virtual void prepare_submodels(void) override;

    //! Initialize
    virtual void do_init(void) override;

    virtual void do_reinit(void) override;

    //! Calculate band properties
    virtual void do_calculate(void) {};

    //! Print some info
    virtual void do_print_info(void) override;


    // ! Get band edge
    //double band_edge(void) const { return _dos_model->reference_energy(); }


    //! Get all band edges
    std::vector<double>& band_edges(void) { return _dos_model->reference_energy(); }


    //! Get effective mass
    double& effective_mass(void) { return _effective_mass; }


    // ! Get degeneracy
    //unsigned int& degeneracy(void) { return _degeneracy; }

    //! Get the DOS factor
    double get_dos_factor(void) const { return _dos_factor; }


  private:

    CarrierProperties(const CarrierProperties& cp);

    //! The different conductivity models available
    enum ConductivityModel
    {
      CONST = 0,        //< Constant conductivity, e.g. for metals
      STD = 1,          //< The standard model, \f[\sigma = \mu n \f]
      GEN = 2           //< The generalized model, \f[\sigma = D \partial n / \partial E_F \f]
    };

    //! The particle this band is describing
    /*! \obsolete
     *
     */
    char _particle;

    //! The logical name of the particle
    std::string _particle_name;

    //! The carrier's unique ID in the module
    ID _carrier_id;


    //! The band edge
    double _band_edge;


    //! The effective mass for the DOS
    /*!
     * It includes any degeneration, i.e. also spin
     */
    double _effective_mass;


    //! The particle charge in units of \c e
    double _charge;

    //! The particle spin in units of \c h/2pi
    double _spin;

    //! Tell if the carrier is a dopant
    bool _is_dopant;

    // ! Charge carrier forming the exciton
    //std::vector<std::string> _exciton_carriers;

    //! The temperature in eV
    double _temperature;

    //! Background conductivity
    /*!
     * Sometimes, the conductivity becomes so small that
     * no solution can be found. In these cases, a small background
     * conductivity can recover convergence.
     */
    double _background_conductivity;

    //! Which conductivity model we use
    ConductivityModel _conductivity_model;

    //! The calculated conductivity
    double _conductivity;

    //! The conductivity derivative w.r.t the quasi Fermi level
    double _conductivity_derivative_qFermi;

    //! The conductivity derivative w.r.t the quasi Fermi level gradient
    libMesh::RealGradient _conductivity_derivative_gradqFermi;

    //! The conductivity derivative w.r.t the electrostatic potential gradient
    libMesh::RealGradient _conductivity_derivative_gradPotential;

    //! The lattice temperature interface
    TemperatureInterface _lattice_temp;

    //! The parabolic band DOS factor
    const double _dos_factor;

    //! The DOS model
    DensityOfStates* _dos_model;

    //! The mobility model
    MobilityModelInterface* _mobility_model;

};

inline
void
CarrierProperties::set_band_edge(double)
{}

inline
CarrierProperties*
CarrierProperties::create(const ModelOptions& options)
{
  return new CarrierProperties(options);
}


inline
ID
CarrierProperties::get_carrier_id(void) const
{
  return(_carrier_id);
}


inline
void
CarrierProperties::set_temperature_interface(const TemperatureInterface& temp)
{
  _lattice_temp = temp;
}



inline
double
CarrierProperties::get_effective_DOS(void) const
{
  return _dos_model->get_effective_dos();
}


inline
double
CarrierProperties::get_conductivity(void) const
{
  return(_conductivity);
}


inline
double
CarrierProperties::get_electric_conductivity(void) const
{
  return(_conductivity*_charge);
}


inline
double
CarrierProperties::get_conductivity_and_derivatives(double& deriv_qFermi,
    libMesh::RealGradient& deriv_grad_qFermi,
    libMesh::RealGradient& deriv_grad_potential) const
{
  deriv_qFermi = _conductivity_derivative_qFermi;
  deriv_grad_qFermi = _conductivity_derivative_gradqFermi;
  deriv_grad_potential = _conductivity_derivative_gradPotential;
  return(_conductivity);
}


inline
double
CarrierProperties::get_thermal_velocity(double temp) const
{
  const double fac = 3 * Constants::e / Constants::me;
  double vth = fac * temp / get_effective_mass();
  //std::cout<<"effective_mass = "<<get_effective_mass()<<std::endl;
  //std::cout<<"effective_mass = "<<_effective_mass<<std::endl;
  return (100.0 * std::sqrt(vth));
}


inline
double
CarrierProperties::get_charge(void) const
{
  return(_charge);
}

inline
int
CarrierProperties::get_charge_sign(void) const
{
  return((_charge > 0) ? 1 : -1);
}


inline
double
CarrierProperties::get_spin(void) const
{
  return(_spin);
}

inline
std::string
CarrierProperties::get_particle_name(void) const
{
  return(_particle_name);
}


/*
inline
std::vector<std::string>
CarrierProperties::get_exciton_carriers(void) const
{
  return(_exciton_carriers);
}
*/

inline
const bool
CarrierProperties::is_dopant(void) const
{
  return(_is_dopant);
}


#endif // _CarrierProperties_H_
