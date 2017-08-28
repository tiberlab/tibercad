// $Id$

#ifndef _CARRIERPROPERTIES_H_
#define _CARRIERPROPERTIES_H_

#include "DriftDiffusionModelInterface.h"
#include "TemperatureInterface.h"
#include "DensityOfStates.h"
#include "Constants.h"

class ModelOptions;


//! Base class for band parameter models
class CarrierProperties : public DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~CarrierProperties(void);

    //! Creator method
    static CarrierProperties* create(const ModelOptions& options);

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
    std::pair<double, double> get_density_and_derivative(double Ef, double Epot) const;

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

    // ! Return the charge carriers names associated to an exciton
    //std::vector<std::string> get_exciton_carriers(void) const;


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

    //! The lattice temperature interface
    TemperatureInterface _lattice_temp;

    //! The parabolic band DOS factor
    const double _dos_factor;

    //! The DOS model
    DensityOfStates* _dos_model;

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
