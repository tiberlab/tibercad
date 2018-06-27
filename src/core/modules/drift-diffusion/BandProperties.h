// $Id$

#ifndef _BANDPROPERTIES_H_
#define _BANDPROPERTIES_H_

#include "DriftDiffusionModelInterface.h"
#include "TemperatureInterface.h"
#include "DensityOfStates.h"
#include "Constants.h"

class ModelOptions;
//class DensityOfStates;
class ParticleDensity;
class Elem;
class Point;


//! Base class for band parameter models
class BandProperties : public DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~BandProperties(void);

    //! Creator method
    static BandProperties* create(const ModelOptions& options);

    //! Set the temperature interface
    void set_temperature_interface(const TemperatureInterface& temp);

    //! Calculate for a given temperature
    void calculate(double temperature);

    //! Set the temperature
    void set_temperature(double temperature) { _temperature = temperature; }

    //! Get the band edge
    double get_band_edge(void) const;

    //! Get all band energies
    void get_bands(std::vector<double>& bands) const
      { bands = _dos_model->get_reference_energy(); }

    //! Get the effective mass
    double get_effective_mass(void) const;

    //! Get the effective DOS
    double get_effective_DOS(void) const;

    // ! Get degeneracy
    //unsigned int get_degeneracy(void) const { return _degeneracy; }

    //! Get thermal velocity in cm/s
    /*!
     * \param temp the temperature in eV
     */
    double get_thermal_velocity(double temp) const;


    //! Get the particle density and its derivative
    /*!
     * calculate() has to be called before
     */
    void get_density_and_derivative(void) const;

    //! Get the particle density and its derivative
    void get_density_and_derivative(std::vector<double>& den_and_der, double Ef, double Epot) const;

    //! Get the \f$\gamma\f$ factor
    double get_gamma(void) const;

    //! Set the band edge
    void set_band_edge(double band_edge);

    //! Set the effective mass
    void set_effective_mass(double effective_mass) { _effective_mass = effective_mass; }


    //! Do we have quantum density?
    /*!
     * \note this is only for the time being
     */
    bool has_quantum(void) const;

    //! Should we use quantum?
    void use_quantum(bool use_quantum = true);


  protected:

    //! Constructor
    BandProperties(const ModelOptions& options);

    //! Get the temperature in eV
    double get_temperature(void) const { return _temperature; }

    //! Get the lattice temperature in eV
    double get_lattice_temperature(void) const;

    //! Prepare submodels
    virtual void prepare_submodels(void);

    //! Initialize
    virtual void do_init(void) {};

    //! Calculate band properties
    virtual void do_calculate(void) {};

    //! Print some info
    virtual void do_print_info(void);


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
    char _particle;


    //! The band edge
    double _band_edge;


    //! The effective mass for the DOS
    /*!
     * It includes any degeneration, i.e. also spin
     */
    double _effective_mass;

    // ! The total degeneracy, including spin
    //unsigned int _degeneracy;

    //! The temperature in eV
    double _temperature;

    //! The lattice temperature interface
    TemperatureInterface _lattice_temp;

    //! The parabolic band DOS factor
    const double _dos_factor;

    //! The DOS model
    DensityOfStates* _dos_model;

    //! The particle density
    ParticleDensity* _density;

};


inline
BandProperties*
BandProperties::create(const ModelOptions& options)
{
  return new BandProperties(options);
}


inline
void
BandProperties::set_temperature_interface(const TemperatureInterface& temp)
{
  _lattice_temp = temp;
}



inline
double
BandProperties::get_effective_DOS(void) const
{
  return _dos_model->get_effective_dos();
}



inline
double
BandProperties::get_thermal_velocity(double temp) const
{
  const double fac = 3 * Constants::e / Constants::me;
  double vth = fac * temp / get_effective_mass();
  //std::cout<<"effective_mass = "<<get_effective_mass()<<std::endl;
  //std::cout<<"effective_mass = "<<_effective_mass<<std::endl;
  return (100.0 * std::sqrt(vth));
}

inline
void
BandProperties::set_band_edge(double )
{
}


#endif // _BANDPROPERTIES_H_
