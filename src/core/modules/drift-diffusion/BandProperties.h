// $Id$

#ifndef _BANDPROPERTIES_H_
#define _BANDPROPERTIES_H_

#include "DriftDiffusionModelInterface.h"
#include "Constants.h"

class ModelOptions;


//! Base class for band parameter models
class BandProperties : public DriftDiffusionModelInterface
{

  public:

    // temporary
    BandProperties();

    //! Destructor
    virtual ~BandProperties(void) {}

    //! Calculate for a given temperature
    void calculate(double temperature);

    //! Set the temperature
    void set_temperature(double temperature) { _temperature = temperature; }

    //! Get the band edge
    double get_band_edge(void) const { return _band_edge; }

    //! Get all band energies
    void get_bands(std::vector<double>& bands) const {};

    //! Get the effective mass
    double get_effective_mass(void) const { return _effective_mass; }

    //! Get the effective DOS
    double get_effective_DOS(void) const;

    //! Get degeneracy
    unsigned int get_degeneracy(void) const { return _degeneracy; }

    //! Get thermal velocity in cm/s
    /*!
     * \param temp the temperature in eV
     */
    double get_thermal_velocity(double temp) const;



  protected:

    //! Constructor
    BandProperties(const ModelOptions& options);

    //! Get the temperature
    double get_temperature(void) const { return _temperature; }

    //! Initialize
    virtual void do_init(void) {};

    //! Calculate band properties
    virtual void do_calculate(void) {};


    //! Get band edge
    double& band_edge(void) { return _band_edge; }


    //! Get effective mass
    double& effective_mass(void) { return _effective_mass; }


    //! Get degeneracy
    unsigned int& degeneracy(void) { return _degeneracy; }

    //! Get the DOS factor
    double get_dos_factor(void) const { return _dos_factor; }


  //private:
  public:

    //! The band edge
    double _band_edge;

    //! The effective mass for the DOS
    /*!
     * It includes any degeneration, i.e. also spin
     */
    double _effective_mass;

    //! The effective density of states
    //double _effective_DOS;

    //! The degeneracy
    unsigned int _degeneracy;

    //! The temperature in eV
    double _temperature;

    const static double _dos_factor;
};


inline
void
BandProperties::calculate(double temperature)
{
  set_temperature(temperature);
  do_calculate();
}

inline
double
BandProperties::get_effective_DOS(void) const
{
  return _dos_factor * std::pow(_temperature * get_effective_mass(), 1.5);
}


inline
double
BandProperties::get_thermal_velocity(double temp) const
{
  const double fac = 3 * Constants::e / Constants::me;
  double vth = fac * temp / _effective_mass;
  return (100.0 * std::sqrt(vth));
}

#endif // _BANDPROPERTIES_H_
