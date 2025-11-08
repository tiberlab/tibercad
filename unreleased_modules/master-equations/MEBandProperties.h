// $Id$

#ifndef _MEBANDPROPERTIES_H_
#define _MEBANDPROPERTIES_H_

#include "MasterEquationsModelInterface.h"
#include "TemperatureInterface.h"
#include "DensityOfStates.h"
#include "Constants.h"

class ModelOptions;
//class DensityOfStates;
class ParticleDensity;
class Elem;
class Point;


//! Base class for band parameter models
class MEBandProperties : public MasterEquationsModelInterface // ??
{

  public:

    //! Destructor
    virtual ~MEBandProperties(void);

    //! Creator method
    static MEBandProperties* create(const ModelOptions& options);

    //! Set the temperature interface
    void set_temperature_interface(const TemperatureInterface& temp);

    //! Calculate for a given temperature
    void calculate(double temperature);

    //! Set the temperature
    void set_temperature(double temperature) { _temperature = temperature; }

    //! Get the band edge
    double get_band_edge(void) const;

    //! Get the effective mass
    double get_effective_mass(void) const;

    //! Get the effective DOS
    double get_effective_DOS(void) const;

    //! Set the band edge
    void set_band_edge(double band_edge) {  }

    //! Set the effective mass
    void set_effective_mass(double effective_mass) { _effective_mass = effective_mass; }

    //! Get thermal velocity in cm/s
    /*!
     * \param temp the temperature in eV
     */
    double get_thermal_velocity(double temp) const;


  protected:

    //! Constructor
    MEBandProperties(const ModelOptions& options);

    //! Get the temperature in eV
    double get_temperature(void) const { return _temperature; }

    //! Get the lattice temperature in eV
    double get_lattice_temperature(void) const;

    //! Initialize
    virtual void do_init(void) {};

    //! Calculate band properties
    virtual void do_calculate(void) {};

    //! Print some info
    virtual void do_print_info(void);


    //! Get all band edges
    std::vector<double>& band_edges(void) { return _dos_model->reference_energy(); }


    //! Get effective mass
    double& effective_mass(void) { return _effective_mass; }


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
MEBandProperties*
MEBandProperties::create(const ModelOptions& options)
{
  return new MEBandProperties(options);
}


inline
void
MEBandProperties::set_temperature_interface(const TemperatureInterface& temp)
{
  _lattice_temp = temp;
}



inline
double
MEBandProperties::get_effective_DOS(void) const
{
  return _dos_model->get_effective_dos();
}


inline
double
MEBandProperties::get_thermal_velocity(double temp) const
{
  const double fac = 3 * Constants::e / Constants::me;
  double vth = fac * temp / get_effective_mass();
  //std::cout<<"effective_mass = "<<get_effective_mass()<<std::endl;
  //std::cout<<"effective_mass = "<<_effective_mass<<std::endl;
  return (100.0 * std::sqrt(vth));
}


#endif // _MEBANDPROPERTIES_H_
