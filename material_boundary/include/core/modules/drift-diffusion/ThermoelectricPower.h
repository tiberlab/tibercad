// $Id$



#ifndef _THERMOELECTRICPOWER_H_
#define _THERMOELECTRICPOWER_H_


//#ifndef TIBER_MODULE_NAME
//# define TIBER_MODULE_NAME dd_thelpow
//#endif



#include "DriftDiffusionModelInterface.h"
#include "vector_value.h"

//! This class computes the thermoelectric power
/*!
 *
 * The thermoelectric power must be in [V/K]
 *
 * \f$ P_n = -\frac{k_b}{q}\left(\frac{5}{2} + \alpha_n +
 *    E_c + q \varphi_n \right) \f$
 *
 * \f$ P_p = +\frac{k_b}{q}\left(\frac{5}{2} + \alpha_p -
 *    q \varphi_n - E_v \right) \f$
 *
 * where
 *  \f$  \alpha = \frac{T}{\mu_n} \frac {\partial \mu_n}{\partial T}\f$
 *  \f$  \alpha = \frac{T}{\mu_p} \frac {\partial \mu_p}{\partial T}\f$
 *
 */
class ThermoelectricPower : public DriftDiffusionModelInterface
{

  public:

    //! Constructor 
    ThermoelectricPower(void);

    //! Destructor
    ~ThermoelectricPower(void) { };

    //! set the electro-chemical potential for electrons and holes
    void set_potentials(double eQfermi, double hQfermi, double ElPot);

    ////! set the electrostatic potential 
    // void set_electric_potential(double phi);

    //! set the band edge of conduction band and valence band
    void set_band_edges(double Ec, double Ev);

    //! set the local temperature
    void set_temperature(double Tloc);

    //!provides electrons thermoelectric power [V/K]
    double get_electrons_thermoelectric_power(void) const;

    //!provides holes thermoelectric power [V/K]
    double get_holes_thermoelectric_power(void) const;

    //!provides holes thermoelectric power [V/K]
    RealGradient get_electron_thermoelectric_power_gradient(void) const;

    //!provides holes thermoelectric power [V/K]
    RealGradient get_hole_thermoelectric_power_gradient(void) const;


    //! Calculate the value of the thermoelectric power
    void calculate(void);

    //! Calculate the thermoelectric power derivatives
    void calculate_derivatives(void);


    static ThermoelectricPower* create_model(const std::string& model,
        const ModelOptions& options = ModelOptions());

    //! Set the electron and hole charge density
    //  void set_charge_densities(double n, double p);

    //! Set the electron and hole charge density
    void set_potential_gradients(RealGradient eFermiGrad, 
        RealGradient hFermiGrad, 
        RealGradient ElectricField);




  protected:

    virtual void do_init (void);

    virtual void read_database(void);

    virtual void do_init_alloy(const PhysicalModelInterface *comp_A,
        const PhysicalModelInterface *comp_B, double xa); 

    virtual PhysicalModelInterface* create_new(void) const;


  private:

    enum TEPModel
    {
      UNKNOWN = 0, /*!< unknown model */
      CONSTANT,    /*!< constant model */
      DIFFUSIVITY, /*!< diffusivity model */
    };

    //! Lattice temperature in eV
    double _Tloc; 

    double _eQfermi;

    double _hQfermi;

    double _Ec;

    double _Ev;

    double _ElPot;

    RealGradient _eFermiGrad;

    RealGradient _hFermiGrad;

    RealGradient _ElectricField;

    RealGradient _eTEpowerGrad;

    RealGradient _hTEpowerGrad;


    //! Model for thermoelectric power 
    TEPModel _TEmodel;

    //! Electron thermoelectric power in simulation system. Units V/K
    double _eTEpower;

    //! Hole thermoelectric power in simulation system. Units V/K
    double _hTEpower;

};



inline
double
ThermoelectricPower::get_electrons_thermoelectric_power(void) const
{
  return _eTEpower;
}


inline
double
ThermoelectricPower::get_holes_thermoelectric_power(void) const
{
  return _hTEpower;
}

inline
RealGradient
ThermoelectricPower::get_electron_thermoelectric_power_gradient(void) const
{
  return _eTEpowerGrad;
}


inline
RealGradient
ThermoelectricPower::get_hole_thermoelectric_power_gradient(void) const
{
  return _hTEpowerGrad;
}


inline
void
ThermoelectricPower::set_potentials(double eQfermi, double hQfermi, double ElPot)
{
  _eQfermi = eQfermi;
  _hQfermi = hQfermi;
  _ElPot   = ElPot;
}


inline
void 
ThermoelectricPower::set_potential_gradients(RealGradient eFermiGrad, 
    RealGradient hFermiGrad, 
    RealGradient ElectricField )
{

  _eFermiGrad = eFermiGrad;
  _hFermiGrad = hFermiGrad;
  _ElectricField = ElectricField;
}




inline
void
ThermoelectricPower::set_band_edges(double Ec, double Ev)
{
  _Ec = Ec;

  _Ev = Ev;

}


inline
void
ThermoelectricPower::set_temperature(double Tloc)
{

  _Tloc = Tloc;

}



inline
PhysicalModelInterface*
ThermoelectricPower::create_new(void) const
{
  return (new ThermoelectricPower()); 
}


#endif
