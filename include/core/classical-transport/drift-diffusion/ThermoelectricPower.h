// $Id$



#ifndef _THERMOELECTRICPOWER_H_
#define _THERMOELECTRICPOWER_H_



#include "DriftDiffusionModelInterface.h"




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
    void set_fermi_potential(double eQfermi, double hQfermi);

    //! set the electrostatic potential 
    void set_electric_potential(double phi);

    //! set the band edge of conduction band and valence band
    void set_band_edges(double Ec, double Ev);

    //! set the local temperature
    void set_temperature(double Tloc);

    //!provides electrons thermoelectric power [V/K]
    double get_electrons_thermoelectric_power(void) const;

    //!provides holes thermoelectric power [V/K]
    double get_holes_thermoelectric_power(void) const;

    //! Calculate the value of the thermoelectric power
    void calculate(void);

    static ThermoelectricPower* create();



  protected:

    virtual void do_init (void);

    virtual void copy_from(const PhysicalModelInterface *rhs);

    virtual void read_database(void);

    //virtual void read_bowing_parameters(void) {};

    virtual void calculate_VCA (const PhysicalModelInterface *comp_A,
        const PhysicalModelInterface *comp_B, double xa); 


    virtual PhysicalModelInterface* create_new (void) const;


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

    double _phi;

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
void
ThermoelectricPower::set_fermi_potential(double eQfermi, double hQfermi)
{
  _eQfermi = eQfermi;
  _hQfermi = hQfermi;
}


inline
void
ThermoelectricPower::set_electric_potential(double phi)
{
  _phi = phi;
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
ThermoelectricPower*
ThermoelectricPower::create()
{
  return (new ThermoelectricPower());
}



inline
PhysicalModelInterface*
ThermoelectricPower::create_new() const
{
  return (new ThermoelectricPower()); 
}




#endif
