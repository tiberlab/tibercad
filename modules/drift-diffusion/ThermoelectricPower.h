/*  
 * This file is part of the tiberCAD module driftdiffusion.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file ThermoelectricPower.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */




#ifndef TC_THERMOELECTRICPOWER_H
#define TC_THERMOELECTRICPOWER_H


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
class TBDLLOCAL ThermoelectricPower : public DriftDiffusionModelInterface
{

  public:

    //! Constructor 
    ThermoelectricPower(const ModelOptions& options);

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
    libMesh::RealGradient get_electron_thermoelectric_power_gradient(void) const;

    //!provides holes thermoelectric power [V/K]
    libMesh::RealGradient get_hole_thermoelectric_power_gradient(void) const;


    //! Calculate the value of the thermoelectric power
    void calculate(void);

    //! Calculate the thermoelectric power derivatives
    void calculate_derivatives(void);


    static ThermoelectricPower* create_model(const std::string& model,
        const Material* mat, const ModelOptions& options = ModelOptions());

    //! Set the electron and hole charge density
    //  void set_charge_densities(double n, double p);

    //! Set the electron and hole charge density
    void set_potential_gradients(libMesh::RealGradient eFermiGrad,
        libMesh::RealGradient hFermiGrad,
        libMesh::RealGradient ElectricField);




  protected:

    virtual void do_init (void);

    virtual void read_database(void);

    virtual PhysicalModel* create_new(void) const;


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

    libMesh::RealGradient _eFermiGrad;

    libMesh::RealGradient _hFermiGrad;

    libMesh::RealGradient _ElectricField;

    libMesh::RealGradient _eTEpowerGrad;

    libMesh::RealGradient _hTEpowerGrad;


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
libMesh::RealGradient
ThermoelectricPower::get_electron_thermoelectric_power_gradient(void) const
{
  return _eTEpowerGrad;
}


inline
libMesh::RealGradient
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
ThermoelectricPower::set_potential_gradients(libMesh::RealGradient eFermiGrad,
    libMesh::RealGradient hFermiGrad,
    libMesh::RealGradient ElectricField )
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
PhysicalModel*
ThermoelectricPower::create_new(void) const
{
  return (new ThermoelectricPower(get_options()));
}


#endif
