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
 * \file ElectricalContact.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_ELECTRICALCONTACT_H
#define TC_ELECTRICALCONTACT_H


#include "DDInterfaceModel.h"
#include "tibercad/base/ModelOptions.h"

// C++ includes
#include <string>

class DriftDiffusionProperties;
class FowlerNordheim;

//! Base class for any kind of electrical contact models
/*!
 * Derive from this class to define contact models. Examples can be found
 * in OhmicContact and SchottkyContact.
 *
 * This class is derived also from Variable to be able to make
 * a voltage sweep.
 */
class TC_DLEXPORT ElectricalContact : public DDInterfaceModel
{
  public:

    //! Empty destructor
    virtual ~ElectricalContact(void) {};


    //! Set the simulation voltage for this contact
    void set_simulation_voltage(double voltage);

    //! Get the simulation voltage for this contact
    double get_simulation_voltage(void) const;

    //! Get the inner voltage for this contact
    /*!
     * This takes surface resistance into account
     */
    double get_inner_voltage(void) const;


    //! Set a BC to homogeneous von Neumann type
    /*!
     * This is only useful if one has a Dirichlet type boundary condition
     * and wants to change it into homogeneous von Neumann.
     * It can be used to have different BCs for electron and hole
     * electro-chemical potentials in small devices, or to change to zero
     * electric field BC
     */
    //void set_zero_derivative_bc(DriftDiffusionDefs::DDVariable variable);



  protected:

    //! This class is not intended for direct use
    ElectricalContact(const ModelOptions& options);

    //! Initialize
    virtual void do_init(void);


    //! Setup the coefficients
    /*!
     * This implementation works for ohmic and Schottky contacts,
     * if the barrier has been set before.
     */
    virtual void do_compute(void);


    //! Get the contact voltage drop
    /*!
     * A positive value means: the inner voltage is lower than
     * the outer contact voltage
     */
    //double get_contact_voltage_drop() const;


    /*!
     * \brief Set the contact fermilevel
     */
    void set_contact_fermilevel(double Ef);



    //! Set the recombination velocities
    /*!
     * A value of -1 will not change the current value.
     */
    void set_recombination_velocities(double vn, double vp);


    //! Get the recombination velocities
    /*!
     * A value of -1 will not change the current value.
     */
    void get_recombination_velocities(double& vn, double& vp);



  private:

    //! The boundary value (eg. applied voltage)
    double _voltage;

    //! The contact surface resistance
    double _surfres;

    //! The contact fermilevel
    double _contact_fermilevel;

    //! The electron surface recombination velocity
    double _vrec_n;

    //! The hole surface recombination velocity
    double _vrec_p;

    //! vrec_n has been overridden
    bool _fixed_vrec_n;

    //! vrec_p has been overridden
    bool _fixed_vrec_p;

    //! tunneling option
    bool _tunneling;
    
    //! tunneling parameters for electrons and holes
    double _bias;
    double _alpha_n;
    double _mass_n;
    double _phi_n;
    double _alpha_p;
    double _phi_p;
    double _mass_p;
    double _Vbi;
    double _length;
};




//
// inline members
//



inline
void
ElectricalContact::set_simulation_voltage(double voltage)
{
  _voltage = voltage;
}


inline
double
ElectricalContact::get_simulation_voltage(void) const
{
  return _voltage;
}



inline
double
ElectricalContact::get_inner_voltage(void) const
{
  double v = get_simulation_voltage();
  //if (_surfres > 1e-12)
  //  v -= get_contact_voltage_drop();

  return v;
}


inline
void
ElectricalContact::set_contact_fermilevel(double Ef)
{
  _contact_fermilevel = Ef;
}


inline
void
ElectricalContact::set_recombination_velocities(double vn, double vp)
{
  if ((vn > 0) && !_fixed_vrec_n)
  {
    _vrec_n = vn;
    set_type(1, NEUMANN);
  }
  if ((vp > 0) && !_fixed_vrec_p)
  {
    _vrec_p = vp;
    set_type(2, NEUMANN);
  }
}


inline
void
ElectricalContact::get_recombination_velocities(double& vn, double& vp)
{
  vn = _vrec_n;
  vp = _vrec_p;
}



#endif // TC_ELECTRICALCONTACT_H
