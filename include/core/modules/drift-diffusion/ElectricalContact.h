// $Id$

#ifndef _ELECTRICALCONTACT_H_
#define _ELECTRICALCONTACT_H_


#include "DDInterfaceModel.h"
#include "ModelOptions.h"

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
class TBDLEXPORT ElectricalContact : public DDInterfaceModel
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
     * This is only useful if on has a Dirichlet type boundary condition
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
     * \brief Set the barrier (difference between metal fermi
     *   level and conduction band edge)
     */
    void set_barrier(double barrier);


    //! Set the recombination velocities
    /*!
     * A value of -1 will not change the current value.
     */
    void set_recombination_velocities(double vn, double vp);



  private:

    //! The boundary value (eg. applied voltage)
    double _voltage;


    //! The contact surface resistance
    double _surfres;


    //! The barrier (difference between metal fermi level and conduction band edge)
    double _barrier;


    //! The electron surface recombination velocity
    double _vrec_n;


    //! The hole surface recombination velocity
    double _vrec_p;


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
ElectricalContact::set_barrier(double barrier)
{
  _barrier = barrier;
}


inline
void
ElectricalContact::set_recombination_velocities(double vn, double vp)
{
  if (vn > 0) _vrec_n = vn;
  if (vp > 0) _vrec_p = vp;
}


#endif // _ELECTRICALCONTACT_H_
