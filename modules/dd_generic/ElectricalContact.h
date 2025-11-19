// $Id: ElectricalContact.h 4135 2015-09-25 10:19:38Z maufder $

#ifndef _ELECTRICALCONTACT_H_
#define _ELECTRICALCONTACT_H_


#include "DDInterfaceModel.h"
#include "tibercad/base/ModelOptions.h"
#include "tibercad/base/ModelErrorException.h"

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


    /*!
     * \brief Set the contact fermilevel
     */
    void set_contact_fermilevel(double Ef);



    //! Set the recombination velocity for all carriers
    /*!
     * A value of -1 will not change the current value.
     */
    void set_recombination_velocities(double vrec);


    //! Set the recombination velocity for one carrier
    /*!
     * A value of -1 will not change the current value.
     *
     * \param carrier the global carrier ID
     */
    void set_recombination_velocities(ID carrier, double vrec);


    //! Get the recombination velocities
    /*!
     * The ordering in the vector is the same as the global ordering of
     * carrier quasi Fermi variables
     */
    void get_recombination_velocities(std::vector<double>& vrec) const;



  private:

    //! The boundary value (eg. applied voltage)
    double _voltage;

    //! The contact surface resistance
    double _surfres;

    //! The contact fermilevel
    double _contact_fermilevel;

    //! The carriers surface recombination velocities
    std::vector<double> _vrec;

    //! vrec has been overridden
    std::vector<bool> _fixed_vrec;

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
ElectricalContact::set_recombination_velocities(double vrec)
{
  for (unsigned int i = 0; i < n_known_carriers(); i++)
  {
    if ((vrec > 0) && !_fixed_vrec[i])
    {
      _vrec[i] = vrec;
      set_type(i, NEUMANN);
    }
  }
}

inline
void
ElectricalContact::set_recombination_velocities(ID carrier, double vrec)
{
  if (carrier >= _vrec.size())
    throw ModelErrorException("ElectricalContact: invalid carrier ID provided");

  if ((vrec > 0) && !_fixed_vrec[carrier])
  {
    _vrec[carrier] = vrec;
    set_type(carrier, NEUMANN);
  }
}


inline
void
ElectricalContact::get_recombination_velocities(std::vector<double>& vrec) const
{
  vrec = _vrec;
}



#endif // _ELECTRICALCONTACT_H_
