// $Id$

#ifndef _FOWLERNORDHEIM_H_
#define _FOWLERNORDHEIM_H_


//! Implements the Fowler-Nordheim field emission model
/*!
 * Formulas from [Zuber, Jensen, Sullivan; JAP, vol. 91, no. 11, 2002]
 */
class FowlerNordheim
{

  public:

    //! The Constructor
    FowlerNordheim(double workfunction = 1.0);

    //! Destructor
    ~FowlerNordheim(void) {};


    //! Set the metal workfunction
    void set_workfunction(double workfunction);


    //! Get the emission current density for a given field strength
    double get_emission_current(double F);


    //! Get the emission particle velocity
    void set_velocity(double velocity);

    //! Get the emission particle velocity
    double get_velocity(void) const;


  private:

    //! The metal work function in eV
    double _workfunction;


    //! Emission velocity in cm/s
    double _velocity;

};


//
// inline members
//

inline
FowlerNordheim::FowlerNordheim(double workfunction)
  : _workfunction(workfunction)
{
}


inline
void
FowlerNordheim::set_workfunction(double workfunction)
{
  _workfunction = workfunction;
}


inline
double
FowlerNordheim::get_velocity(void) const
{
  return _velocity;
}

inline
void
FowlerNordheim::set_velocity(double velocity)
{
  _velocity = velocity;
}

#endif // _FOWLERNORDHEIM_H_
