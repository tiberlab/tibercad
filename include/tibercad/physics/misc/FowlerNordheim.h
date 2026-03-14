/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file FowlerNordheim.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_FOWLERNORDHEIM_H
#define TC_FOWLERNORDHEIM_H


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

#endif // TC_FOWLERNORDHEIM_H
