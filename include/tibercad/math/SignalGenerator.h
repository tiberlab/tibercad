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
 * \file SignalGenerator.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#include "tibercad/module/TiberModelObject.h"

/*!
 * \brief a class to define signals
 *
 * Objects derived from this class can be used to define input signals,
 * in particular time dependent data for boundary conditions. Practically,
 * the implementation calculates a dependent quantity from the value of an
 * input variable (as default the global time, \c $time), and provides it
 * as a possibly user defined variable.
 * Note that the base class only defines the input variable, which must be
 * a real value. The output variable has to be defined in the derived
 * classes since it's type might change.
 */
class SignalGenerator : public TiberModelObject
{

  public:

    //! Create a SignalGenerator object
    static SignalGenerator* create(const ModelOptions& options);

    //! Destructor
    virtual ~SignalGenerator(void);

    //! Initialize
    void init(void);

    //! Calculate new value of dependent variables
    void update_dependent_variables(void);

  protected:

    //! Constructor
    SignalGenerator(const ModelOptions& options);

    //! Calculate new value of dependent variables
    virtual void do_update_dependent_variables(void) = 0;

    //! Implementation specific initialization
    virtual void do_init(void) = 0;

    //! Get the value of the input variable
    double get_input(void) const;

  private:

    //! The input variable
    double _input;
};


inline
double
SignalGenerator::get_input(void) const
{
  return(_input);
}
