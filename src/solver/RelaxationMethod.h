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
 * \file RelaxationMethod.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef _RELAXATIONMETHOD_H_
#define _RELAXATIONMETHOD_H_

#include "SelfconsistentSolver.h"


//! Interface for self-consistent calculations
class TBDLLOCAL RelaxationMethod : public SelfconsistentSolver
{

  public:

    //! The destructor
    virtual ~RelaxationMethod(void);

    //! Create a Sweep object
    static RelaxationMethod* create(const ModelOptions& options);


  protected:

    //! The empty Constructor
    RelaxationMethod(const ModelOptions& options);

    
    /*! \copydoc SimulationInterface::do_solve() */
    virtual void do_solve(void);

    
    /*! \copydoc SimulationInterface::parse_options() */
    virtual void parse_options(void);



  private:

    //! The relaxation factor to be used
    double _relax;


};


//
// inline methods
//

inline
RelaxationMethod::RelaxationMethod(const ModelOptions& options)
  : SelfconsistentSolver(options),
    _relax(1.0)
{
}

inline
RelaxationMethod::~RelaxationMethod(void)
{
}

inline
RelaxationMethod*
RelaxationMethod::create(const ModelOptions& options)
{
  return new RelaxationMethod(options);
}



#endif // _RELAXATIONMETHOD_H_
