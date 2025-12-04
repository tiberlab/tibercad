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
 * \file TiberNonlinTR.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef _TIBERNONLINTR_H_
#define _TIBERNONLINTR_H_


#include "tibercad/solver/TiberLineSearch.h"


class TiberLinearSolver;


//! An implementation of line search to solve nonlinear systems
class TBDLLOCAL TiberNonlinTR : public TiberLineSearch
{

  public:

    //! Constructor
    TiberNonlinTR(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number);

    //! Destructor
    virtual ~TiberNonlinTR(void);


    /*! \copydoc TiberNonlinearSystem::system_type() */
    virtual std::string system_type(void) const;


  protected:


    /*! \copydoc TiberNonlinearSystem::do_solve() */
    virtual void do_solve(void);



  private:

    //! The parent class type
    typedef TiberLineSearch Parent;

};



//
// inline methods
//


inline
std::string
TiberNonlinTR::system_type(void) const
{
  return "TiberNonlinTR";
}


#endif // _TIBERNONLINTR_H_
