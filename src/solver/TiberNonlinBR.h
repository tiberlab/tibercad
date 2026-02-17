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
 * \file TiberNonlinBR.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef TC_TIBERNONLINBR_H
#define TC_TIBERNONLINBR_H


#include "tibercad/solver/TiberLineSearch.h"


class TiberLinearSolver;


//! An implementation of line search to solve nonlinear systems
class TBDLLOCAL TiberNonlinBR : public TiberLineSearch
{

  public:

    //! Constructor
    TiberNonlinBR(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number);

    //! Destructor
    virtual ~TiberNonlinBR(void);


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
TiberNonlinBR::system_type(void) const
{
  return "TiberNonlinBR";
}


#endif // TC_TIBERNONLINBR_H
