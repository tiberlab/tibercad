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
 * \file TiberNonlinPetsc.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef _TIBERNONLINPETSC_H_
#define _TIBERNONLINPETSC_H_


#include "tibercad/solver/TiberNonlinearSystem.h"


class TiberPetscNonlinearSolver;


//! An implementation of line search to solve nonlinear systems
class TBDLLOCAL TiberNonlinPetsc : public TiberNonlinearSystem
{

  public:

    //! Constructor
    TiberNonlinPetsc(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number);

    //! Destructor
    virtual ~TiberNonlinPetsc(void);


    /*! \copydoc TiberNonlinearSystem::clear() */
    virtual void clear(void);


    /*! \copydoc TiberNonlinearSystem::reinit() */
    virtual void reinit(void);


    /*! \copydoc TiberNonlinearSystem::system_type() */
    virtual std::string system_type(void) const;


    //! Get the solution vector
    virtual libMesh::NumericVector<double>& get_solution_vector(void);



  protected:

    /*! \copydoc TiberNonlinearSystem::do_solve() */
    virtual void do_solve(void);


    
  private:

    //! The parent class type
    typedef TiberNonlinearSystem Parent;

    //! The nonlinear solver to be used
    TiberPetscNonlinearSolver* _solver;

};



//
// inline methods
//


inline
std::string
TiberNonlinPetsc::system_type(void) const
{
  return "TiberNonlinPetsc";
}


inline
libMesh::NumericVector<double>&
TiberNonlinPetsc::get_solution_vector(void)
{
  return *(solution);
}



#endif // _TIBERNONLINPETSC_H_
