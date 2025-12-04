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
 * \file TiberLinearSystem.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _TIBERLINEARSYSTEM_H_
#define _TIBERLINEARSYSTEM_H_

#include "tibercad/solver/TiberEqSystem.h"
#include "tibercad/base/libMeshDefs.h"

#include "libmesh/linear_implicit_system.h"


class TiberLinearSolver;



//! A base class for linear systems
class TiberLinearSystem : public TiberEqSystem, public libMesh::LinearImplicitSystem
{

  public:

    //! Constructor
    TiberLinearSystem(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number);


    //! Destructor
    virtual ~TiberLinearSystem(void) { };


    //! Create a linear system
    /*!
     * \param es the EquationSystems object where the new system will be added
     * \param sysname the name of the new system
     * \param options the options for this system
     * \return a pointer to the newly created system
     */
    static TiberLinearSystem* create(libMesh::EquationSystems& es,
        const std::string& sysname, const ModelOptions& options);


    /*! \copydoc ImplicitSystem::solve() */
    virtual void solve(void);


    /*! \copydoc ImplicitSystem::system_type() */
    virtual std::string system_type(void) const;


    /*! \copydoc System:user_initialization() */
    virtual void user_initialization(void);


    //! Get the solution vector including ghost values
    virtual libMesh::NumericVector<double>& get_solution_vector(void);


    //! Get the local solution vector without ghost values
    virtual libMesh::NumericVector<double>& get_local_solution_vector(void);



  private:


};
    

//
// inline methods
// 

inline
std::string
TiberLinearSystem::system_type(void) const
{
  return "TiberLinear";
}


inline
libMesh::NumericVector<double>&
TiberLinearSystem::get_solution_vector(void)
{
  return *current_local_solution;
}


inline
libMesh::NumericVector<double>&
TiberLinearSystem::get_local_solution_vector(void)
{
  return *solution;
}


#endif // _TIBERLINEARSYSTEM_H_
