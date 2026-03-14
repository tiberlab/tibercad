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
 * \file TiberTransientSystem.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_TIBERTRANSIENTSYSTEM_H
#define TC_TIBERTRANSIENTSYSTEM_H

#include "tibercad/solver/TiberEqSystem.h"
#include "tibercad/base/libMeshDefs.h"

#include "libmesh/transient_system.h"




/*!
 * \brief First order transient system in tibercad
 *
 * This system is used to describe/solve a first order ODE of type
 * a*\partial u \partial t + Lu = f
 * 
 * The following methods are implemented:
 * - Forward Euler
 * - Backward Euler
 * - trapezoidal
 * - adaptive time step using trapezoidal and backward Euler
 * 
 * Matrices and vectors are assembled on the final time coordinate in all
 * cases for now.
 */
template <class Base>
class TiberTransientSystem : public TiberEqSystem, public libMesh::TransientSystem<Base>
{

  public:

    //! Constructor
    TiberTransientSystem(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number);


    //! Destructor
    virtual ~TiberTransientSystem(void);


    //! Create a linear system
    /*!
     * \param es the EquationSystems object where the new system will be added
     * \param sysname the name of the new system
     * \param options the options for this system
     * \return a pointer to the newly created system
     */
    static TiberTransientSystem* create(libMesh::EquationSystems& es,
        const std::string& sysname, const ModelOptions& options);


    /*! \copydoc ImplicitSystem::solve() */
    virtual void solve(void) override;


    /*! \copydoc ImplicitSystem::system_type() */
    virtual std::string system_type(void) const override;


    /*! \copydoc System:user_initialization() */
    virtual void user_initialization(void) override;


    //! Get the solution vector including ghost values
    virtual libMesh::NumericVector<double>& get_solution_vector(void) override;


    //! Get the local solution vector without ghost values
    virtual libMesh::NumericVector<double>& get_local_solution_vector(void) override;


    //! Set the target time
    void set_target_time(double time);


    //! Pass a set with Dirichlet DoF IDs
    /*!
     * We use this information in the explicit Euler implementation,
     * in order to set the residual to 0, so that the boundary nodes
     * do not change their value even for large timesteps
     */
    void set_dirichlet_dofs(const std::set<libMesh::dof_id_type>& dirichlet_dofs);


  private:

    typedef libMesh::TransientSystem<Base> parent_type;

    //! Forward Euler time stepping
    void forward_euler(void);

    //! Backward Euler time stepping
    void backward_euler(void);

    //! Backward Euler time stepping
    void trapezoidal(void);

    //! Backward Euler time stepping
    void adaptive(void);

    //! The next target time
    double _target_time = 0.;

    //! The current time step
    double _time_step = 1e12;

    //! The minimum allowed time step
    double _min_time_step = 1e-12;

    //! The Dirichlet DoFs
    std::set<libMesh::dof_id_type> _dirichlet_dofs;

};
    

typedef TiberTransientSystem<libMesh::LinearImplicitSystem> TiberTransientLinSystem;

#endif // TC_TIBERTRANSIENTSYSTEM_H
