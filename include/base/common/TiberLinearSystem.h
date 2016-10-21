// $Id$

#ifndef _TIBERLINEARSYSTEM_H_
#define _TIBERLINEARSYSTEM_H_

#include "TiberEqSystem.h"
#include "libMeshDefs.h"

#include "linear_implicit_system.h"


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


    /*! \copydoc ImplicitSystem::clear() */
    //virtual void clear(void);


    /*! \copydoc ImplicitSystem::reinit() */
    //virtual void reinit(void);


    /*! \copydoc ImplicitSystem::solve() */
    virtual void solve(void);


    /*! \copydoc ImplicitSystem::system_type() */
    virtual std::string system_type(void) const;


    /*! \copydoc System:user_initialization() */
    virtual void user_initialization(void);


    //! Get the solution vector
    virtual libMesh::NumericVector<double>& get_solution_vector(void);



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


#endif // _TIBERLINEARSYSTEM_H_
