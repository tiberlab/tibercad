// $Id$

#ifndef _TIBERNONLINLS_H_
#define _TIBERNONLINLS_H_


#include "TiberNonlinearSystem.h"


template<typename> class TiberPetscLinearSolver;


//! An implementation of line search to solve nonlinear systems
class TiberNonlinLS : public TiberNonlinearSystem
{

  public:

    //! Constructor
    TiberNonlinLS(EquationSystems& es,
        const std::string& name,
        const unsigned int number);

    //! Destructor
    virtual ~TiberNonlinLS(void);


    /*! \copydoc TiberNonlinearSystem::clear() */
    virtual void clear(void);


    /*! \copydoc TiberNonlinearSystem::reinit() */
    virtual void reinit(void);


    /*! \copydoc TiberNonlinearSystem::solve() */
    virtual void solve(void);


    /*! \copydoc System:user_initialization() */
    virtual void user_initialization(void);


    /*! \copydoc TiberNonlinearSystem::system_type() */
    virtual std::string system_type(void) const;


    //! Get the solution vector
    virtual NumericVector<double>& get_solution_vector(void);



  private:

    //! The parent class type
    typedef TiberNonlinearSystem Parent;

    //! The linear solver to be used for the Newton iteration
    TiberPetscLinearSolver<double>* _solver;

};



//
// inline methods
//


inline
std::string
TiberNonlinLS::system_type(void) const
{
  return "TiberNonlinLS";
}



inline
NumericVector<double>&
TiberNonlinLS::get_solution_vector(void)
{
  return get_vector("sol");
}


#endif // _TIBERNONLINLS_H_
