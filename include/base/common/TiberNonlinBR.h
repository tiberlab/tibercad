// $Id$

#ifndef _TIBERNONLINBR_H_
#define _TIBERNONLINBR_H_


#include "TiberNonlinearSystem.h"


class TiberLinearSolver;


//! An implementation of line search to solve nonlinear systems
class TiberNonlinBR : public TiberNonlinearSystem
{

  public:

    //! Constructor
    TiberNonlinBR(EquationSystems& es,
        const std::string& name,
        const unsigned int number);

    //! Destructor
    virtual ~TiberNonlinBR(void);


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
    TiberLinearSolver* _solver;

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



inline
NumericVector<double>&
TiberNonlinBR::get_solution_vector(void)
{
  return get_vector("sol");
}


#endif // _TIBERNONLINBR_H_
