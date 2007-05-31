// $Id$

#ifndef _TIBERNONLINPETSC_H_
#define _TIBERNONLINPETSC_H_


#include "TiberNonlinearSystem.h"


template<typename> class TiberPetscNonlinearSolver;


//! An implementation of line search to solve nonlinear systems
class TiberNonlinPetsc : public TiberNonlinearSystem
{

  public:

    //! Constructor
    TiberNonlinPetsc(EquationSystems& es,
        const std::string& name,
        const unsigned int number);

    //! Destructor
    virtual ~TiberNonlinPetsc(void);


    /*! \copydoc TiberNonlinearSystem::clear() */
    virtual void clear(void);


    /*! \copydoc TiberNonlinearSystem::reinit() */
    virtual void reinit(void);


    /*! \copydoc TiberNonlinearSystem::solve() */
    virtual void solve(void);


    /*! \copydoc TiberNonlinearSystem::system_type() */
    virtual std::string system_type(void) const;


    //! Get the solution vector
    virtual NumericVector<double>& get_solution_vector(void);



  private:

    //! The parent class type
    typedef TiberNonlinearSystem Parent;

    //! The nonlinear solver to be used
    TiberPetscNonlinearSolver<double>* _solver;

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
NumericVector<double>&
TiberNonlinPetsc::get_solution_vector(void)
{
  return *(solution);
}



#endif // _TIBERNONLINPETSC_H_
