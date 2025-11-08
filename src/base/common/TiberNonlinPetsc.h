// $Id$

#ifndef _TIBERNONLINPETSC_H_
#define _TIBERNONLINPETSC_H_


#include "TiberNonlinearSystem.h"


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
