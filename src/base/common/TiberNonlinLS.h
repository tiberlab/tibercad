// $Id$

#ifndef _TIBERNONLINLS_H_
#define _TIBERNONLINLS_H_


#include "TiberLineSearch.h"


class TiberLinearSolver;


//! An implementation of line search to solve nonlinear systems
class TBDLLOCAL TiberNonlinLS : public TiberLineSearch
{

  public:

    //! Constructor
    TiberNonlinLS(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number);

    //! Destructor
    virtual ~TiberNonlinLS(void);


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
TiberNonlinLS::system_type(void) const
{
  return "TiberNonlinLS";
}


#endif // _TIBERNONLINLS_H_
