// $Id$

#ifndef _TIBERNONLINTR_H_
#define _TIBERNONLINTR_H_


#include "TiberLineSearch.h"


class TiberLinearSolver;


//! An implementation of line search to solve nonlinear systems
class TBDLLOCAL TiberNonlinTR : public TiberLineSearch
{

  public:

    //! Constructor
    TiberNonlinTR(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number);

    //! Destructor
    virtual ~TiberNonlinTR(void);


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
TiberNonlinTR::system_type(void) const
{
  return "TiberNonlinTR";
}


#endif // _TIBERNONLINTR_H_
