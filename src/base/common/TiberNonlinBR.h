// $Id$

#ifndef _TIBERNONLINBR_H_
#define _TIBERNONLINBR_H_


#include "TiberLineSearch.h"


class TiberLinearSolver;


//! An implementation of line search to solve nonlinear systems
class TBDLLOCAL TiberNonlinBR : public TiberLineSearch
{

  public:

    //! Constructor
    TiberNonlinBR(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number);

    //! Destructor
    virtual ~TiberNonlinBR(void);


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
TiberNonlinBR::system_type(void) const
{
  return "TiberNonlinBR";
}


#endif // _TIBERNONLINBR_H_
