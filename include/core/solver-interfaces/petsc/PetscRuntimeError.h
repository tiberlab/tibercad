// $Id$

#ifndef _PETSCRUNTIMEERROR_H_
#define _PETSCRUNTIMEERROR_H_

#include "SolverException.h"

class PetscRuntimeError : public SolverException
{

  public:

    PetscRuntimeError(int reason)
      : SolverException("Internal PETSc error."),
        _reason(reason) {};

    int get_reason(void) const { return _reason; };

    virtual const char* what(void) const throw();


  private:

    int _reason;
};


#endif // _PETSCRUNTIMEERROR_H_
