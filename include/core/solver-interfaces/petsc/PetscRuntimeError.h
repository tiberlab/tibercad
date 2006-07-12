// $Id$

#ifndef _PETSCRUNTIMEERROR_H_
#define _PETSCRUNTIMEERROR_H_

#include <stdexcept>

class PetscRuntimeError : public std::runtime_error
{

  public:
    PetscRuntimeError(int reason)
      : std::runtime_error("Internal PETSc error."),
        _reason(reason) {};

    int get_reason(void) { return _reason; };

  private:

    int _reason;
};


#endif // _PETSCRUNTIMEERROR_H_
