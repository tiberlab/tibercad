// $Id$

#ifndef _PETSCRUNTIMEERROR_H_
#define _PETSCRUNTIMEERROR_H_

#include "tibercad/solver/SolverException.h"

#include <string>

class PetscRuntimeError : public SolverException
{

  public:

    PetscRuntimeError(int reason);

    virtual ~PetscRuntimeError(void) throw() {};

    int get_reason(void) const { return _reason; };

    virtual const char* what(void) const throw();

  protected:

    void set_message(const std::string& msg);

  private:

    std::string _msg;

    int _reason;
};


#endif // _PETSCRUNTIMEERROR_H_
