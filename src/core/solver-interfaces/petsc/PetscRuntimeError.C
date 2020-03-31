// $Id$

#include "PetscRuntimeError.h"

#include <sstream>


PetscRuntimeError::PetscRuntimeError(int reason)
  : SolverException("Internal PETSc error."),
    _reason(reason)
{
  std::ostringstream os;
  os << "Internal PETSc error: " << get_reason();
  _msg = os.str();
}


const char*
PetscRuntimeError::what(void) const throw()
{
  return _msg.c_str();
}


void
PetscRuntimeError::set_message(const std::string& msg)
{
  _msg = msg;
}
