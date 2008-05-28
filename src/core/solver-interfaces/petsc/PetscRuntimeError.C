// $Id$

#include "PetscRuntimeError.h"

#include <sstream>

const char*
PetscRuntimeError::what(void) const throw()
{
  std::ostringstream os;
  os << "Internal PETSc error : " << get_reason();
  return os.str().c_str();
}
