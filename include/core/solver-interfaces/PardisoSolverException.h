// $Id: LinearSolverException.h 365 2007-06-11 16:15:32Z maufder $


#ifndef _PARDISOSOLVEREXCEPTION_H_
#define _PARDISOSOLVEREXCEPTION_H_

#include "LinearSolverException.h"

#include <stdexcept>
#include <string>

//! An exception class for the linear solver interfaces
class PardisoSolverException : public LinearSolverException
{

 public:
  
  PardisoSolverException(int error)
    :LinearSolverException("Error in Pardiso solver"){};
    



  private:

};





#endif // _LINEARSOLVEREXCEPTION_H_
