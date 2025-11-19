// $Id$


#ifndef _PARDISOSOLVEREXCEPTION_H_
#define _PARDISOSOLVEREXCEPTION_H_

#include "tibercad/solver/SolverException.h"

#include <stdexcept>
#include <string>

//! An exception class for the PARDISO solver
class PardisoSolverException : public SolverException
{

 public:
  
  PardisoSolverException(int error)
    : SolverException("Error in Pardiso solver") { };
    


  private:

};





#endif // _PARDISOSOLVEREXCEPTION_H_
