// $Id$

#include "tiber_config.h"
#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"
#ifdef ENABLE_PARDISO
# include "PardisoLinearSolver.h"
#endif
#include "InitFailedException.h"


TiberLinearSolver*
TiberLinearSolver::create(const std::string& type)
{
  TiberLinearSolver* solver = NULL;

  if (type == "petsc")
    solver = new TiberPetscLinearSolver();
#ifdef ENABLE_PARDISO
  else if (type == "pardiso")
    solver = new PardisoLinearSolver();
#endif
  
  if (solver == NULL)
  {
    std::string msg = "TiberLinearSolver: no such solver ";
    msg += type;
    throw InitFailedException(msg);
  }
#ifdef DEBUG
  std::cerr << "Created linear solver type" << type << std::endl;
#endif

  return solver;
}
        
