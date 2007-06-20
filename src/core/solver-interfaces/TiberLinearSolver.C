// $Id$


#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"
#include "PardisoLinearSolver.h"
#include "InitFailedException.h"


TiberLinearSolver*
TiberLinearSolver::create(const std::string& type)
{
  TiberLinearSolver* solver = NULL;

  if (type == "petsc")
    solver = new TiberPetscLinearSolver();
  else if (type == "pardiso")
    solver = new PardisoLinearSolver();
    //solver = NULL;
  
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
        
