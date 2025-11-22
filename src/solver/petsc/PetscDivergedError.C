// $Id$

#include "PetscDivergedError.h"

#include <sstream>

using namespace std;

PetscDivergedError::PetscDivergedError(int reason, int iteration,
    double fnorm, PETScSolverType type)
  : PetscRuntimeError(PETSC_ERR_CONV_FAILED),
    _iteration(iteration), _fnorm(fnorm),
    _reason(reason), _type(type)
{
  ostringstream os;
  switch (type)
  {
    case KSP:
      os << "Linear PETSc solver (KSP) error in iteration "
         << iteration << ": ";
      switch (reason)
      {
        case KSP_DIVERGED_ITS:
          os << "reached maximum number of iterations";
          break;

        case KSP_DIVERGED_DTOL:
          os << "residual norm increased by a factor of divtol";
          break;

        case KSP_DIVERGED_NANORINF:
          os << "residual norm became Not-a-number or Inf";
          break;

        case KSP_DIVERGED_BREAKDOWN:
          os << "generic breakdown in solver method";
          break;

        case KSP_DIVERGED_BREAKDOWN_BICG:
          os << "Initial residual is orthogonal to preconditioned initial residual";
          break;

        default:
          os << "unknown error " << reason;
          break;

      }
      break;

    case SNES:
      os << "Nonlinear PETSc solver (SNES) error in iteration "
         << iteration << ": ";

      switch (reason)
      {
        case SNES_DIVERGED_FUNCTION_DOMAIN:
          os << "new trial solution is utside the valid domain";
          break;

        case SNES_DIVERGED_FUNCTION_COUNT:
          os << "reached maximum allowed number of evaluations of residual";
          break;

        case SNES_DIVERGED_LINEAR_SOLVE:
          os << "the linear solver failed";
          break;

        case SNES_DIVERGED_FNORM_NAN:
          os << "the resiual 2-norm produces Not-a-number";
          break;

        case SNES_DIVERGED_MAX_IT:
          os << "reached the maximum number of iterations requested";
          break;

        case SNES_DIVERGED_LINE_SEARCH:
          os << "the line search failed";
          break;

        case SNES_DIVERGED_INNER:
          os << "the inner solve failed";
          break;

        case SNES_DIVERGED_LOCAL_MIN:
          os << "the algorithm seems to have stagnated at a local minimum that is not zero";
          break;

        /*
        case SNES_DIVERGED_DTOL:
          os << "the norm of the residual increased by more than the allowed factor";
          break;

        case SNES_DIVERGED_JACOBIAN_DOMAIN:
          os << "Jacobian cannot be calculated";
          break;

        case SNES_DIVERGED_TR_DELTA:
          os << "trust region method diverged";
          break;
        */

        default:
          os << "unknown error " << reason;
          break;
      }
      break;

    default:
      os << "Error in unknown PETSc solver.";
      break;
  }

  set_message(os.str());
}
