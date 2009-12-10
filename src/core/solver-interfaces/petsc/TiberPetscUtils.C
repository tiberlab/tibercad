// $Id$

#include "TiberPetscUtils.h"
#include "ModelOptions.h"


#include <string>
#include <iostream>


KSPType
TiberPetscUtils::extract_KSPType(const ModelOptions& options)
{
  const std::string& ksp(options.get_option("ksp_type", "bcgs"));
  if (ksp == "bcgs")
    return KSPBCGS;
  else if (ksp == "bcgsl")
    return KSPBCGSL;
  else if (ksp == "gmres")
    return KSPGMRES;
  else if (ksp == "bicg")
    return KSPBICG;
  else if (ksp == "cg")
    return KSPCG;
  else if (ksp == "cgs")
    return KSPCGS;
  else if (ksp == "richardson")
    return KSPRICHARDSON;
  else if (ksp == "pconly")
    return KSPPREONLY;
  else 
  {
    std::cerr << "PETSc: unknown Krylov method \'"
      << ksp << "\'. Falling back to \'bcgs\'" << std::endl;
  }

  return KSPBCGS;
}


PCType
TiberPetscUtils::extract_PCType(const ModelOptions& options)
{
  const std::string& pc(options.get_option("pc_type", "ilu"));
  if (pc == "ilu")
    return PCILU;
  else if (pc == "lu")
    return PCLU;
  else if (pc == "jacobi")
    return PCJACOBI;
  else if (pc == "composite")
    return PCCOMPOSITE;
  else if (pc == "none")
    return PCNONE;
  else if (pc == "cholesky")
    return PCCHOLESKY;
  else
  {
    std::cerr << "PETSc: unknown preconditioner \'"
      << pc << "\'. Falling back to \'ilu\'" << std::endl;
  }

  return PCILU;
}



int
TiberPetscUtils::extract_LSType(const ModelOptions& options)
{
  int ls_type = 3;
  const std::string lstype(options.get_option("ls_type", "cubic"));
  if (lstype == "cubic")
    ls_type = 3;
  else if (lstype == "none")
    ls_type = 1;
  else if (lstype == "quadratic")
    ls_type = 2;
  else
  {
    std::cerr << "PETSc: unknown linesearch \'" << lstype
      << "\', falling back to \'cubic\'." << std::endl;
  }

  return ls_type;
}


void
TiberPetscUtils::_checkerr(int errorcode, int line, const char* file)
{
  if (errorcode != 0)
  {
    std::cerr << "Error in " << file << ", line " << line << ":\n" << std::flush;
    throw PetscRuntimeError(errorcode);
  }
}


