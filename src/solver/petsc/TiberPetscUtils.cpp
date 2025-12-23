/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file TiberPetscUtils.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */


#include "TiberPetscUtils.h"
#include "tibercad/base/ModelOptions.h"


#include <string>
#include <iostream>


const KSPType
TiberPetscUtils::extract_KSPType(const ModelOptions& options)
{
  const std::string& ksp(options.get_option("method", "bcgs"));
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
  else if (ksp == "lsqr")
    return KSPLSQR;
  else if (ksp == "pconly")
    return KSPPREONLY;
  else 
  {
    std::cerr << "PETSc: unknown Krylov method \'"
      << ksp << "\'. Falling back to \'bcgs\'" << std::endl;
  }

  return KSPBCGS;
}


const PCType
TiberPetscUtils::extract_PCType(const ModelOptions& options)
{
  const std::string& pc(options.get_option("preconditioner", "ilu"));
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


