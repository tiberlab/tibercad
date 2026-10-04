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
 * \file EigenSolver.C
 * \brief tiberCAD API implementation.
 */


#include <iostream>
#include <cassert>
#include <string>

#include "tibercad/solver/slepc/EigenSolver.h"
#include "tibercad/io/Messages.h"
#include "tibercad/base/RuntimeException.h"
#include "solver/petsc/TiberPetscUtils.h"

#include "slepceps.h"

// Private PETSc dense-matrix header: needed only to force-clear the
// matinuse flag on DS matrices that SLEPc leaves "checked out" when
// KSPSolve diverges mid-Krylov iteration (EPSDestroy would crash otherwise).
// This is an internal workaround tied to PETSc 3.23 + SLEPc 3.23.
#include "../../../external/petsc-3.23.4/src/mat/impls/dense/seq/dense.h"

// Private SLEPc DS header: needed to access omat[] array in _p_DS struct.
#include "../../../external/slepc-3.23.1/include/slepc/private/dsimpl.h"


#include "petsc/private/matimpl.h"


using namespace std;

namespace
{
  Mat A; // Hamiltonian
  Mat B; // S-matrix
  Mat M; // Shell matrix for spectrum folding
  Mat P; // preconditioner matrix
  EPS eps; // EigenSolver
  MPI_Comm slepc_comm;
  // Set to true the first time a retry-with-LU succeeds, so that subsequent
  // calls in the same k-point loop also use LU directly instead of the
  // iterative KSP that already proved unreliable for this problem.
  bool lu_fallback_active = false;
  double shift; // could be stored in ST but lapack does not apply any shift
  
  vector<Vec> _deflation_space;

  // A context structure for spectrum folding
  typedef struct
  {
    Mat       A;
    Vec       w;
    PetscReal target;
  } STFoldCtx;

  // An instance of STFoldCtx
  STFoldCtx *fold_ctx = nullptr;

  // The spectrum folding
  PetscErrorCode MatMult_Fold(Mat M, Vec x, Vec y)
  {
    STFoldCtx   *ctx;
    PetscScalar    sigma;

    MatShellGetContext(M, &ctx);
    sigma = -ctx->target;
    MatMult(ctx->A, x, ctx->w);
    VecAXPY(ctx->w, sigma,x);
    MatMult(ctx->A, ctx->w,y);
    VecAXPY(y,sigma, ctx->w);

    return(0);
  }


  // For spectrum folding, we need to recover the real eigenvalues
  PetscErrorCode RayleighQuotient(Mat A, Vec x, PetscScalar *r)
  {
    Vec Ax;

    VecDuplicate(x, &Ax);
    MatMult(A, x, Ax);
    VecDot(Ax, x, r);
    VecDestroy(&Ax);

    return(0);
  }
}


static int set_ksp_and_pc(ST st, const EigenSolver::SLEPCoptions& opt);

static void set_sub_pc(PC pc, PCType pc_type);

// ---------------------------------------------------------------------------
// Force-release any DS matrices that SLEPc left "checked out" (matinuse != 0)
// after a failed Krylov iteration.  KSPSolve divergence causes BVMatLanczos /
// BVMatArnoldi to unwind without calling DSRestoreMat, so the underlying
// Mat_SeqDense objects have matinuse set.  EPSDestroy -> DSDestroy ->
// DSReset -> MatDestroy then crashes with "Need to call
// MatDenseRestoreSubMatrix() first".
//
// We iterate over all DS matrix slots and, for any that have matinuse != 0,
// directly zero the flag and null out the pointer – the same minimal work that
// MatDenseRestoreSubMatrix_SeqDense() does, minus the array-reset step that
// isn't needed here since we are about to destroy the whole DS anyway.
//
// This is intentionally a private, fire-and-forget helper: it is only called
// immediately before EPSDestroy inside the retry path, so safety is not a
// concern for any subsequent use of the same objects.
// ---------------------------------------------------------------------------
static void force_release_ds_matrices(EPS eps)
{
  DS ds;
  EPSGetDS(eps, &ds);
  if (!ds) return;

  for (int m = 0; m < DS_NUM_MAT; ++m)
  {
    // Access omat[] through the public DS struct layout. DS is typedef'd as
    // struct _p_DS* (defined in slepc/private/dsimpl.h).
    Mat omat = ds->omat[m];
    if (!omat) continue;

    Mat_SeqDense *a = (Mat_SeqDense *)omat->data;
    if (!a) continue;

    if (a->matinuse)
    {
      // Zero the flag – this is exactly what MatDenseRestoreSubMatrix does
      // after resetting the column array. We skip MatDenseResetArray because
      // we are about to destroy everything; skipping it is safe.
      a->matinuse = 0;
      a->cmat     = nullptr;
    }
  }
}


int EigenSolver::_size_of_matrix;
//-------------------------------------------------------------//
void EigenSolver::slepc_init(int argc1, char** argv1, MPI_Comm comm)
{

  slepc_comm = comm;

  SlepcInitialize(&argc1,&argv1,NULL,NULL);
  PetscPopSignalHandler();


}

//--------------------------------------------------------------//
void  EigenSolver::slepc_done()
{
 SlepcFinalize();
}


//--------------------------------------------------------------//
int EigenSolver::eig_value_problem_general(const EigenSolver::SLEPCoptions& opt)
{
  return(eig_value_problem(opt, GENERALIZED));
}




int EigenSolver::eig_value_problem(const EigenSolver::SLEPCoptions& opt,
                                   EigenSolver::EVPType evp_type)
{

  EPSType     type;
  PetscReal   error, tol, re, im;
  PetscScalar kr, ki;
  int         nev, ierr, maxit, i, its, lits, nconv;
  char        filename[256];
  PetscViewer viewer, viewer_out, viewer_eigvals;
  //PetscBool  flg;
  //PetscMPIInt    rank,size;
  ST st;
  KSP ksp;
  PC pc;



  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //print_options(opt);


  //ierr = MPI_Comm_size(slepc_comm,&size);TiberPetscUtils::checkerr(ierr);
  //ierr = MPI_Comm_rank(slepc_comm,&rank);TiberPetscUtils::checkerr(ierr);

  if (opt.read_matrix_from_file)
  {
    ierr = PetscViewerBinaryOpen(slepc_comm,opt.H_file_name.c_str(),FILE_MODE_READ,&viewer);TiberPetscUtils::checkerr(ierr); //their

    //ierr = MatLoad(viewer,MATAIJ,&A);TiberPetscUtils::checkerr(ierr);
    ierr = MatLoad(A,viewer);TiberPetscUtils::checkerr(ierr);
    ierr = PetscViewerDestroy(&viewer);TiberPetscUtils::checkerr(ierr);

    ierr = MatGetSize(A, &_size_of_matrix, NULL);

    if (evp_type == GENERALIZED)
    {
      ierr = PetscViewerBinaryOpen(slepc_comm,opt.S_file_name.c_str(),FILE_MODE_READ,&viewer);TiberPetscUtils::checkerr(ierr); //their

      //ierr = MatLoad(viewer,MATAIJ,&B);TiberPetscUtils::checkerr(ierr);
      ierr = MatLoad(B,viewer);TiberPetscUtils::checkerr(ierr);
      ierr = PetscViewerDestroy(&viewer);TiberPetscUtils::checkerr(ierr);
    }
  }



  if (opt.matrix_output)
  {

    ierr = PetscViewerASCIIOpen(slepc_comm,"matA.m",&viewer_out); TiberPetscUtils::checkerr(ierr);
    ierr = PetscViewerPushFormat(viewer_out,PETSC_VIEWER_ASCII_MATLAB);
    ierr = MatView(A, viewer_out); TiberPetscUtils::checkerr(ierr);
    ierr = PetscViewerPopFormat(viewer_out);
    ierr = PetscViewerDestroy(&viewer_out);TiberPetscUtils::checkerr(ierr);


    if (evp_type == GENERALIZED)
    {
      ierr = PetscViewerASCIIOpen(slepc_comm,"matB.m",&viewer_out); TiberPetscUtils::checkerr(ierr);
      ierr = PetscViewerPushFormat(viewer_out,PETSC_VIEWER_ASCII_MATLAB);
      ierr = MatView(B, viewer_out); TiberPetscUtils::checkerr(ierr);
      ierr = PetscViewerPopFormat(viewer_out);
      ierr = PetscViewerDestroy(&viewer_out);TiberPetscUtils::checkerr(ierr);
    }
  }


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  shift = opt.spectrum_shift.real();

  if (opt.spectral_trans == "folding")
  {
    if (fold_ctx == nullptr)
    {
      PetscInt N, nloc, mloc;
      MatGetSize(A, &N, nullptr); // we only have square matrices
      MatGetLocalSize(A, &nloc, &mloc);

      PetscNew(&fold_ctx);
      fold_ctx->A = A;
      MatCreateVecs(A, &fold_ctx->w, nullptr);
      //VecDuplicate(x, &ctx->w);

      MatCreateShell(slepc_comm, nloc, mloc, N, N, fold_ctx, &M);
      MatShellSetOperation(M, MATOP_MULT, (void(*)(void))MatMult_Fold);
    }

    fold_ctx->target = shift;

  }

  // Set operators.
  if (evp_type == GENERALIZED)
  {
    if (opt.spectral_trans == "folding")
      ierr = EPSSetOperators(eps, M, B);
    else
      ierr = EPSSetOperators(eps, A, B);

    TiberPetscUtils::checkerr(ierr);

    ierr = EPSSetProblemType(eps, EPS_GHEP);TiberPetscUtils::checkerr(ierr);
  }
  else
  {
    if (opt.spectral_trans == "folding")
      ierr = EPSSetOperators(eps, M, PETSC_NULLPTR);
    else
      ierr = EPSSetOperators(eps, A, PETSC_NULLPTR);

    TiberPetscUtils::checkerr(ierr);

    ierr = EPSSetProblemType(eps,EPS_HEP);TiberPetscUtils::checkerr(ierr);
  }



  ierr = EPSSetTolerances(eps, opt.eps_tolerance, opt.eps_max_it);  TiberPetscUtils::checkerr(ierr);

  if (opt.spectral_trans == "folding")
  {
    ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL); TiberPetscUtils::checkerr(ierr);
    ierr = EPSSetTarget(eps, 0.0);TiberPetscUtils::checkerr(ierr);
  }
  else
  {
    ierr = EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); TiberPetscUtils::checkerr(ierr);
    //ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);TiberPetscUtils::checkerr(ierr);
    ierr = EPSSetTarget(eps, opt.spectrum_shift);TiberPetscUtils::checkerr(ierr);
  }


  if (opt.solver_type == "arnoldi" || opt.solver_type == "krylovschur")
  {
    //ierr = EPSSetProblemType(eps,EPS_GNHEP);TiberPetscUtils::checkerr(ierr);

    if (opt.solver_type == "arnoldi")
      ierr = EPSSetType(eps, EPSARNOLDI);
    else if (opt.solver_type == "krylovschur")
      ierr = EPSSetType(eps, EPSKRYLOVSCHUR);


    ierr = EPSGetST(eps,&st); TiberPetscUtils::checkerr(ierr);
    
   
    if (opt.spectral_trans == "folding")
    {
      ierr = STSetType(st, STSHIFT); TiberPetscUtils::checkerr(ierr);
    }
    else
    {
      ierr = STSetType(st, STSINVERT); TiberPetscUtils::checkerr(ierr);
    }

    set_ksp_and_pc(st, opt);

  }
  else if (opt.solver_type == "lapack")
  {
    ierr = EPSSetType(eps, EPSLAPACK);
    //ierr = EPSGetST(eps,&st); TiberPetscUtils::checkerr(ierr);

  }
  else if (opt.solver_type == "arpack")
  {
    ierr = EPSSetType(eps, EPSARPACK);

    if (std::abs(opt.spectrum_shift) >1e-8)
    {
      ierr = EPSGetST(eps,&st); TiberPetscUtils::checkerr(ierr);
      ierr = STSetShift(st, opt.spectrum_shift);TiberPetscUtils::checkerr(ierr);

      ierr = STSetType(st,STSHIFT); TiberPetscUtils::checkerr(ierr);

      set_ksp_and_pc(st, opt);
    }
  }
  else if (opt.solver_type == "jd")
  {
    ierr = EPSSetType(eps, EPSJD);
    TiberPetscUtils::checkerr(ierr);

    // EPSJD applies its own preconditioner, so PETSc rejects every ST type
    // other than STPRECOND for it: the shift-and-invert spectral transform used
    // by the other EPS types is not available here. Build the diagonal-of-A^2
    // preconditioner matrix and hand it to the ST directly.
    PetscInt N;
    MatGetSize(A, &N, nullptr); // we only have square matrices

    MatCreate(PETSC_COMM_WORLD, &P);
    MatSetSizes(P, PETSC_DECIDE, PETSC_DECIDE, N, N);

    Vec dgv = nullptr;
    if (fold_ctx != nullptr)
      ierr = VecDuplicate(fold_ctx->w, &dgv);
    else
      ierr = MatCreateVecs(A, &dgv, nullptr);
    TiberPetscUtils::checkerr(ierr);

    // we abuse of the working vector in the folding context to calculate
    // the diagonal of A^2
    PetscInt start, stop;
    MatGetOwnershipRange(A, &start, &stop);

    for (PetscInt i = start; i < stop; ++i)
    {
      PetscInt nvals;
      const PetscScalar *vals;
      const PetscInt *cols;
      MatGetRow(A, i, &nvals, &cols, &vals);

      double sum = 0;
      for (PetscInt j = 0; j < nvals; ++j)
      {
        double norm = std::abs(vals[j]);
        //if (i == cols[j])
        //  norm -= shift;

        sum += norm*norm;
      }

      MatRestoreRow(A, i, &nvals, &cols, &vals);

      VecSetValue(dgv, i, sum, INSERT_VALUES);
    }
    VecAssemblyBegin(dgv);
    VecAssemblyEnd(dgv);

    MatSetUp(P);
    MatDiagonalSet(P, dgv, INSERT_VALUES);
    VecDestroy(&dgv);

    EPSGetST(eps, &st);
    STSetType(st, STPRECOND);
    STSetPreconditionerMat(st, P);
    MatDestroy(&P);
    //set_ksp_and_pc(st, opt);
    STGetKSP(st, &ksp);
    KSPSetType( ksp, KSPMINRES);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetTolerances(ksp,opt.spectrum_inversion_tolerance, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

    /*if (opt.monitor)
    {
      PetscViewerAndFormat *vf;
      ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT, &vf);
      TiberPetscUtils::checkerr(ierr);
      ierr = KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorResidual, vf, 0);
      TiberPetscUtils::checkerr(ierr);
    }*/
  }
  else if (opt.solver_type == "gd")
  {
    ierr = EPSSetType(eps, EPSGD);

    if ((opt.spectral_trans == "folding") && (fold_ctx != nullptr))
    {
      PetscInt N;
      MatGetSize(A, &N, nullptr); // we only have square matrices

      MatCreate(PETSC_COMM_WORLD, &P);
      MatSetSizes(P, PETSC_DECIDE, PETSC_DECIDE, N, N);

      // we abuse of the woking vector in the folding context to calculate
      // the diagonal of A^2
      PetscInt start, stop;
      MatGetOwnershipRange(A, &start, &stop);

      for (PetscInt i = start; i < stop; ++i)
      {
        PetscInt nvals;
        const PetscScalar *vals;
        const PetscInt *cols;
        MatGetRow(A, i, &nvals, &cols, &vals);

        double sum = 0;
        for (unsigned int j = 0; j < nvals; ++j)
        {
          double norm = std::abs(vals[j]);
          //if (i == cols[j])
          //  norm -= shift;

          sum += norm*norm;
        }

        MatRestoreRow(A, i, &nvals, &cols, &vals);

        VecSetValue(fold_ctx->w, i, sum, INSERT_VALUES);
      }
      VecAssemblyBegin(fold_ctx->w);
      VecAssemblyEnd(fold_ctx->w);


      MatSetUp(P);
      MatDiagonalSet(P, fold_ctx->w, INSERT_VALUES);

      EPSGetST(eps, &st);
      STSetType(st, STPRECOND);
      STSetPreconditionerMat(st, P);
      MatDestroy(&P);
      STGetKSP(st, &ksp);
      KSPSetType( ksp, KSPPREONLY);
      KSPGetPC(ksp, &pc);
      PCSetType(pc, PCJACOBI);
    }
  }
  else if (opt.solver_type == "feast")
  {
    ierr = EPSSetType(eps, EPSFEAST);
    TiberPetscUtils::checkerr(ierr);

    ierr = EPSSetInterval(eps, shift - 0.5, shift + 0.5);
    TiberPetscUtils::checkerr(ierr);

    ierr = EPSSetWhichEigenpairs(eps, EPS_ALL);
    TiberPetscUtils::checkerr(ierr);
  }
  else
  {
    throw RuntimeException("Unknown SLEPc solver: " + opt.solver_type);
  }


  ierr = do_solve(opt);
  TiberPetscUtils::checkerr(ierr);


  return ierr;
}

//------------------------------------------------------------------------------//

void EigenSolver::print_options(const EigenSolver::SLEPCoptions& opt)
{
  std::cout<<std::endl;
  std::cout<< "   (ES) read matrix from file: " << opt.read_matrix_from_file << std::endl;
  std::cout<< "   (ES) matrix output: " << opt.matrix_output << std::endl;
  std::cout<< "   (ES) ev_number: " << opt.ev_number << std::endl;
  std::cout<< "   (ES) solver type: " << opt.solver_type << std::endl;
  std::cout<< "   (ES) tolerance: " << opt.eps_tolerance << std::endl;
  std::cout<< "   (ES) max iterations: " << opt.eps_max_it << std::endl;
  std::cout<< "   (ES) spectral trans: " << opt.spectral_trans << std::endl;
  std::cout<< "   (ES) shift: " << opt.spectrum_shift << std::endl;
  std::cout<< "   (ES) linear system solver: " << opt.st_ksp_type << std::endl;
  std::cout<< "   (ES) preconditioner: " << opt.pc_type << std::endl;
  std::cout<< "   (ES) tolerance: " << opt.spectrum_inversion_tolerance << std::endl;

}




int set_ksp_and_pc(ST st, const EigenSolver::SLEPCoptions& opts)
{

  EigenSolver::SLEPCoptions opt(opts);

  // If a previous call already proved that the iterative KSP is unreliable
  // for this problem (lu_fallback_active was set by the retry path), upgrade
  // to a direct LU solve immediately instead of repeating the same failure.
  if (lu_fallback_active &&
      opt.solver_package != "mumps" &&
      opt.solver_package != "mkl_pardiso")
  {
    opt.pc_type      = "lu";
    opt.st_ksp_type  = "preonly";
  }

  int ierr;

  KSP ksp;
  ierr = STGetKSP(st, &ksp);

  PC pc;
  ierr = KSPGetPC(ksp, &pc);

#if ((SLEPC_VERSION_MAJOR < 3) || \
    ((SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR < 5)))
  PCSetOperators(pc, A, A, SAME_NONZERO_PATTERN);
#else
  PCSetOperators(pc, A, A);
#endif

  // if MUMPS or PARDISO is used as solver package, then we
  // want to use LU or Cholesky decomposition (direct solve through KSPPREONLY).
  // With solver_package=petsc the user keeps full control through the
  // ``pc_type`` / ``ksp_type`` options: the default (ilu + bcgsl) is an
  // iterative solve suited to large sparse systems, while ``pc_type = lu``
  // gives a direct solve (recommended for small/dense, e.g. coarse-grained,
  // Hamiltonians where ILU does not converge).
  if ((opt.solver_package == "mumps") ||
      (opt.solver_package == "mkl_pardiso"))
  {
    if (opt.pc_type != "cholesky")
      opt.pc_type = "lu";

    opt.st_ksp_type = "preonly";
  }

  if (opt.st_ksp_type == "bcgsl")
    ierr = KSPSetType( ksp, KSPBCGSL);
  else if (opt.st_ksp_type == "gmres" )
    ierr = KSPSetType( ksp, KSPGMRES);
  else if (opt.st_ksp_type == "bcgs" )
    ierr = KSPSetType( ksp, KSPBCGS);
  else if (opt.st_ksp_type == "cg" )
    ierr = KSPSetType( ksp, KSPCG);
  else if (opt.st_ksp_type == "richardson" )
    ierr = KSPSetType( ksp, KSPCG);
  else if (opt.st_ksp_type == "preonly")
    ierr = KSPSetType( ksp, KSPPREONLY);
  else
    throw RuntimeException("KSP type \'" + opt.st_ksp_type +
        "\' not supported in EigenSolver.");

  PetscMPIInt comm_size;
  MPI_Comm_size(slepc_comm, &comm_size);

  if (opt.pc_type == "cholesky")
    ierr = PCSetType(pc, PCCHOLESKY);
  else if (opt.pc_type == "jacobi" )
    ierr =  PCSetType(pc, PCJACOBI);
  else if (opt.pc_type == "icc" )
    ierr =  PCSetType(pc, PCICC);
  else if (opt.pc_type == "ilu" )
  {
    // in principle, this should have worked with this if, but for some reason even
    // when running with a single process it complains about mpiaij matrix type
    if (comm_size > 1)
    {
      ierr = PCSetType(pc, PCBJACOBI);
      ierr = PCSetUp(pc);
      set_sub_pc(pc, PCILU);
    }
    else
    {
      ierr = PCSetType(pc, PCILU);
      // ILU(0) with no reordering breaks down easily when the shift sigma is
      // close to an eigenvalue of H (which is exactly the shift-and-invert
      // scenario), because (H - sigma*I) is deliberately ill-conditioned and
      // the zero-fill ILU pivot can become near-zero.
      // ILU(1) with natural reordering is a cheap improvement: one level of
      // fill-in typically reduces the condition number of the preconditioned
      // system enough for BCGSL to converge reliably at modest extra cost.
      PCFactorSetLevels(pc, 1);
      PCFactorSetMatOrderingType(pc, MATORDERINGNATURAL);
    }
  }
  else if (opt.pc_type == "lu" )
  {
    // in principle, this should have worked with this if, but for some reason even
    // when running with a single process it complains about mpiaij matrix type
    if ((comm_size > 1) && (opt.solver_package != "mumps"))
    {
      ierr = PCSetType(pc, PCBJACOBI);
      ierr = PCSetUp(pc);
      set_sub_pc(pc, PCLU);
    }
    else
      ierr =  PCSetType(pc, PCLU);
  }
  else if (opt.pc_type == "redundant" )
    ierr =  PCSetType(pc,PCREDUNDANT);
  else if (opt.pc_type == "composite" )
    ierr =  PCSetType(pc,PCCOMPOSITE);
  else
    throw RuntimeException("Preconditioner \'" + opt.pc_type +
        "\' not supported in EigenSolver.");

  PCFactorSetMatSolverType(pc, opt.solver_package.c_str());


  ierr = KSPSetTolerances(ksp,opt.spectrum_inversion_tolerance, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

  if (opt.monitor)
  {
    PetscViewerAndFormat *vf;
    ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT, &vf); TiberPetscUtils::checkerr(ierr);
    ierr = KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorResidual, vf, 0);
  }

  return ierr;

}

//--------------------------------------------------------------//
int EigenSolver::number_of_converged_eigenvalues()
{
  int ierr, nconv;
  ierr =  EPSGetConverged(eps,&nconv);TiberPetscUtils::checkerr(ierr);



  return(nconv);
}


//--------------------------------------------------------------//
double EigenSolver::get_eigenvalue(int i)
{
  int ierr;
  PetscScalar ev, ev_i;


  // in case of spectrum folding, we recover the real eigenvalue
  // by calculating the Rayleigh quotient
  if (fold_ctx != nullptr)
  {
    Vec eigen_vector;
    ierr = MatCreateVecs(A, PETSC_NULLPTR, &eigen_vector);
    TiberPetscUtils::checkerr(ierr);

    ierr = EPSGetEigenpair(eps, i, &ev, &ev_i, eigen_vector, PETSC_NULLPTR);
    TiberPetscUtils::checkerr(ierr);

    ierr = RayleighQuotient(A, eigen_vector, &ev);
    VecDestroy(&eigen_vector);
  }
  else
    ierr = EPSGetEigenvalue(eps, i, &ev,  &ev_i);

  TiberPetscUtils::checkerr(ierr);

  double eigen_value = PetscRealPart(ev);

  return(eigen_value);

}
//----------------------------------------------------------------//


void EigenSolver::get_eigen_vector(int i, std::vector<Complex>& eigen_vector_out)
{
  int ierr, vec_size;
  PetscScalar kr, ki;
  Vec eigen_vector;



  MatCreateVecs(A,PETSC_NULLPTR,&eigen_vector);

  EPSGetEigenpair(eps,i,&kr,&ki,eigen_vector,PETSC_NULLPTR);

  PetscScalar *loc_part;
  VecGetArray(eigen_vector, &loc_part);
  VecGetLocalSize(eigen_vector, &vec_size);

  eigen_vector_out.resize(vec_size);

  for (int j= 0; j < vec_size; j++)
  {
    eigen_vector_out[j] = loc_part[j];
  }

  VecRestoreArray(eigen_vector, &loc_part);


  ierr = VecDestroy(&eigen_vector);

}



//-----------------------------------------------------------------------------//
void EigenSolver::set_initial_vector( const std::vector<Complex>& in_vector)
{

  int ierr;

  Vec v0;



  MatCreateVecs(A,PETSC_NULLPTR,&v0);

  PetscScalar y[_size_of_matrix];
  PetscInt ix[_size_of_matrix];

  for (int j= 0; j < _size_of_matrix ; j++)
  {
    y[j]  = in_vector[j];

    ix[j] = j;
  }


  VecSetValues(v0,_size_of_matrix,ix,y,INSERT_VALUES);


  //EPSSetInitialVector(eps, v0);


  ierr = VecDestroy(&v0);
}


//-----------------------------------------------------------------------------//

int EigenSolver::prepare_slepc(MPI_Comm comm)
{
  /*
     Create eigensolver context
  */
  int ierr;

  slepc_comm = comm;

  // Reset the LU fallback flag: each new k-point starts fresh and the
  // iterative KSP gets a fair chance again (it may work fine at other shifts).
  lu_fallback_active = false;

//  if (eps == NULLPTR)
  {
    ierr = EPSCreate(slepc_comm,&eps);TiberPetscUtils::checkerr(ierr);
  }


  return(ierr);
}


void
EigenSolver::set_deflation_space(
    const std::vector<const std::vector<Complex>*>& solutions)
{
  //if (solutions.size() == 0)
    return;


  Vec* defl = new Vec[solutions.size()];
  //vector<PetscVector<Complex>*> vecs(solutions.size(), NULLPTR);

  PetscScalar y[_size_of_matrix];
  PetscInt ix[_size_of_matrix];

  for (int j= 0; j < _size_of_matrix ; j++)
    ix[j] = j;


  for (int i = 0; i < solutions.size(); ++i)
  {
    //vecs[i] = new PetscVector<Complex>(libMesh::CommWorld);
    //vecs[i]->init(solutions[i]->size());
    //*vecs[i] = *solutions[i];
    MatCreateVecs(A, PETSC_NULLPTR, &defl[i]);
    for (int j= 0; j < _size_of_matrix ; j++)
      y[j]  = (*solutions[i])[j];
    VecSetValues(defl[i], _size_of_matrix, ix, y, INSERT_VALUES);
  }

  PetscInt defl_dim = solutions.size();
  EPSSetDeflationSpace(eps, defl_dim, defl);
}

//-----------------------------------------------------------------------------//

int EigenSolver::clear_slepc()
{
  /*
    Free memory
   */
  //EPSRemoveDeflationSpace(eps);

  int ierr;
  ierr = MatDestroy(&A);TiberPetscUtils::checkerr(ierr);
  {
    PetscBool generalized;
    ierr = EPSIsGeneralized(eps,&generalized); TiberPetscUtils::checkerr(ierr);

    if ( generalized)  ierr = MatDestroy(&B);TiberPetscUtils::checkerr(ierr);
  }

  if (fold_ctx != nullptr)
  {
    VecDestroy(&fold_ctx->w);
    delete fold_ctx;
    fold_ctx = nullptr;
  }


  // NOTE: with real MPI this leads to too many communicators
  // in a future version of SLEPc one could maybe use EPSReset()
  ierr = EPSDestroy(&eps);TiberPetscUtils::checkerr(ierr);
  //eps = NULLPTR;


  for (int i = 0; i < _deflation_space.size(); ++i)
    VecDestroy(&_deflation_space[i]);

  _deflation_space.clear();


  return(ierr);
}

//-------------------------------------------------------------//
int EigenSolver::do_solve(const SLEPCoptions& opt)
{

  int ierr;

  int ncv, nconv;


  //print_options(opt);


  if (opt.ev_number > 8)
    ncv =  4*opt.ev_number;
  else
    ncv = 32;



  if (ncv > _size_of_matrix) ncv = _size_of_matrix;

#if ((SLEPC_VERSION_MAJOR == 2) && (SLEPC_VERSION_MINOR == 3) && \
    (SLEPC_VERSION_SUBMINOR <= 2))
  if (opt.monitor) EPSSetMonitor(eps, EPSDefaultMonitor, PETSC_NULLPTR);
#else
  if (opt.monitor)
  {
    ierr = EPSMonitorCancel(eps); TiberPetscUtils::checkerr(ierr);
    PetscViewerAndFormat *vf;
    ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT, &vf); TiberPetscUtils::checkerr(ierr);
    ierr = EPSMonitorSet(eps, (PetscErrorCode (*)(EPS, PetscInt, PetscInt, PetscScalar*, PetscScalar*, PetscReal*, PetscInt, void*))EPSMonitorAll,
                         vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); TiberPetscUtils::checkerr(ierr);
  }
#endif



#if (SLEPC_VERSION_MAJOR >= 3)
  //ierr = EPSSetDimensions(eps,opt.ev_number, PETSC_DECIDE, PETSC_DECIDE); TiberPetscUtils::checkerr(ierr);
  ierr = EPSSetDimensions(eps,opt.ev_number, ncv, PETSC_DECIDE); TiberPetscUtils::checkerr(ierr);
#else
  ierr = EPSSetDimensions(eps,opt.ev_number, ncv); TiberPetscUtils::checkerr(ierr);
#endif

  ST st;
  EPSGetST(eps, &st);
  KSP ksp;
  STGetKSP(st, &ksp);
  PC pc;
  KSPGetPC(ksp, &pc);
  MatSolverType ms_type;
  PCFactorGetMatSolverType(pc, &ms_type);
  ostringstream os;
  os << "Solver package used in Eigensolver: " << ms_type << endl;
  Messages::info(os.str());


  // ---------------------------------------------------------------------
  // Safety net for the shift-and-invert linear solves.
  // The default (ilu + bcgsl) is an iterative solve that is fine for large
  // sparse Hamiltonians but can stall (DIVERGED_ITS) on small, dense-like or
  // strongly indefinite (H - shift) matrices, e.g. coarse-grained ones. If the
  // eigensolver failed *because* the inner KSP diverged, retry once with a
  // direct LU factorization instead of aborting the whole simulation.
  //
  // PETSc's default error handler aborts the process when an error propagates
  // up through EPSSolve (e.g. DIVERGED_ITS from KSPSolve). To allow the first
  // EPSSolve to fail gracefully and return an error code instead of aborting,
  // we temporarily install PetscReturnErrorHandler which converts fatal errors
  // into non-zero return codes. The original handler is restored immediately
  // after, whether the call succeeded or failed.
  //
  // Restrictions, so that no working configuration changes behaviour:
  //  - only reached if EPSSolve already failed;
  //  - only for the shift-and-invert transform whose KSP is still iterative;
  //  - only on a single MPI process (PETSc's own LU is serial; with several
  //    processes the user must choose mumps/mkl_pardiso explicitly).
  // ---------------------------------------------------------------------
  PetscPushErrorHandler(PetscReturnErrorHandler, NULL);
  ierr = EPSSolve(eps);
  PetscPopErrorHandler();

  // Also treat insufficient convergence as a failure: if EPSSolve returned
  // success but converged fewer eigenvalues than requested, the iterative
  // inner solve likely produced inaccurate (H-shift)^{-1} actions, causing
  // the eigensolver to converge to wrong eigenpairs rather than diverge
  // outright.  Treat this the same way as DIVERGED_ITS.
  if (ierr == 0)
  {
    PetscInt nconv_check = 0;
    EPSGetConverged(eps, &nconv_check);
    if (nconv_check < static_cast<PetscInt>(opt.ev_number))
      ierr = PETSC_ERR_CONV_FAILED;
  }

  if (ierr != 0)
  {
    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    KSPGetConvergedReason(ksp, &reason);

    KSPType ksp_type_now = nullptr;
    KSPGetType(ksp, &ksp_type_now);

    STType st_type_now = nullptr;
    STGetType(st, &st_type_now);

    PetscMPIInt comm_size_now = 0;
    MPI_Comm_size(slepc_comm, &comm_size_now);

    const bool ksp_diverged = (reason < 0);
    const bool ksp_iterative = (ksp_type_now != nullptr) &&
                               (std::string(ksp_type_now) != std::string(KSPPREONLY));
    const bool is_sinvert = (st_type_now != nullptr) &&
                            (std::string(st_type_now) == std::string(STSINVERT));

    if (ksp_diverged && ksp_iterative && is_sinvert && comm_size_now == 1)
    {
      // info (not Messages::warning): a warning may require acknowledgment when
      // stop_on_warning is set, which would block a batch run on a recoverable event.
      Messages::info("  (slepc) WARNING: iterative shift-and-invert linear solve did not converge; "
                     "retrying with a direct LU factorization "
                     "(set pc_type = lu in the Solver block to skip this retry)");

      PetscBool generalized = PETSC_FALSE;
      EPSIsGeneralized(eps, &generalized);

      // In SLEPc >= 3.23, a failed Krylov-Schur iteration leaves multiple
      // internal objects (BV, DS and its dense matrices) in dirty states that
      // cannot be cleanly reset through the public API: EPSReset clears BV but
      // not DS, and DSReset triggers MatDestroy on a matrix that still has an
      // unreleased submatrix (matinuse != 0), causing a second fatal error.
      // The only safe way to get a clean EPS before the retry is to destroy
      // the existing context and create a fresh one. Before destroying, we
      // must force-release any DS matrices left checked out by the failed solve.
      // The MPI communicator concern noted in clear_slepc() ("too many
      // communicators with real MPI") does not apply here because this branch
      // is already guarded to single-process runs (comm_size_now == 1).
      force_release_ds_matrices(eps);
      EPSDestroy(&eps);
      EPSCreate(slepc_comm, &eps);

      // Re-supply operators to the new EPS.
      EPSSetOperators(eps, A, (generalized == PETSC_TRUE) ? B : PETSC_NULLPTR);
      EPSSetProblemType(eps, (generalized == PETSC_TRUE) ? EPS_GHEP : EPS_HEP);

      // Reapply the user-level EPS settings that were on the old context.
      EPSSetTolerances(eps, opt.eps_tolerance, opt.eps_max_it);
      EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
      EPSSetTarget(eps, opt.spectrum_shift);
      EPSSetType(eps, (opt.solver_type == "arnoldi") ? EPSARNOLDI : EPSKRYLOVSCHUR);

      {
        ST st_retry;
        EPSGetST(eps, &st_retry);
        STSetType(st_retry, STSINVERT);

        KSP ksp_retry;
        STGetKSP(st_retry, &ksp_retry);
        PC pc_retry;
        KSPGetPC(ksp_retry, &pc_retry);

        // Switch the inner linear solver to a direct LU factorization.
        KSPSetType(ksp_retry, KSPPREONLY);
        PCSetType(pc_retry, PCLU);
        PCFactorSetMatSolverType(pc_retry, opt.solver_package.c_str());
      }

      // Restore the subspace size that was set in do_solve().
      {
        PetscInt ncv_retry = (opt.ev_number > 8) ? 4 * opt.ev_number : 32;
        if (ncv_retry > _size_of_matrix) ncv_retry = _size_of_matrix;
        EPSSetDimensions(eps, opt.ev_number, ncv_retry, PETSC_DECIDE);
      }

      // Re-apply any deflation space accumulated from previous rounds.
      if (opt.use_deflation_space && !_deflation_space.empty())
      {
        PetscInt defl_dim = _deflation_space.size();
        EPSSetDeflationSpace(eps, defl_dim, _deflation_space.data());
      }

      ierr = EPSSolve(eps);

      // Retry succeeded: remember that iterative KSP is unreliable for this
      // problem so that subsequent rounds in the same loop use LU directly.
      if (ierr == 0)
        lu_fallback_active = true;
    }
  }
  if (ierr == 0 && opt.use_deflation_space)
  {
    ierr = EPSGetConverged(eps, &nconv);  TiberPetscUtils::checkerr(ierr);

    // Re-query ncv from the EPS context: after a retry the EPS object may be
    // a freshly created one with a different ncv than the variable computed
    // at the top of this function.
    PetscInt ncv_actual, mpd_actual, nev_actual;
    EPSGetDimensions(eps, &nev_actual, &ncv_actual, &mpd_actual);

    // EPSGetInvariantSubspace fills the whole ncv-dimensional invariant
    // subspace, not only the nconv converged eigenvectors, so the buffer has to
    // be sized with ncv.
    Vec* v = new Vec[ncv_actual];
    for (int i = 0; i < ncv_actual; ++i)
      MatCreateVecs(A,PETSC_NULLPTR,&v[i]);

    EPSGetInvariantSubspace(eps, v);

    // _deflation_space owns the vectors, they are released by clear_slepc(),
    // so the ones handed to the deflation space must be copies.
    for (int i = 0; i < nconv; ++i)
    {
      Vec d;
      ierr = VecDuplicate(v[i], &d); TiberPetscUtils::checkerr(ierr);
      _deflation_space.push_back(d);
    }

    PetscInt defl_dim = _deflation_space.size();
    ierr = EPSSetDeflationSpace(eps, defl_dim, _deflation_space.data());
    TiberPetscUtils::checkerr(ierr);

    for (int i = 0; i < ncv_actual; ++i)
      VecDestroy(&v[i]);
    delete [] v;
  }

  return ierr;

}


//--------------------------------------------------------------//
int EigenSolver::init_H_matrix(unsigned int n)
{

  int ierr;

  ierr = MatCreate(slepc_comm, &A);
  TiberPetscUtils::checkerr(ierr);
  //ierr = MatSetType(A, MATAIJ);
  //TiberPetscUtils::checkerr(ierr);



  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);

  TiberPetscUtils::checkerr(ierr);

  return(ierr);

}


//----------------------------------------------------------//
int EigenSolver::init_S_matrix(unsigned int n)
{

  int ierr;

  ierr = MatCreate(slepc_comm,&B);
  TiberPetscUtils::checkerr(ierr);
  //ierr = MatSetType(B, MATAIJ);
  //TiberPetscUtils::checkerr(ierr);


  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n);
  TiberPetscUtils::checkerr(ierr);


  return(ierr);

}



void EigenSolver::finalize_matrix_assembly(const char matrix)
{
  int ierr;

  if (matrix == 'H')
  {
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    TiberPetscUtils::checkerr(ierr);

    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    TiberPetscUtils::checkerr(ierr);
  }
  else if (matrix == 'S')
  {
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    TiberPetscUtils::checkerr(ierr);

    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
    TiberPetscUtils::checkerr(ierr);
  }

}


void EigenSolver::insert_matrix_row(const char matrix, int row,
    const std::vector<unsigned int>& colums,
    const std::vector<Complex>& value_vector, int indexing_base)
{
  int ierr;
  int number_of_columns =  colums.size();
  PetscInt col[number_of_columns];
  PetscScalar value[number_of_columns];

  for (unsigned int i = 0; i < number_of_columns; i++)
  {
    col[i] = colums[i] - indexing_base;
    value[i] = value_vector[i];
  }


  if (matrix == 'H')
    ierr = MatSetValues(A,1,&row,number_of_columns,col,value,INSERT_VALUES);
  else if (matrix == 'S')
    ierr = MatSetValues(B,1,&row,number_of_columns,col,value,INSERT_VALUES);

  TiberPetscUtils::checkerr(ierr);
}

//------------------------------------------------------------------------------------//
int EigenSolver::preallocate_matrix(const char matrix, unsigned int N, unsigned int n, vector<int>& d_nnz, vector<int>& o_nnz)
{

  int ierr;
  
  if (matrix == 'H')
  {
    ierr = MatCreateAIJ(slepc_comm, n, n, N, N,
                        d_nnz.size(), d_nnz.data(),
                        o_nnz.size(), o_nnz.data(), &A);
  }
  else if (matrix == 'S')
  {
    ierr = MatCreateAIJ(slepc_comm, n, n, N, N,
                        d_nnz.size(), d_nnz.data(),
                        o_nnz.size(), o_nnz.data(), &B);
  }

  _size_of_matrix = N;

  return(ierr);
}

//------------------------------------------------------------------------------------//
double  EigenSolver::get_shift(void)
{
  //PetscErrorCode ierr;
  //ST st;
  //PetscScalar shift;             //

  //ierr = EPSGetST(eps, &st);     // ierr = EPSGetST(eps, st);

  //TiberPetscUtils::checkerr(ierr);

  //ierr = STGetShift(st, &shift); // ierr = STGetShift(*st, &shift);

  //TiberPetscUtils::checkerr(ierr);

  //return(real(shift));
  return(shift);

}
//------------------------------------------------------------------------------------//
bool EigenSolver::check_matrices(double tol, bool verbose)
{
  PetscErrorCode ierr;
  PetscBool is;
  bool ans;

  ans = false;
  if(A->hermitian)
  {
    if (verbose)
      std::cout<<"   (ES) Hamiltonian is defined Hermitian"<<std::endl;
    ans = true;
  }
  else
  {
    ierr = MatIsHermitian(A, tol, &is);

    if (is==PETSC_TRUE)
    {
      if (verbose)
        std::cout<<"   (ES) Hamiltonian is Hermitian within "<<tol<<std::endl;
    }
    else
    {
      if (verbose)
        std::cout<<"   (ES) Hamiltonian is NOT Hermitian!"<<std::endl;
    }
    ans = (is==PETSC_TRUE);
  }



  if(B->hermitian)
  {
    if (verbose)
      std::cout<<"   (ES) Overlap is defined Hermitian"<<std::endl;
    is = PETSC_TRUE;
  }
  else
  {
    ierr = MatIsHermitian(B, tol, &is);

    if (is==PETSC_TRUE)
    {
      if (verbose)
        std::cout<<"   (ES) Overlap is Hermitian within "<<tol<<std::endl;
    }
    else
    {
      if (verbose)
        std::cout<<"   (ES) Overlap is NOT Hermitian!"<<std::endl;
    }
  }

  return (is==PETSC_TRUE) && ans;

}

bool EigenSolver::check_matrices(void)
{
  PetscErrorCode ierr;
  PetscBool is;
  bool ans;

  ans = false;

  ierr = EPSIsHermitian(eps, &is);

  ans = is==PETSC_TRUE;

  //ierr = SlepcIsHermitian(B, &is);

  //return (is==PETSC_TRUE) && ans;
  return(ans);

}



void set_sub_pc(PC pc, PCType pc_type)
{

  PetscErrorCode ierr;

  // To store array of local KSP contexts on this processor
  KSP* subksps;

  // the number of blocks on this processor
  PetscInt n_local;

  // Fill array of local KSP contexts
  ierr = PCBJacobiGetSubKSP(pc, &n_local, PETSC_NULLPTR, &subksps);

  // Loop over sub-ksp objects, set ILU preconditioner
  for (PetscInt i = 0; i < n_local; ++i)
  {
    // Get pointer to sub KSP object's PC
    PC subpc;
    ierr = KSPGetPC(subksps[i], &subpc);

    // Set requested type on the sub PC
    ierr = PCSetType(subpc, pc_type);
  }
}
