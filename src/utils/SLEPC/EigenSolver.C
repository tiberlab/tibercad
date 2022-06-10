// $Id$

#include <iostream>
#include <cassert>
#include <string>

#include "EigenSolver.h"
#include "RuntimeException.h"
#include "TiberPetscUtils.h"
#include "slepceps.h"
//#include "petscsys.h"

//#include <eigen_solver.h>

#include "petsc/private/matimpl.h"

//#include "libmesh/petsc_vector.h"

using namespace std;

namespace
{
  Mat A; // Hamiltonian
  Mat B; // S-matrix
  Mat M; // Shell matrix for spectrum folding
  Mat P; // preconditioner matrix
  EPS eps; // EigenSolver
  MPI_Comm slepc_comm;
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


int EigenSolver::_size_of_matrix;
//-------------------------------------------------------------//
void EigenSolver::slepc_init(int argc1, char** argv1, MPI_Comm comm)
{

  slepc_comm = comm;

  //Seems to work ^^
  //  TODO Looks poor, but in current version of petsc there is no methods for it. (In later releases there is ...)
  // It fix problem when during factorization MUMPS does not have have enouph memory.
/*
  int __empty_argc = 3;
  char** __empty_argv = new char*[3];
  __empty_argv[0] = "tibercad";
  __empty_argv[1] = "-mat_mumps_icntl_14";
  __empty_argv[2] = "100";

  SlepcInitialize(&__empty_argc,&__empty_argv,NULL,NULL);
*/
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
      MatCreateVecs(A, &fold_ctx->w, NULL);
      //VecDuplicate(x, &ctx->w);

      MatCreateShell(slepc_comm, nloc, mloc, N, N, fold_ctx, &M);
      MatShellSetOperation(M, MATOP_MULT, (void(*)(void))MatMult_Fold);
    }

    fold_ctx->target = shift;

  }

  // Set operators.
  if (evp_type == GENERALIZED)
  {
    ierr = EPSSetOperators(eps,A,B);TiberPetscUtils::checkerr(ierr);
    ierr = EPSSetProblemType(eps,EPS_GHEP);TiberPetscUtils::checkerr(ierr);
  }
  else
  {
    if (opt.spectral_trans == "folding")
      ierr = EPSSetOperators(eps, M, PETSC_NULL);
    else
      ierr = EPSSetOperators(eps, A, PETSC_NULL);

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


  if (opt.solver_type == "arnoldi" || opt.solver_type == "krylovshur")
  {
    //ierr = EPSSetProblemType(eps,EPS_GNHEP);TiberPetscUtils::checkerr(ierr);

    if (opt.solver_type == "arnoldi")
      ierr = EPSSetType(eps, EPSARNOLDI);
    else if (opt.solver_type == "krylovshur")
      ierr = EPSSetType(eps, EPSKRYLOVSCHUR);


    ierr = EPSGetST(eps,&st); TiberPetscUtils::checkerr(ierr);
    
   
    if (opt.spectral_trans == "folding")
    {
      //ierr = STSetType(st,STFOLD); TiberPetscUtils::checkerr(ierr);

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

    //set_ksp_and_pc(st, opt);
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
      //set_ksp_and_pc(st, opt);
      STGetKSP(st, &ksp);
      KSPSetType( ksp, KSPMINRES);
      KSPGetPC(ksp, &pc);
      PCSetType(pc, PCJACOBI);
      KSPSetTolerances(ksp,opt.spectrum_inversion_tolerance, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

      if (opt.monitor)
      {
        PetscViewerAndFormat *vf;
        ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT, &vf);
        TiberPetscUtils::checkerr(ierr);
        ierr = KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorResidual, vf, 0);
        TiberPetscUtils::checkerr(ierr);
      }
    }
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




int set_ksp_and_pc(ST st, const EigenSolver::SLEPCoptions& opt)
{
  int ierr;

  KSP ksp;
  ierr = STGetKSP(st, &ksp);

  PC pc;
  ierr = KSPGetPC(ksp,&pc);

#if ((SLEPC_VERSION_MAJOR < 3) || \
    ((SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR < 5)))
  PCSetOperators(pc, A, A, SAME_NONZERO_PATTERN);
#else
  PCSetOperators(pc, A, A);
#endif

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
      ierr =  PCSetType(pc, PCILU);

  }
  else if (opt.pc_type == "lu" )
  {
    // in principle, this should have worked with this if, but for some reason even
    // when running with a single process it complains about mpiaij matrix type
    if (comm_size > 1)
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
    ierr = MatCreateVecs(A, PETSC_NULL, &eigen_vector);
    TiberPetscUtils::checkerr(ierr);

    ierr = EPSGetEigenpair(eps, i, &ev, &ev_i, eigen_vector, PETSC_NULL);
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


void EigenSolver::get_eigen_vector( int i, std::vector<Complex>& eigen_vector_out)
{
  int ierr, vec_size;
  PetscScalar kr, ki;
  Vec eigen_vector;



  MatCreateVecs(A,PETSC_NULL,&eigen_vector);

  EPSGetEigenpair(eps,i,&kr,&ki,eigen_vector,PETSC_NULL);

  //VecGetSize(eigen_vector, &vec_size);
  

  PetscScalar *loc_part;
  VecGetArray(eigen_vector, &loc_part);
  VecGetLocalSize(eigen_vector, &vec_size);

  eigen_vector_out.resize(vec_size);

  for (int j= 0; j < vec_size; j++)
  {
    eigen_vector_out[j] = loc_part[j];
  }

  VecRestoreArray(eigen_vector, &loc_part);

  //PetscMPIInt rank;
  //MPI_Comm_rank(slepc_comm,&rank);
/*
  {

    PetscInt ix[vec_size];

    for (int j= 0; j < vec_size; j++) ix[j] = j;

    PetscScalar y[vec_size];

    VecGetValues(eigen_vector,vec_size,ix,y);

    for (int j= 0; j < vec_size; j++)
    {
      eigen_vector_out[j] = Complex(  PetscRealPart(y[j]), PetscImaginaryPart(y[j]) );
      //if (rank == 0) cerr << eigen_vector_out[j] << endl;
    }

  }
*/

  ierr = VecDestroy(&eigen_vector);

}
//-----------------------------------------------------------------------------//
void EigenSolver::set_initial_vector( const std::vector<Complex>& in_vector)
{

  int ierr;

  Vec v0;



  MatCreateVecs(A,PETSC_NULL,&v0);

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


//  if (eps == NULL)
  {
    ierr = EPSCreate(slepc_comm,&eps);TiberPetscUtils::checkerr(ierr);
//    ierr = EPSCreate(comm ,&eps);TiberPetscUtils::checkerr(ierr);
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
  //vector<PetscVector<Complex>*> vecs(solutions.size(), NULL);

  PetscScalar y[_size_of_matrix];
  PetscInt ix[_size_of_matrix];

  for (int j= 0; j < _size_of_matrix ; j++)
    ix[j] = j;


  for (int i = 0; i < solutions.size(); ++i)
  {
    //vecs[i] = new PetscVector<Complex>(libMesh::CommWorld);
    //vecs[i]->init(solutions[i]->size());
    //*vecs[i] = *solutions[i];
    MatCreateVecs(A, PETSC_NULL, &defl[i]);
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
  //eps = NULL;


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
  if (opt.monitor) EPSSetMonitor(eps, EPSDefaultMonitor, PETSC_NULL);
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


  ierr = EPSSolve(eps);

  if (ierr == 0 && opt.use_deflation_space)
  {
    ierr = EPSGetConverged(eps, &nconv);  TiberPetscUtils::checkerr(ierr);

    Vec* v = new Vec[nconv];
    for (int i = 0; i < nconv; ++i)
      MatCreateVecs(A,PETSC_NULL,&v[i]);

    EPSGetInvariantSubspace(eps, v);

    for (int i = 0; i < nconv; ++i)
      _deflation_space.push_back(v[i]);

    PetscInt defl_dim = _deflation_space.size();
    EPSSetDeflationSpace(eps, defl_dim, _deflation_space.data());
    //for (int i = 0; i < nconv; ++i)
    //  VecDestroy(&v[i]);
    //delete [] v;
    //VecDestroyVecs(v, nconv);
  }

  return ierr;

}


//--------------------------------------------------------------//
int EigenSolver::init_H_matrix(unsigned int n)
{

  int ierr;

  ierr = MatCreate(slepc_comm,&A);
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
  if(A->hermitian_set)
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



  if(B->hermitian_set)
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
  ierr = PCBJacobiGetSubKSP(pc, &n_local, PETSC_NULL, &subksps);

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
