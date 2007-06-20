// $Id$

#include "PardisoLinearSolver.h"


#include "petsc_vector.h"
#include "petsc_matrix.h"

#include "petscmat.h"

#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscversion.h>
# include <petscksp.h>
}
#else
# include <petscversion.h>
# include <petscksp.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif

extern "C" {int PARDISO
	(void *, int *, int *, int *, int *, int *,
	double *, int *, int *, int *, int *, int *,
	int *, double *, double *, int *);
}


extern "C"{ int omp_get_max_threads();}



void
PardisoLinearSolver::clear(void)
{

}


void
PardisoLinearSolver::init(void)
{

    mtype = 11; /* Real unsymmetric matrix */

    nrhs = 1;   /* Number of right hand sides. */     

    maxfct = 1;  /* Maximum number of numerical factorizations. */

    mnum =  1;  /* Which factorization to use. */

    error = 0;  /* Initialize error flag */

    msglvl = 1; /* Print statistical information in file */


/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */



	/* Pardiso control parameters. */

	for (int i = 0; i < 64; i++) {
		iparm[i] = 0;
	}
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
	/* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = omp_get_max_threads();
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 0; /* Not in use */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Not in use */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */


}


std::pair<unsigned int, Real> 
PardisoLinearSolver::solve(SparseMatrix<Number>&  matrix_in,
			     SparseMatrix<Number>&  precond_in,
			     NumericVector<Number>& solution_in,
			     NumericVector<Number>& rhs_in,
			     const double tol,
			     const unsigned int m_its)
{


  PetscMatrix<Number>* matrix   = dynamic_cast<PetscMatrix<Number>*>(&matrix_in); 
  PetscMatrix<Number>* precond  = dynamic_cast<PetscMatrix<Number>*>(&precond_in);
  PetscVector<Number>* solution = dynamic_cast<PetscVector<Number>*>(&solution_in);
  PetscVector<Number>* rhs      = dynamic_cast<PetscVector<Number>*>(&rhs_in);


   // We cast to pointers so we can be sure that they succeeded
  // by comparing the result against NULL.

  assert(matrix   != NULL);
  assert(precond  != NULL);
  assert(solution != NULL);
  assert(rhs      != NULL);



  PetscErrorCode ierr;
  PetscInt       *ia, *ja;
  PetscScalar    *a, *b, *x;
  PetscInt       n;
  PetscTruth     done;
  Mat            C;
  Vec            u,sol;            



  //Initialize Petsc-------------------------------
  int argc = 0;
  char** args[0];
  static char help[] = "Test";
  PetscInitialize(&argc,args,(char *)0,help);
  //------------------------------------------------ 


  //Get matrix
  C = matrix->mat();
  ierr =   MatGetArray(C, &a);
  _checkerr(ierr);
  

  //Get rsh
  u = rhs->vec();
  ierr = VecGetArray(u, &b);
  _checkerr(ierr);

  //Get sol
  sol = solution->vec();
  ierr = VecGetArray(sol, &x);
  _checkerr(ierr);

  //Get n ia and ja
  ierr = MatGetRowIJ(C, 1, PETSC_FALSE, &n, &ia, &ja, &done);
  _checkerr(ierr);
  

  if (done)
  {
    
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    
    
    /* Auxiliary variables. */
    double ddum; /* Double dummy */
    int idum; /* Integer dummy. */
    int  phase;
    
    
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    for (int i = 0; i < 64; i++) {
      pt[i] = 0;
    }
    


    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
  
    phase = 11;

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, a, ia, ja, &idum, &nrhs,
	     iparm, &msglvl, &ddum, &ddum, &error);



    if (error != 0) {
      printf("\nERROR during symbolic factorization: %d\n", error);
      exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors = %d", iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, a, ia, ja, &idum, &nrhs,
	     iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
      printf("\nERROR during numerical factorization: %d", error);
      exit(2);
    }
    printf("\nFactorization completed ... ");
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2; /* Max numbers of iterative refinement steps. */
    /* Set right hand side to one. */
    for (int i = 0; i < n; i++) {
      b[i] = 1;
    }
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, a, ia, ja, &idum, &nrhs,
	     iparm, &msglvl, b, x, &error);
    if (error != 0) {
      printf("\nERROR during solution: %d", error);
      exit(3);
    }
    printf("\nSolve completed ... ");
    printf("\nThe solution of the system is: ");
    for (int i = 0; i < n; i++) {
      printf("\n x [%d] = % f", i, x[i] );
    }
    printf ("\n");
    
    
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, &ddum, ia, ja, &idum, &nrhs,
	     iparm, &msglvl, &ddum, &ddum, &error);
    
    
  }

  
  ierr = VecRestoreArray(sol, &x);//CHKERRQ(ierr); 
  
    
  //Insert the vec object in the solution
  PetscVector<Number> sol_temp(sol);
  *solution = sol_temp;


}



