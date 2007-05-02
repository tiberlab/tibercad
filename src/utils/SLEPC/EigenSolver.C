
static char help[] = "Solves a generalized eigensystem Ax=kBx with matrices loaded from a file.\n\n"
  "This example works for both real and complex numbers.\n\n"
  "The command line options are:\n\n"
  "  -f1 <filename>, where <filename> = matrix (A) file in PETSc binary form.\n"
  "  -f2 <filename>, where <filename> = matrix (B) file in PETSc binary form.\n\n";

#include <iostream>

#include <string.h>

#include "EigenSolver.h"
Mat EigenSolver::A;
Mat EigenSolver::B;
EPS EigenSolver::eps;

//-------------------------------------------------------------//
void EigenSolver::slepc_init()
{
  int argc1 = 0;
  char **argv1;
  argv1 = (char**) malloc(sizeof(char*));
  argv1[0] = (char*) malloc(0);
  

  SlepcInitialize(&argc1,&argv1,NULL,NULL);


  free(argv1[0]);
  free(argv1);
}

//--------------------------------------------------------------//
void  EigenSolver::slepc_done()
{
 SlepcFinalize();
}
//--------------------------------------------------------------//
int EigenSolver::eig_value_problem_general(const EigenSolver::SLEPCoptions& opt )
{
  
  EPSType     type;
  PetscReal   error, tol, re, im;
  PetscScalar kr, ki;
  int         nev, ierr, maxit, i, its, lits, nconv;
  char        filename[256];
  PetscViewer viewer, viewer_out, viewer_eigvals;
  PetscTruth  flg;
  Vec         eigen_vector;
  Vec         eig_vals;
  PetscMPIInt    rank,size;
  ST st;
  KSP ksp;
  
  //SlepcInitialize(NULL,NULL,(char*)0,help);
/*
  int argc1 = 0;
  char **argv1;
  
  argv1 = (char**) malloc(sizeof(char*));
  argv1[0] = (char*) malloc(0);
  //memset (argv1[0],'\0',256);
  //argv1[255] = '\0';
  
  //SlepcInitialize(&argc1,&argv1,(char*)0,help);

   SlepcInitialize(&argc1,&argv1,NULL,NULL);
 */
 

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 
 

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  //rr = PetscPrintf(PETSC_COMM_WORLD,"Number of processors = %d, rank = %d\n",size,rank);CHKERRQ(ierr);
  


#if defined(PETSC_USE_COMPLEX)
  // ierr = PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrices from binary files...\n");CHKERRQ(ierr);
#else
  // ierr = PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrices from binary files...\n");CHKERRQ(ierr);
#endif

   
  //
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,opt.H_file_name.c_str(),PETSC_FILE_RDONLY,&viewer);CHKERRQ(ierr); //their
  ierr = MatLoad(viewer,MATAIJ,&A);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);


  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,opt.S_file_name.c_str(),PETSC_FILE_RDONLY,&viewer);CHKERRQ(ierr); //their
  ierr = MatLoad(viewer,MATAIJ,&B);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  MatGetVecs(A,PETSC_NULL,&eigen_vector);

  

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /* 
     Set operators. In this case, it is a generalized eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  // ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);
  
  ierr = EPSSetProblemType(eps,EPS_GNHEP);CHKERRQ(ierr);

 
  ierr = EPSSetTolerances(eps,opt.eps_tolerance,opt.eps_max_it);  CHKERRQ(ierr);

  ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE);CHKERRQ(ierr);

  //ierr = EPSSetType(eps, EPSLAPACK);
  
  ierr = EPSSetType(eps, EPSARNOLDI); CHKERRQ(ierr);
	
  ierr = EPSGetST(eps,&st); CHKERRQ(ierr);

  ierr = STSetShift(st, opt.spectrum_shift);CHKERRQ(ierr); 

  ierr = STSetType(st,STSINV); CHKERRQ(ierr);


  ierr = STGetKSP(st, &ksp);CHKERRQ(ierr);

  ierr = KSPSetType( ksp, KSPBCGS);CHKERRQ(ierr);


  //rtol, abstol, dtol, maxits
  ierr = KSPSetTolerances(ksp,1e-10, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);

 

  ierr = EPSSetDimensions(eps,opt.ev_number, PETSC_DECIDE); CHKERRQ(ierr);
  //ierr = EPSSetDimensions(eps,opt.ev_number,opt.ev_number * 10 ); CHKERRQ(ierr); //test
  /*
     Set solver parameters at runtime
  */
  //ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetIterationNumber(eps, &its);CHKERRQ(ierr);
  //  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
  ierr = EPSGetNumberLinearIterations(eps, &lits);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %d\n",lits);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,PETSC_NULL);CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
     Get number of converged eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate eigenpairs: %d\n\n",nconv);CHKERRQ(ierr);

  if (nconv>0) {//there are converged solutions

     ierr = VecCreate(PETSC_COMM_WORLD,&eig_vals);CHKERRQ(ierr);
     ierr = VecSetSizes(eig_vals,PETSC_DECIDE,nconv);CHKERRQ(ierr);
     ierr = VecSetFromOptions(eig_vals);CHKERRQ(ierr);
   
    
    /*
       Display eigenvalues and relative errors
    */

    
    

     ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"eigvects_SLEPC.out",PETSC_FILE_CREATE,&viewer_out);CHKERRQ(ierr); //their
   
     ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"eigvals_SLEPC.out",PETSC_FILE_CREATE,&viewer_eigvals);CHKERRQ(ierr); //their
    
    /*
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k             ||Ax-kx||/||kx||\n"
         "  --------------------- ------------------\n" );CHKERRQ(ierr);

    */

    for( i=0; i<nconv; i++ ) {//eigen solutions loop
      /* 
         Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
         ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,eigen_vector,PETSC_NULL);CHKERRQ(ierr);
 
      /*
         Compute the relative error associated to each eigenpair
      */
      ierr = EPSComputeRelativeError(eps,i,&error);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX) 
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      

      if( im != 0.0 ) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," % 20.14f %+12f i \n",re,im);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"       % 6f      ",re); CHKERRQ(ierr);
      }
      // ierr = PetscPrintf(PETSC_COMM_WORLD," % 12f\n",error);CHKERRQ(ierr);

       

      ierr = VecAssemblyBegin(eigen_vector);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(eigen_vector);CHKERRQ(ierr);

      VecView(eigen_vector,viewer_out); //save eigen vector to disk //should be!!! 
      ierr =  VecSetValue( eig_vals, i , kr, INSERT_VALUES);  CHKERRQ(ierr);//write eivals into the vector eig_vals

     

    }
    ierr = VecAssemblyBegin(eig_vals);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(eig_vals);CHKERRQ(ierr);
   
    ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);
   
    VecView(eig_vals,viewer_eigvals); //save eigvals to disk  

    ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);

    ierr = PetscViewerDestroy(viewer_out);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Saved \n" );CHKERRQ(ierr);


    
  }
  
  /* 
     Free work space
  */

  ierr = EPSDestroy(eps);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = MatDestroy(B);CHKERRQ(ierr);
  ierr = VecDestroy(eigen_vector); CHKERRQ(ierr);
  ierr = VecDestroy(eig_vals); CHKERRQ(ierr);
  //  ierr = SlepcFinalize();CHKERRQ(ierr);
 
 

  return 0;
}

//------------------------------------------------------------------------------//
