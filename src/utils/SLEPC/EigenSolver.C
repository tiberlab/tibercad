
static char help[] = "Solves a generalized eigensystem Ax=kBx with matrices loaded from a file.\n\n"
  "This example works for both real and complex numbers.\n\n"
  "The command line options are:\n\n"
  "  -f1 <filename>, where <filename> = matrix (A) file in PETSc binary form.\n"
  "  -f2 <filename>, where <filename> = matrix (B) file in PETSc binary form.\n\n";

#include <iostream>

#include <string.h>
#include "slepceps.h"
#include "EigenSolver.h"

namespace
{
Mat A;
Mat B;
EPS eps;
}

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
  PetscMPIInt    rank,size;
  ST st;
  KSP ksp;
  


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 
 

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
 

   
  //
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,opt.H_file_name.c_str(),PETSC_FILE_RDONLY,&viewer);CHKERRQ(ierr); //their
  ierr = MatLoad(viewer,MATAIJ,&A);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);


  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,opt.S_file_name.c_str(),PETSC_FILE_RDONLY,&viewer);CHKERRQ(ierr); //their
  ierr = MatLoad(viewer,MATAIJ,&B);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);


  

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
     Set operators. In this case, it is a generalized eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  // ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);
  
 

 
  ierr = EPSSetTolerances(eps,opt.eps_tolerance,opt.eps_max_it);  CHKERRQ(ierr);


  if (opt.solver_type == "arnoldi")
  {
    ierr = EPSSetProblemType(eps,EPS_GNHEP);CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSARNOLDI); CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE);CHKERRQ(ierr);
 

    ierr = EPSGetST(eps,&st); CHKERRQ(ierr);
    ierr = STSetShift(st, opt.spectrum_shift);CHKERRQ(ierr); 
    
   
    ierr = STSetType(st,STSINV); CHKERRQ(ierr);
    ierr = STGetKSP(st, &ksp);CHKERRQ(ierr);

    ierr = KSPSetType( ksp, KSPBCGS);CHKERRQ(ierr);


    //rtol, abstol, dtol, maxits
    ierr = KSPSetTolerances(ksp,1e-10, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
   
  }
  else if (opt.solver_type == "lapack")
  {
    ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSLAPACK);
    ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);
  }
  else if (opt.solver_type == "arpack")
  {
    ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSARPACK);
    ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);

    if (std::abs(opt.spectrum_shift) >1e-8)
    {
      ierr = EPSGetST(eps,&st); CHKERRQ(ierr);
      ierr = STSetShift(st, opt.spectrum_shift);CHKERRQ(ierr); 
    
   
      ierr = STSetType(st,STSHIFT); CHKERRQ(ierr);

      ierr = STGetKSP(st, &ksp);CHKERRQ(ierr);

      //      ierr = KSPSetType( ksp, KSPBCGS);CHKERRQ(ierr);

      //rtol, abstol, dtol, maxits
      ierr = KSPSetTolerances(ksp,1e-10, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);

    }
  }
 





  ierr = do_solve(opt);

 

  return ierr;
}

//------------------------------------------------------------------------------//



int EigenSolver::eig_value_problem(const EigenSolver::SLEPCoptions& opt )
{
  
  EPSType     type;
  PetscReal   error, tol, re, im;
  PetscScalar kr, ki;
  int         nev, ierr, maxit, i, its, lits, nconv;
  char        filename[256];
  PetscViewer viewer, viewer_out, viewer_eigvals;
  PetscTruth  flg;
  PetscMPIInt    rank,size;
  ST st;
  KSP ksp;
  


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 
 

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
 

   
  //
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,opt.H_file_name.c_str(),PETSC_FILE_RDONLY,&viewer);CHKERRQ(ierr); //their
  ierr = MatLoad(viewer,MATAIJ,&A);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);



  

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
     Set operators. In this case, it is a generalized eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,PETSC_NULL);CHKERRQ(ierr);
 
  
 

 
  ierr = EPSSetTolerances(eps,opt.eps_tolerance,opt.eps_max_it);  CHKERRQ(ierr);


  if (opt.solver_type == "arnoldi")
  {
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSARNOLDI); CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE);CHKERRQ(ierr);
 

    ierr = EPSGetST(eps,&st); CHKERRQ(ierr);
    ierr = STSetShift(st, opt.spectrum_shift);CHKERRQ(ierr); 
    
   
    ierr = STSetType(st,STSINV); CHKERRQ(ierr);
    ierr = STGetKSP(st, &ksp);CHKERRQ(ierr);

    ierr = KSPSetType( ksp, KSPBCGS);CHKERRQ(ierr);


    //rtol, abstol, dtol, maxits
    ierr = KSPSetTolerances(ksp,1e-10, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
   
  }
  else if (opt.solver_type == "lapack")
  {
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSLAPACK);
    ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);
  }
  else if (opt.solver_type == "arpack")
  {
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSARPACK);
    ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);

    if (std::abs(opt.spectrum_shift) >1e-8)
    {
      ierr = EPSGetST(eps,&st); CHKERRQ(ierr);
      ierr = STSetShift(st, opt.spectrum_shift);CHKERRQ(ierr); 
    
   
      ierr = STSetType(st,STSHIFT); CHKERRQ(ierr);

      ierr = STGetKSP(st, &ksp);CHKERRQ(ierr);

      //      ierr = KSPSetType( ksp, KSPBCGS);CHKERRQ(ierr);

      //rtol, abstol, dtol, maxits
      ierr = KSPSetTolerances(ksp,1e-10, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);

    }

  }
 
  ierr = do_solve(opt);

 

  return ierr; 

 
}

//--------------------------------------------------------------//
int EigenSolver::number_of_converged_eigenvalues()
{
  int ierr, nconv; 
  ierr =  EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  return(nconv);
}
//--------------------------------------------------------------//
double EigenSolver::get_eigenvalue( int i)
{
  int ierr;
  PetscScalar ev, ev_i;
 
  ierr = EPSGetValue(eps, i, &ev,  &ev_i);

  double eigen_value = PetscRealPart(ev);

  return(eigen_value);
 
}
//----------------------------------------------------------------//


 void EigenSolver::get_eigen_vector( int i, std::vector<Complex>& eigen_vector_out)
{
  int ierr, vec_size;
  PetscScalar kr, ki;
  Vec eigen_vector;

 

  MatGetVecs(A,PETSC_NULL,&eigen_vector);

  EPSGetEigenpair(eps,i,&kr,&ki,eigen_vector,PETSC_NULL);

 
  VecGetSize(eigen_vector, &vec_size);
  
  eigen_vector_out.resize(vec_size);

  {

    PetscInt ix[vec_size];

    for (int j= 0; j < vec_size; j++) ix[j] = j;

    PetscScalar y[vec_size];

    VecGetValues(eigen_vector,vec_size,ix,y);

    for (int j= 0; j < vec_size; j++)
    {
      eigen_vector_out[j] = Complex(  PetscRealPart(y[j]), PetscImaginaryPart(y[j]) );  
    }

  }
}
//-----------------------------------------------------------------------------//

int EigenSolver::prepare_slepc()
{
  /* 
     Create eigensolver context
  */
  int ierr;
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  return(ierr);
}

//-----------------------------------------------------------------------------//

int EigenSolver::clear_slepc()
{
  /*
    Free memory
   */
  int ierr;
  ierr = MatDestroy(A);CHKERRQ(ierr);
  {
    PetscTruth generalized;
    ierr = EPSIsGeneralized(eps,&generalized); CHKERRQ(ierr);

    if ( generalized)  ierr = MatDestroy(B);CHKERRQ(ierr);
  }

  ierr = EPSDestroy(eps);CHKERRQ(ierr);
 

  

  return(ierr);
}

//-------------------------------------------------------------//
int EigenSolver::do_solve(const SLEPCoptions& opt)
{

  int ierr;

  int ncv;
  

  if (opt.ev_number > 8)
    ncv =  4*opt.ev_number;
  else
    ncv = 32;
   
  // ierr = EPSSetDimensions(eps,opt.ev_number, PETSC_DECIDE); CHKERRQ(ierr);

  ierr = EPSSetDimensions(eps,opt.ev_number, ncv); CHKERRQ(ierr);


  ierr = EPSSolve(eps);CHKERRQ(ierr);
  
 

  return ierr;

}


//--------------------------------------------------------------//
int EigenSolver::init_H_matrix(unsigned int n)
{

  int ierr;

  ierr = MatCreate(PETSC_COMM_WORLD,&A);
  CHKERRQ(ierr);
  
  if (ierr = 0)
  {
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
    CHKERRQ(ierr);
  
    if (ierr = 0)
    {
      ierr = MatSetFromOptions(A);
      CHKERRQ(ierr);
    }
  }

  return(ierr);
  
}


//----------------------------------------------------------//
int EigenSolver::init_S_matrix(unsigned int n)
{

  int ierr;

  ierr = MatCreate(PETSC_COMM_WORLD,&B);
  CHKERRQ(ierr);
  
  if (ierr = 0)
  {
    ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n);
    CHKERRQ(ierr);
  
    if (ierr = 0)
    {
      ierr = MatSetFromOptions(B);
      CHKERRQ(ierr);
    }
  }

  return(ierr);
  
}
//-----------------------------------------------------------------//
int EigenSolver::finalize_H_assembly(void)
{
  int ierr;

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  if (ierr = 0)
  {
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
  }

  return(ierr);
}
//------------------------------------------------------------------//
int EigenSolver::finalize_S_assembly(void)
{
  int ierr;

  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  if (ierr = 0)
  {
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
  }

  return(ierr);
}
//-------------------------------------------------------------------//

int EigenSolver::insert_H_row( int row, const std::vector<unsigned int>& colums, const std::vector<Complex>& value_vector)
{
  int ierr;
  int number_of_columns =  colums.size();
  PetscInt col[number_of_columns];
  PetscScalar value[number_of_columns];
  
  for (unsigned int i = 0; i < number_of_columns; i++)
  {
    col[i] = colums[i];
    value[i] = value_vector[i];
  }


  ierr = MatSetValues(A,1,&row,number_of_columns,col,value,INSERT_VALUES);

  CHKERRQ(ierr);
  return(ierr);
}

//-------------------------------------------------------------------//

int EigenSolver::insert_S_row( int row, const std::vector<unsigned int>& colums, const std::vector<Complex>& value_vector)
{
  int ierr;
  int number_of_columns =  colums.size();
  PetscInt col[number_of_columns];
  PetscScalar value[number_of_columns];
  
  for (unsigned int i = 0; i < number_of_columns; i++)
  {
    col[i] = colums[i];
    value[i] = value_vector[i];
  }


  ierr = MatSetValues(B,1,&row,number_of_columns,col,value,INSERT_VALUES);

  CHKERRQ(ierr);
  return(ierr);
}
