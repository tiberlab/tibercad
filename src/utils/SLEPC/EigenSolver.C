// $Id$
 
#include <iostream>
#include <cassert>
#include <string>

#include "slepceps.h"
#include "EigenSolver.h"

#include "private/matimpl.h"

namespace
{
  Mat A; //Hamiltonian
  Mat B; //S-matrix
  EPS eps; //EigenSolver 
  double shift; //could be stored in ST but lapack does not apply any shift
}

int EigenSolver::_size_of_matrix;
//-------------------------------------------------------------//
void EigenSolver::slepc_init(int argc1, char** argv1)
{

  SlepcInitialize(&argc1,&argv1,NULL,NULL);

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
  PC pc;
  


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //print_options(opt);
  

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
 

  if (opt.read_matrix_from_file)
  {
  
  
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,opt.H_file_name.c_str(),FILE_MODE_READ,&viewer);CHKERRQ(ierr); //their

    ierr = MatLoad(viewer,MATAIJ,&A);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

    ierr = MatGetSize(A, &_size_of_matrix, NULL);
   
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,opt.S_file_name.c_str(),FILE_MODE_READ,&viewer);CHKERRQ(ierr); //their

    ierr = MatLoad(viewer,MATAIJ,&B);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  

  }
  
 

  if (opt.matrix_output)
  {//test of the matrix
   

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matA.m",&viewer_out); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer_out,PETSC_VIEWER_ASCII_MATLAB);
    ierr = MatView(A, viewer_out); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer_out);CHKERRQ(ierr);


    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matB.m",&viewer_out); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer_out,PETSC_VIEWER_ASCII_MATLAB);
    ierr = MatView(B, viewer_out); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer_out);CHKERRQ(ierr);

  }


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 
  /* 
     Set operators. In this case, it is a generalized eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  // ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);
  
 
  shift = opt.spectrum_shift;
 
  ierr = EPSSetTolerances(eps,opt.eps_tolerance,opt.eps_max_it);  CHKERRQ(ierr);


  if (opt.solver_type == "arnoldi" || opt.solver_type == "krylovshur" )
  {
    ierr = EPSSetProblemType(eps,EPS_GNHEP);CHKERRQ(ierr);
    
    //ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);

    if (opt.solver_type == "arnoldi")
      ierr = EPSSetType(eps, EPSARNOLDI);
    else
      ierr = EPSSetType(eps, EPSKRYLOVSCHUR);

    
    ierr = EPSGetST(eps,&st); CHKERRQ(ierr);
  
    if (opt.spectral_trans == "folding")
    {
      ierr = STSetType(st,STFOLD); CHKERRQ(ierr);
      ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);
    }
    else
    {
      ierr = STSetType(st,STSINV); CHKERRQ(ierr);
      ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE);CHKERRQ(ierr);
    }
      
    ierr = STSetShift(st, opt.spectrum_shift);CHKERRQ(ierr); 

  
    ierr = STGetKSP(st, &ksp);CHKERRQ(ierr);

    if (opt.st_ksp_type == "bcgsl")
      ierr = KSPSetType( ksp, KSPBCGS);
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

    

    ierr = KSPGetPC( ksp,&pc);

    if (opt.pc_type == "cholesky")
      ierr = PCSetType(pc,PCCHOLESKY);
    else if (opt.pc_type == "jacobi" )
      ierr =  PCSetType(pc,PCJACOBI);
    else if (opt.pc_type == "ilu" )
      ierr =  PCSetType(pc,PCILU);
    else if (opt.pc_type == "composite" )
      ierr =  PCSetType(pc,PCCOMPOSITE);
  

   
   
   
    

    //   ierr = KSPSetTolerances(ksp,1e-10, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
    
    ierr = KSPSetTolerances(ksp,opt.spectrum_inversion_tolerance, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);

  }
  else if (opt.solver_type == "lapack")
  {
    ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSLAPACK);
    ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);
    //ierr = EPSGetST(eps,&st); CHKERRQ(ierr);
    //ierr = STSetShift(st, opt.spectrum_shift);CHKERRQ(ierr); 

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


      ierr = KSPSetType( ksp, KSPBCGS);CHKERRQ(ierr);


      ierr = STGetKSP(st, &ksp);CHKERRQ(ierr);

      if (opt.st_ksp_type == "bcgsl")
	ierr = KSPSetType( ksp, KSPBCGS);
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
       
    

       ierr = KSPGetPC( ksp,&pc);
       
       if (opt.pc_type == "cholesky")
	 ierr = PCSetType(pc,PCCHOLESKY);
       else if (opt.pc_type == "jacobi" )
	 ierr =  PCSetType(pc,PCJACOBI);
       else if (opt.pc_type == "ilu" )
	 ierr =  PCSetType(pc,PCILU);
       else if (opt.pc_type == "composite" )
	 ierr =  PCSetType(pc,PCCOMPOSITE);
       
     
       //rtol, abstol, dtol, maxits
       ierr = KSPSetTolerances(ksp,opt.spectrum_inversion_tolerance, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
       
    }
  }
 


  ierr = do_solve(opt);  
  
  CHKERRQ(ierr);

 
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
 

   if (opt.read_matrix_from_file)
   { 
  
    
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,opt.H_file_name.c_str(),FILE_MODE_READ,&viewer);CHKERRQ(ierr); //their
    ierr = MatLoad(viewer,MATAIJ,&A);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  ierr = MatGetSize(A, &_size_of_matrix, NULL);
  

  }
  
 

  if (opt.matrix_output)
  {//test of the matrix
   

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matA.m",&viewer_out); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer_out,PETSC_VIEWER_ASCII_MATLAB);
    ierr = MatView(A, viewer_out); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer_out);CHKERRQ(ierr);


   

  } 
 

  shift = opt.spectrum_shift;

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
      ierr = KSPSetTolerances(ksp,opt.spectrum_inversion_tolerance, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
   
  }
  else if (opt.solver_type == "lapack")
  {
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSLAPACK);
    ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);
    //ierr = EPSGetST(eps,&st); CHKERRQ(ierr);
    //ierr = STSetShift(st, opt.spectrum_shift);CHKERRQ(ierr); 

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
      ierr = KSPSetTolerances(ksp,opt.spectrum_inversion_tolerance, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);

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


  ierr = VecDestroy(eigen_vector);

}
//-----------------------------------------------------------------------------//
void EigenSolver::set_initial_vector( const std::vector<Complex>& in_vector)
{

  int ierr;

  Vec v0;

  

  MatGetVecs(A,PETSC_NULL,&v0);
  
  PetscScalar y[_size_of_matrix];
  PetscInt ix[_size_of_matrix];

  for (int j= 0; j < _size_of_matrix ; j++)
  {
    y[j]  = in_vector[j];

    ix[j] = j;
  }


  VecSetValues(v0,_size_of_matrix,ix,y,INSERT_VALUES);

 

  EPSSetInitialVector(eps, v0);
  
 
  ierr = VecDestroy(v0);
}


//-----------------------------------------------------------------------------//

int EigenSolver::prepare_slepc()
{
  /* 
     Create eigensolver context
  */
  int ierr;
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);


  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);


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

  int ncv, nconv;
  



  if (opt.ev_number > 8)
    ncv =  4*opt.ev_number;
  else
    ncv = 32;
   

 
  if (ncv > _size_of_matrix) ncv = _size_of_matrix;
 
#if ((SLEPC_VERSION_MAJOR == 2) && (SLEPC_VERSION_MINOR == 3) && \
    (SLEPC_VERSION_SUBMINOR <= 2))
  if (opt.monitor) EPSSetMonitor(eps, EPSDefaultMonitor, PETSC_NULL);
#else
  if (opt.monitor) EPSMonitorSet(eps, EPSMonitorDefault, PETSC_NULL, PETSC_NULL);
#endif

  

#if (SLEPC_VERSION_MAJOR >= 3)
  ierr = EPSSetDimensions(eps,opt.ev_number, ncv, PETSC_DECIDE); CHKERRQ(ierr);
#else
  ierr = EPSSetDimensions(eps,opt.ev_number, ncv); CHKERRQ(ierr);
#endif

   
  ierr = EPSSolve(eps);
 
  //ierr =  EPSGetConverged(eps,&nconv);  CHKERRQ(ierr);
 
 
  return ierr;

}


//--------------------------------------------------------------//
int EigenSolver::init_H_matrix(unsigned int n)
{

  int ierr;

  ierr = MatCreate(PETSC_COMM_WORLD,&A);
 
  CHKERRQ(ierr);
  
   
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
   
  CHKERRQ(ierr);
  
   
  ierr = MatSetFromOptions(A);
     

  CHKERRQ(ierr);
   
  
 

 
  return(ierr);

  
}


//----------------------------------------------------------//
int EigenSolver::init_S_matrix(unsigned int n)
{

  int ierr;

  ierr = MatCreate(PETSC_COMM_WORLD,&B);
  CHKERRQ(ierr);
  
  
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n);
  CHKERRQ(ierr);
  
   
  ierr = MatSetFromOptions(B);
  CHKERRQ(ierr);
   

  return(ierr);
  
}
//-----------------------------------------------------------------//
int EigenSolver::finalize_H_assembly(void)
{
  int ierr;

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
 
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);


  return(ierr);
}
//------------------------------------------------------------------//
int EigenSolver::finalize_S_assembly(void)
{
  int ierr;

  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  

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
//------------------------------------------------------------------------------------//
int  EigenSolver::preallocate_H_matrix(unsigned int matrix_size,  int*  non_zeros)
{

  int ierr;
  
  ierr = MatCreateSeqAIJ (PETSC_COMM_WORLD, matrix_size, matrix_size,0, non_zeros, &A);

  _size_of_matrix = matrix_size;       
  
  return(ierr);
}
//------------------------------------------------------------------------------------//
int  EigenSolver::preallocate_S_matrix(unsigned int matrix_size,  int*  non_zeros)
{
  int ierr;
  
  ierr = MatCreateSeqAIJ (PETSC_COMM_WORLD, matrix_size, matrix_size,0, non_zeros, &B);
 
  return(ierr);
}

//------------------------------------------------------------------------------------//
double  EigenSolver::get_shift(void)
{
  //PetscErrorCode ierr;
  //ST st;
  //PetscScalar shift;             //

  //ierr = EPSGetST(eps, &st);     // ierr = EPSGetST(eps, st); 

  //CHKERRQ(ierr);   

  //ierr = STGetShift(st, &shift); // ierr = STGetShift(*st, &shift);
 
  //CHKERRQ(ierr); 

  //return(real(shift));
  return(shift);

}
//------------------------------------------------------------------------------------//
bool EigenSolver::check_matrices(double tol, bool verbose)
{
  PetscErrorCode ierr;
  PetscTruth is;
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
  PetscTruth is;
  bool ans;

  ans = false;    

  ierr = SlepcIsHermitian(A, &is);

  ans = is==PETSC_TRUE;

  ierr = SlepcIsHermitian(B, &is);

  return (is==PETSC_TRUE) && ans;  

}
