// $Id$

#include "FEMEigenvalueProblem.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "SimulationOptions.h"
#include "EigenSolver.h"
#include "Messages.h"

#include "equation_systems.h"
#include "linear_implicit_system.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
//#include "fe.h"
#include "fe_interface.h"
#include "petsc_matrix.h"

using namespace std;

FEMEigenvalueProblem::FEMEigenvalueProblem(const ModelOptions& options)
 : EigenvalueProblem(options)
{

  es = NULL;

  mesh = NULL;

  system = NULL;

  _hamiltonian_size = 0;

}

//===============================================================//
void FEMEigenvalueProblem::make_new_dofs( )
{
  new_dofs.clear();


  DofMap& dof_map = system->get_dof_map();

  
  number_of_all_dofs  =   dof_map.n_dofs();
  
  new_dofs.resize(number_of_all_dofs);
  const std::set<unsigned int> :: const_iterator  n_begin = dirichlet_dofs.begin();
  const std::set<unsigned int> :: const_iterator  n_end   = dirichlet_dofs.end();
  std::set<unsigned int> :: const_iterator n_it;

  unsigned int number_it = 0;


  for (unsigned int i = 0; i < number_of_all_dofs ; i++)
    {
      if ( !( dof_map.is_constrained_dof(i) ) && (find(n_begin, n_end, i) == n_end))
	{
	    new_dofs[i].independent = true;
	    new_dofs[i].new_number = number_it;
	    number_it++;

	  }
	else
	  {
	  
	    new_dofs[i].independent = false;
	  }
    }

  
  number_of_new_dofs = number_it;

}
//=====================================================================//
void FEMEigenvalueProblem::make_constraints(void)
{
 
  DofMap& dof_map = system->get_dof_map();
  
 
  //----------------------------------------------------------------------------//
  //I recalculate my copy of the dof constraints because I need them explicitely!
  //  my_dof_constraints.clear(); not clear!!!

  // Look at all the variables in the system
  for (unsigned int variable_number=0; variable_number < dof_map.n_variables();
        ++variable_number)
  {
      

    MeshBase::const_element_iterator       elem_it  = mesh->elements_begin();
    const MeshBase::const_element_iterator elem_end = mesh->elements_end(); 
      
    for ( ; elem_it != elem_end; ++elem_it)
      FEInterface::compute_constraints (my_dof_constraints,
					dof_map,
					variable_number,
					*elem_it);
  }
 
  //-----------------------------------------------------------------------//
  //TODO: add periodic boundary conditions constraints

}

//=========================================================================//
//======================================================================//
void  FEMEigenvalueProblem::create_dirichlet_dofs( )
{
  
  

  SimulationEnvironment& se = get_environment(); 

  DofMap& dof_map = system->get_dof_map();

  MeshBase::const_element_iterator it = mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh->active_local_elements_end();

  dirichlet_dofs.clear();

  unsigned int number_of_variables = dof_map.n_variables();
 
  std::vector<unsigned int> dof_indices;

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    { 
      unsigned int  node_id =  elem->node(n);
      
      const Node* nd = elem->get_node(n); 

      Boundary* bd = se.get_boundary(nd); 
      
      //does a node belong to a a dirichlet boundary condition
      
      //   if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) )
       if (  (bd != NULL  ) )
	if ( bd->get_name() == "Dirichlet" || bd->get_name() == "dirichlet" )
	{
	  
	  for (short band = 0 ; band < number_of_variables ; band++)
	  {
	    dof_map.dof_indices (elem, dof_indices,band); 
	    dirichlet_dofs.insert(dof_indices[n]);
	  }
	     
	}

	  
    }
      
  }
  

}

//=======================================================================//
void FEMEigenvalueProblem::apply_dirichlet_at_all_boundaries()
{
  MeshBase::const_element_iterator it = mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh->active_local_elements_end();
 
 

  dirichlet_dofs.clear();

  

  DofMap& dof_map = system->get_dof_map();

  unsigned int number_of_variables = dof_map.n_variables();

  std::vector<unsigned int> dof_indices;

  

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    unsigned int n_sides;

    if ( dim > 1 ) 
      n_sides = elem->n_sides();
    else
      n_sides = elem->n_nodes();
     
	 	
    for (short i = 0; i < n_sides; i++)
    {
	 
      //check if the side is external -------------------------
      
      Elem* el1 = elem->neighbor(i);

      bool side_is_external = false;
      
      if (el1 == NULL) 
      {
	side_is_external = true;

       
      }
      else
      {


	std::vector< const Elem * > active_family;
	
	if ( el1->has_children() )
	{
	  el1->active_family_tree (active_family);
	  
	  if (active_family.size() == 0)
	    side_is_external = true;
	  
	  //TODO
	  // 
	  //	has to be corrected because it may contain active child that does not belong to boundary
	  // 

	}
	else
	{//no children
	  if ( !(el1->active()) )
	    side_is_external = true;
	}



      }
  

      //-------------------------------------------------------
	



      if (   side_is_external   )  
      {
	
	if (dim > 1)
	{//2D/3D
	  for (unsigned int nd = 0; nd < elem->n_nodes(); nd++)
	  {
	    if (elem->is_node_on_side(nd, i))
	    {
	      for (short band = 0 ; band <  number_of_variables; band++)
	      {
		dof_map.dof_indices (elem, dof_indices,band); 
		dirichlet_dofs.insert(dof_indices[nd]);
	      }
	    }
	  }
	  
	}
	else
	{//1D
	  for (short band = 0 ; band < number_of_variables; band++)
	  {
	    dof_map.dof_indices (elem, dof_indices,band);
	    dirichlet_dofs.insert(dof_indices[i]);
	  }
	}
      }      
    }   
  }
}



//=======================================================================//
void FEMEigenvalueProblem::solve_eigen_value_problem(unsigned int ev_number, double st_shift_value)
{

 
  assemble(); //calculate Hamiltonian and S matrix
  
  copy_H_to_solver( );

  if (_haveS)  copy_S_to_solver( );
  
  // for practical reasons, we let ev_number always be even
  if ((ev_number % 2) == 1)
    ev_number += 1;
  
  if (ev_number > _hamiltonian_size)
    throw SolveFailedException("Number of requested eigenvalues is bigger than the Hamiltonian size");

  

  EigenSolver::prepare_slepc();

  EigenSolver::SLEPCoptions slep_opt;

  slep_opt.solver_type = solver_opt.solver;

  slep_opt.H_file_name = "H.out";
    
  slep_opt.S_file_name = "S.out";

  slep_opt.eps_max_it =  solver_opt.max_iteration_number;

  slep_opt.spectral_trans = solver_opt.spectral_trans;

  slep_opt.read_matrix_from_file = false;

  slep_opt.matrix_output = solver_opt.dump_on_file;

  slep_opt.pc_type = solver_opt.preconditioner;

  slep_opt.st_ksp_type = solver_opt.st_ksp_type;

  slep_opt.use_deflation_space =
      get_solver_options().get_option("use_deflation_space", true);
 
 
  slep_opt.monitor = solver_opt.monitor;

  slep_opt.spectrum_inversion_tolerance = solver_opt.spectrum_inversion_tolerance;

  //EigenSolver::check_matrices(1e-10,true);

  slep_opt.eps_tolerance = solver_opt.eigen_solver_tolerance;

  slep_opt.ev_number = ev_number;
  
  slep_opt.spectrum_shift  = st_shift_value;
 
  //std::cout << "  (EFA) Solving using guess (Hartree) " << st_shift_value << endl;

  bool foundall = false;

  while (!foundall)
  {
    int result;
    if (_haveS) 
      result = EigenSolver::eig_value_problem_general(slep_opt);
    else
      result = EigenSolver::eig_value_problem(slep_opt);
      
    if (result !=0 )
      throw SolveFailedException("Eigensolver problem\n");

    foundall = read_SLEPC_solution();

    slep_opt.spectrum_shift = get_new_spectrum_shift();
    slep_opt.ev_number = solver_opt.number_of_eigenstates;
  }
  

  int result = EigenSolver::clear_slepc();
 
  
}


//=====================================================//
void FEMEigenvalueProblem::parse_options()
{
  const ModelOptions& sol_opt = get_solver_options();

  solver_opt.solver = sol_opt.get_option("solver","krylovshur");

  if ( !(solver_opt.solver == "krylovshur" ||
         solver_opt.solver == "arnoldi" ||
         solver_opt.solver == "arpack" ||
         solver_opt.solver == "lapack" ) )
    throw InitFailedException("Invalid solver " +solver_opt.solver);
      

  solver_opt.max_iteration_number = sol_opt.get_option("max_iteration_number",30000);

  solver_opt.eigen_solver_tolerance  = sol_opt.get_option("eigen_solver_tolerance",1e-9);

  solver_opt.spectral_trans = sol_opt.get_option("spectral_transformation","shift_and_invert");

  solver_opt.number_of_eigenstates   = sol_opt.get_option("number_of_eigenstates", 6);

  solver_opt.spectrum_shift = sol_opt.get_option("guess",0.0);

  // only Dirichlet BC works at the moment !
  //solver_opt.Dirichlet_bc_everywhere = sol_opt.get_option("Dirichlet_bc_everywhere", true);
  solver_opt.Dirichlet_bc_everywhere = true;

  solver_opt.monitor = sol_opt.get_option("monitor", false);


  solver_opt.spectrum_inversion_tolerance = sol_opt.get_option("spectrum_inversion_tolerance", 1e-8);

  //cerr <<  solver_opt.Dirichlet_bc_everywhere << "\n";

  {
    std::string  method_name = get_option("discretization_method","FEM");
    if (method_name == "FEM")
      solver_opt.discretization_method = FEM;
    else if (method_name == "BIM")
      solver_opt.discretization_method = BIM;
    else
      throw InitFailedException( "FEMEigenvalueProblem: Incorrect method " + method_name );  

   
  }

  solver_opt.dump_on_file = get_option("dump_HS_on_files",false);

  {

    unsigned int dim = get_mesh().mesh_dimension();
    
    if (dim == 1)
    {
      solver_opt.preconditioner = std::string("redundant");
      solver_opt.st_ksp_type = std::string("preonly");
    }
    else
    {
      solver_opt.preconditioner = std::string("jacobi");
      solver_opt.st_ksp_type = std::string("bcgsl");
    }

  }

  solver_opt.preconditioner =  sol_opt.get_option("pc_type", solver_opt.preconditioner);

  solver_opt.st_ksp_type =  sol_opt.get_option("ksp_type",solver_opt.st_ksp_type);

}

  
//========================================================================================//
void FEMEigenvalueProblem::do_copy_H_to_solver( )
{

 
  int size_matrix = _H_real->n();
  
  PetscMatrix<Number>* H_real_matrix = static_cast<PetscMatrix<Number>* >(_H_real);

  H_real_matrix->close();


  PetscMatrix<Number>* H_imag_matrix = static_cast<PetscMatrix<Number>* >(_H_imag);

  H_imag_matrix->close();


  //----------------------------------------------------------------------------------------------------//

  int non_zeros_number[number_of_new_dofs];

  //preallocate memory for matrix (only for non-parallel calculus)
  for (int row = 0 ; row < size_matrix; row++)
  {
    if (new_dofs[row].independent)
    {
      int ierr = 0;
      const  PetscScalar *petsc_row_vals_real;
      const  PetscInt *petsc_cols_real;
      int n_cols_real = 0;
	  
      const  PetscScalar *petsc_row_vals_imag;
      const  PetscInt *petsc_cols_imag;
      int n_cols_imag = 0;
      
      set<int> real_column, imag_column, complex_column;
      set<int>::iterator com_col_it;
      

      map<int,double> real_values, imag_values;
      map<int, double>::iterator  position;


      insert_iterator<set<int> >  com_ins(complex_column,complex_column.begin() );
      
      ierr = MatGetRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real, &petsc_row_vals_real);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      real_column.clear();
      real_values.clear();
	
      for (int i = 0; i < n_cols_real; i++)
      {
	if (new_dofs[petsc_cols_real[i]].independent && (petsc_row_vals_real[i] != 0.0))
	{
	  real_column.insert(petsc_cols_real[i]);
	  real_values.insert(make_pair(petsc_cols_real[i],petsc_row_vals_real[i] ));
	}
      }

      ierr = MatGetRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
      imag_column.clear();
      imag_values.clear();
      for (int i = 0; i < n_cols_real; i++)
      { 
	if (new_dofs[petsc_cols_imag[i]].independent && (petsc_row_vals_imag[i] != 0.0))
	{
	  imag_column.insert(petsc_cols_imag[i]);
	  imag_values.insert(make_pair(petsc_cols_imag[i],petsc_row_vals_imag[i] ));
	}
      }
      

      set_union(real_column.begin(), real_column.end(), imag_column.begin(), imag_column.end(), com_ins);
	

      non_zeros_number[new_dofs[row].new_number] = complex_column.size();

     
     
      
      
      ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  
	  

    } 
  }

  EigenSolver::preallocate_H_matrix(number_of_new_dofs,  non_zeros_number);
  _hamiltonian_size = number_of_new_dofs;

  //----------------------------------------------------------------------------------------------------//
 
 
  //write data of columns in each row
  
  for (int row = 0 ; row < size_matrix; row++)
    {
      if (new_dofs[row].independent)
      {
	int ierr = 0;
	const  PetscScalar *petsc_row_vals_real;
	const  PetscInt *petsc_cols_real;
	int n_cols_real = 0;
	  
	const  PetscScalar *petsc_row_vals_imag;
	const  PetscInt *petsc_cols_imag;
	int n_cols_imag = 0;
	  
	set<int> real_column, imag_column, complex_column;
	set<int>::iterator com_col_it;


	map<int,double> real_values, imag_values;
	map<int, double>::iterator  position;


	insert_iterator<set<int> >  com_ins(complex_column,complex_column.begin() );

	ierr = MatGetRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	CHKERRABORT(libMesh::COMM_WORLD,ierr);
	real_column.clear();
	real_values.clear();
	
	for (int i = 0; i < n_cols_real; i++)
	{
	  if (new_dofs[petsc_cols_real[i]].independent && (petsc_row_vals_real[i] != 0.0))
	  {
	    real_column.insert(petsc_cols_real[i]);
	    real_values.insert(make_pair(petsc_cols_real[i],petsc_row_vals_real[i] ));
	  }
	}

	ierr = MatGetRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	CHKERRABORT(libMesh::COMM_WORLD,ierr);
	
	imag_column.clear();
	imag_values.clear();
	for (int i = 0; i < n_cols_real; i++)
	{ 
	  if (new_dofs[petsc_cols_imag[i]].independent && (petsc_row_vals_imag[i] != 0.0))
	  {
	    imag_column.insert(petsc_cols_imag[i]);
	    imag_values.insert(make_pair(petsc_cols_imag[i],petsc_row_vals_imag[i] ));
	  }
	}
	

	set_union(real_column.begin(), real_column.end(), imag_column.begin(), imag_column.end(), com_ins);
	
	vector<unsigned int> column_vector;
	vector<Complex> row_values;


	for (com_col_it = complex_column.begin(); com_col_it != complex_column.end(); com_col_it++)
	{
	  int n1 = *com_col_it;

	  double value_r, value_i;

	  //real part------	  
	  position = real_values.find(n1);
	  if (position != real_values.end()) 
	    value_r = position->second;
	  else 
	    value_r = 0.0;
	  
	     

	  //----------------
	  //imag part 
	  position = imag_values.find(n1);
	  if (position != imag_values.end()) 
	    value_i = position->second;
	  else 
	    value_i = 0.0;
	     
	  //----------------

	  column_vector.push_back(new_dofs[n1].new_number);
	  row_values.push_back(Complex(value_r, value_i));
	      
	}
     
	EigenSolver::insert_H_row( new_dofs[row].new_number, column_vector, row_values);
	
	ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	CHKERRABORT(libMesh::COMM_WORLD,ierr);

	ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  
	  

      } 
    }
  //------------------------------------------------------------------------------


  EigenSolver::finalize_H_assembly();

  
}

//============================================================//

void FEMEigenvalueProblem::do_copy_S_to_solver()
{

  int size_matrix = _S_real->n();

  PetscMatrix<double>* p_matrix = static_cast<PetscMatrix<double>* >(_S_real);

  p_matrix->close();


  //----------preallocate memory------------------------------------------------------
  int non_zeros_number[number_of_new_dofs];

  for (int row = 0 ; row < size_matrix; row++)
  {
    if (new_dofs[row].independent)
    {
      int ierr = 0;
      const PetscScalar *petsc_row_vals;
      const PetscInt *petsc_cols;
      int n_cols = 0;

      ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols, &petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      vector<unsigned int> column_vector;
      vector<Complex> row_values;

      non_zeros_number[new_dofs[row].new_number] = 0;

      for (int col = 0; col < n_cols; col++)
      {
	if ((new_dofs[petsc_cols[col]].independent) &&  (petsc_row_vals[col] != 0.0))
	{

	  non_zeros_number[new_dofs[row].new_number]++;


	}
      }

      ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

    }
  }


  EigenSolver::preallocate_S_matrix(number_of_new_dofs,  non_zeros_number);

  //--------------assebmle data--------------------------------------------------------


  for (int row = 0 ; row < size_matrix; row++)
  {
    if (new_dofs[row].independent)
    {
      int ierr = 0;
      const  PetscScalar *petsc_row_vals;
      const PetscInt *petsc_cols;
      int n_cols = 0;

      ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      vector<unsigned int> column_vector;
      vector<Complex> row_values;

      for (int col = 0; col < n_cols; col++)
      {
	if ((new_dofs[petsc_cols[col]].independent) &&  (petsc_row_vals[col] != 0.0))
	{



	  double value = petsc_row_vals[col];
	  double zero = 0.0;

	  column_vector.push_back(new_dofs[petsc_cols[col]].new_number);
	  row_values.push_back(Complex(value, zero));




	}
      }

      EigenSolver::insert_S_row( new_dofs[row].new_number, column_vector, row_values);

      ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

    }
  }

  EigenSolver::finalize_S_assembly();

}

//=======================================================================================/
void FEMEigenvalueProblem::apply_bc()
{
  if (solver_opt.Dirichlet_bc_everywhere)
  {
    apply_dirichlet_at_all_boundaries();
    make_new_dofs();
  }
  else
  {
    create_dirichlet_dofs();
    
    make_constraints(); //creates a copy of them
    
    make_nodes_periodic();
    
    apply_periodic_bc();
    
    make_new_dofs();
  }
  
}
//=================================================================================
void FEMEigenvalueProblem::apply_periodic_bc()
{

  

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.


  PerfLog perf_log ("Periodic bc. Assembly",false);

  DofMap& dof_map = system->get_dof_map();

  unsigned int number_of_variables = dof_map.n_variables();
  
 
  // The dimension that we are running
  //const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving


  vector<unsigned int> uvar(number_of_variables);

  for (unsigned int i = 0; i < number_of_variables; i++) 
  {
    uvar[i] = i;
  }

  unsigned int system_number=system->number();
  
 
  
  FEType fe_type = dof_map.variable_type(uvar[0]);
  
 
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
   

  

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

 
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_component;
  
  const double pos_tol = 1e-10;
  const double func_tol = 1e-10;
 
  //dof_map.print_dof_constraints();  

  for (int i = 0; i < mesh->mesh_dimension(); i++) //Loop over all the mesh directions
    { 
      if  (solver_opt.periodicity[i]) //Check if the periodic b.c. are applied along the direction i
	
	{
	 
	  std :: vector <const Node*>& vec =  nodes_periodic[i];

	  for (unsigned int n = 0; n < vec.size(); n++) // Loop over all the nodes
	    {
	      const Node* node1 = vec[n];
	      
	    
		 
	      for (unsigned int var_index = 0 ; var_index  <  number_of_variables;  var_index ++)
		{//let us find dof for it-----------------
		  

		 // const Node& node = mesh.node(n);
		  
		  const unsigned int  n_dof = node1->dof_number(system_number,uvar[var_index],0);
		
		  //dof is found-------------------------------
		 
		


		  if (! dof_map.is_constrained_dof(n_dof) ) //only if the dof is not constrained do the job
		    {
		      //let us make a  point that lies at the opposite side
		      Point point2(*node1);
		
		      point2(i) = point2(i) + max_coord[i] - min_coord[i];
		  
		      
		      //corresponding point is created
		      
		      
		      //let us find an element this point belongs to and calculate the constraints

		      //the most coarse element first
		      unsigned int refinement_level = 0; 
		      MeshBase::const_element_iterator el3  = mesh->active_elements_begin();
		      MeshBase::const_element_iterator end_el3 = mesh->active_elements_end();
		      
		      const Elem* elem1;
		      bool found = false; 

		      unsigned int el_number = 0;

		      for ( ; ( (el3 != end_el3) ) ; ++el3)  
			{
			  Elem* elem = *el3;
			  
			  if (element_on_boundary(elem))
			    {
			      if (elem->contains_point(point2))
				{
				  elem1 = elem;
				  found = true;
				  break;
				  
				}
			    }
			  el_number++;
			}
		      
		      if (!found)  throw ModelErrorException("EnvelopFunctionApprox: Mesh periproblem");
		  
		      

		      
		      //active elem1 contains the opposite  point, we can constrain it now
		      
		      DofConstraintRow constraint;
		      constraint.clear();
		      
		      dof_map.dof_indices (elem1, dof_indices_component, uvar[var_index]);
		      
		      std::vector<Point> point2_vec(1);
		      
		      point2_vec[0] = point2;
		      
		      std::vector<Point> point2_ref_vec(1);
			  
			  
		      FEInterface::inverse_map (elem1->dim(), fe_type , elem1,  point2_vec,  point2_ref_vec)  ;
		      
		      fe->reinit (elem1, &point2_ref_vec);

		      Point point_temp = point2_ref_vec[0];
		     
		      
		      for (int i1 = 0; i1 < phi.size(); i1++)
			{
			
			  if ( std::abs(phi[i1][0]) >  func_tol )  
			    {
			     
			      constraint[dof_indices_component[i1]] = phi[i1][0];
			    
			    }
			}
		       

		      dof_map.add_constraint_row (n_dof,  constraint); 
                      my_dof_constraints.insert(pair<unsigned int, DofConstraintRow>(n_dof,  constraint));
		      
		    }
		     
		    
		}
	    }
	  
	  
	}
    }  
 
}
//---------------------------------------------------------------------------------------------

void FEMEigenvalueProblem::make_nodes_periodic()
{
  const double pos_tol = 1e-10;
 
  nodes_periodic.clear();

 
  

  for (unsigned dir = 0; dir <=dim-1; dir++)
  {//directions
    std::vector < const Node*> temp_vec;
    temp_vec.clear();
    
    if (solver_opt.periodicity[dir]) 
    {  
      for (unsigned int n = 0; n < mesh->n_nodes(); n++) // Loop over all the nodes
      {
	const Node* node1 = & (mesh->node(n));
	if (node1->active())
	{		
	  if ( std::abs( (*node1)(dir) - min_coord[dir]) < pos_tol)  temp_vec.push_back(node1);
	}
      }
    }
    nodes_periodic.push_back(temp_vec);
  }
}

//--------------------------------------------------------------------------------------//
void FEMEigenvalueProblem::do_init()
{

}

//============================================================//

void FEMEigenvalueProblem::print_H(const std::string& outpath) const
{
  //std::string path = SimulationOptions::scratch_path;

  _H_real->print_matlab(outpath+"/Hr.m");
  _H_imag->print_matlab(outpath+"/Hi.m");
  if (_haveS) _S_real->print_matlab(outpath+"/Sr.m");

}


int FEMEigenvalueProblem::get_H_dim() const 
{ 
  return  _H_real->n();
}
    

int 
FEMEigenvalueProblem::get_H_nnz() const 
{ 
 

  DofMap& dof_map = system->get_dof_map();

  const std::vector<unsigned int>& n_nz = dof_map.get_n_nz();

  unsigned int row_start = _H_real->row_start();
  unsigned int row_stop = _H_real->row_stop();
  unsigned int nnz=0;

  for (unsigned int i = row_start ; i < row_stop; i++)
    nnz += n_nz[i];

  return nnz;

}



void 
FEMEigenvalueProblem::get_H_csr(std::vector<Complex>& A,
                                std::vector<int>& JA,
                                std::vector<int>& IA) const 

{

  PetscMatrix<Number>* H_real_matrix = static_cast<PetscMatrix<Number>* >(_H_real);
  H_real_matrix->close();


  PetscMatrix<Number>* H_imag_matrix = static_cast<PetscMatrix<Number>* >(_H_imag);
  H_imag_matrix->close();

  unsigned int row_start = _H_real->row_start();
  unsigned int row_stop = _H_real->row_stop();
  unsigned int row, col, ind = 0;

  IA[0] = 0;

  for (unsigned int row = row_start ; row < row_stop; row++)
  {
    int ierr = 0;
    const PetscScalar *petsc_row_vals_real;
    const PetscScalar *petsc_row_vals_imag;
    const PetscInt *petsc_cols;
    int n_cols_real = 0;
    int n_cols_imag = 0;
    
    ierr = MatGetRow(H_real_matrix->mat(), row, &n_cols_real, &petsc_cols, &petsc_row_vals_real);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = MatGetRow(H_imag_matrix->mat(), row, &n_cols_imag, &petsc_cols, &petsc_row_vals_imag);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    if (n_cols_real != n_cols_imag) Messages::error("n_cols_real != n_cols_imag");

    for (unsigned int j = 0; j<n_cols_real; j++)
    {
      col = petsc_cols[j];
      
      A[ind] = Complex(petsc_row_vals_real[j], petsc_row_vals_imag[j]); 
      JA[ind] = petsc_cols[j];

      ind++;  
    } 
   
    IA[row+1]= ind;

    ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols, &petsc_row_vals_real);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
    ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols, &petsc_row_vals_imag);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);
 

  }

}



void 
FEMEigenvalueProblem::get_S_csr(std::vector<Complex>& A, 
                                std::vector<int>& JA, 
                                std::vector<int>& IA) const 
{
  PetscMatrix<Number>* S_real_matrix = static_cast<PetscMatrix<Number>* >(_S_real);
  S_real_matrix->close();

  unsigned int row_start = _S_real->row_start();
  unsigned int row_stop = _S_real->row_stop();
  unsigned int row, col, ind = 0;

  IA[0] = 0;

  for (unsigned int row = row_start ; row < row_stop; row++)
  {
    int ierr = 0;
    const PetscScalar *petsc_row_vals;
    const PetscInt *petsc_cols;
    int n_cols = 0;
    
  
    ierr = MatGetRow(S_real_matrix->mat(), row, &n_cols, &petsc_cols, &petsc_row_vals);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);

    for (unsigned int j = 0; j<n_cols; j++)
    {
      col = petsc_cols[j];
      
      A[ind] = Complex(petsc_row_vals[j], 0.0); 
      JA[ind] = petsc_cols[j];
      ind++;  
    } 
   
    IA[row+1]= ind;


    ierr = MatRestoreRow(S_real_matrix->mat(), row ,&n_cols, &petsc_cols, &petsc_row_vals);
    CHKERRABORT(libMesh::COMM_WORLD,ierr);
       

  }


}
