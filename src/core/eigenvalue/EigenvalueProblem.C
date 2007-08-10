#include "EigenvalueProblem.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "EigenSolver.h"

using namespace std;

EigenvalueProblem::EigenvalueProblem(void)
{

  es = NULL;

  mesh = NULL;

  system = NULL;

}

//===============================================================//
void EigenvalueProblem::make_new_dofs( )
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
void EigenvalueProblem::make_constraints(void)
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
void  EigenvalueProblem::create_dirichlet_dofs( )
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
void EigenvalueProblem::apply_diriclet_bc_at_all_boundaries()
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
	      
	      const Point& node = elem->point(nd);
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
void EigenvalueProblem::solve_eigen_value_problem(unsigned int ev_number, double st_shift_value)
{

 
  calculate_Hamiltonian_and_S(); //calculate Hamiltonian and S matrix
 
  EigenSolver::prepare_slepc();

  EigenSolver::SLEPCoptions slep_opt;

  slep_opt.solver_type = solver_opt.solver;

  slep_opt.H_file_name = "H.out";
    
  slep_opt.S_file_name = "S.out";

  slep_opt.eps_max_it =  solver_opt.max_iteration_number;


  slep_opt.read_matrix_from_file = false;

  slep_opt.matrix_output = false;

  slep_opt.preconditioner = solver_opt.preconditioner;

  if (solver_opt.solve_ev_problem_twice)
  {

    st_shift_value = 0.0;


    slep_opt.ev_number = 1;

   
    slep_opt.eps_tolerance =  solver_opt.eigen_solver_tolerance;
   
   
  
    
    slep_opt.spectrum_shift = st_shift_value;


    slep_opt.matrix_output = false;
  
   

    {
      int result;
      if (solver_opt.discretization_method == FEM) 
	result = EigenSolver::eig_value_problem_general(slep_opt);
      else
	result = EigenSolver::eig_value_problem(slep_opt);


      if (result !=0 )
      {
	throw SolveFailedException("Eigensolver problem\n");
      }
    }
   
    st_shift_value = get_new_spectrum_shift();

   

  }


  slep_opt.matrix_output = false;

  
  slep_opt.eps_tolerance = 1e-9;

  slep_opt.ev_number = ev_number;
  
  slep_opt.spectrum_shift  = st_shift_value;
  


  {
      int result;
      if (solver_opt.discretization_method == FEM) 
	result = EigenSolver::eig_value_problem_general(slep_opt);
      else
	result = EigenSolver::eig_value_problem(slep_opt);
      
      

      if (result !=0 )
      {
	cerr << "result of EigenSolver is bad:  " << result << "\n";
	throw SolveFailedException("Eigensolver problem\n");
      }
  }
 



  read_SLEPC_solution(ev_number);

  int result = EigenSolver::clear_slepc();
 
  
}


//=====================================================//
void EigenvalueProblem::parse_options()
{
  const ModelOptions& mod_opt = get_options();

  solver_opt.solver = mod_opt.get_option("solver","arnoldi");

  solver_opt.max_iteration_number = mod_opt.get_option("max_iteration_number",10000);

  solver_opt.eigen_solver_tolerance  = mod_opt.get_option("eigen_solver_tolerance",1e-9);

  solver_opt.solve_ev_problem_twice  = mod_opt.get_option("solve_ev_problem_twice",true);



  solver_opt.number_of_eigenstates   = mod_opt.get_option("number_of_eigenstates", 6);



  solver_opt.Dirichlet_bc_everywhere = mod_opt.get_option("Dirichlet_bc_everywhere",false);


  {
    std::string  method_name = mod_opt.get_option("method","FEM");
    if (method_name == "FEM")
      solver_opt.discretization_method = FEM;
    else if (method_name == "BIM")
      solver_opt.discretization_method = BIM;
    else
      throw InitFailedException( "EigenvalueProblem: Incorrect method " + method_name );  

   
  }

  { 
    std::string prec =  mod_opt.get_option("preconditioner","default");
    if (!(prec == "default" || prec == "cholesky"))
    {
      throw InitFailedException( "EigenvalueProblem: Incorrect preconditioner name" + prec);  
    }
    else
    {
      solver_opt.preconditioner = prec;

    }


  }


}


//========================================================================================//
void EigenvalueProblem::copy_H_matrix_to_solver( )
{

 
  int size_matrix = Ham_real->n();
  

  EigenSolver::init_H_matrix(number_of_new_dofs);

  
  PetscMatrix<Number>* H_real_matrix = static_cast<PetscMatrix<Number>* >(Ham_real);

  H_real_matrix->close();


  PetscMatrix<Number>* H_imag_matrix = static_cast<PetscMatrix<Number>* >(Ham_imag);

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
	

      non_zeros_number[new_dofs[row].new_number] = complex_column.size();

     
     
      
      
      ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  
	  

    } 
  }

  EigenSolver::preallocate_H_matrix(number_of_new_dofs,  non_zeros_number);
  

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


