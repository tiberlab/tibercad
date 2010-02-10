// $Id$

#include "MaxwellEquations.h"
#include "MaxwellPhysicalModel.h"
#include "EigenSolver.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "TiberMath.h"

#include "equation_systems.h"
#include "dense_submatrix.h"

using namespace std;
using namespace Constants;

Device* MaxwellEquations::_device;

MaxwellEquations::MaxwellEquations(const ModelOptions& options)
 : FEMEigenvalueProblem(options)
{

  es = NULL;

  mesh = NULL;

  system = NULL;
}

//=====================================================//
BoundaryProperties* MaxwellEquations::create_boundary_model(const ModelOptions& options) const  throw (ModelErrorException)
{


 
  return NULL;
}

//=======================================================================================================//
void 	MaxwellEquations::build_integrated_quantities (std::vector< double > &values)
{

  values.resize(0);
  if (plot_solution("PhotonEnergy"))
  {
    unsigned int n = solution.size();
    values.resize(n);
    for (unsigned int i = 0; i < n; i++)
    {
    
      
      values[i] = Constants::Hartree * (sqrt(solution[i].k_squared) / (opt.work_units / Constants::bohr_radius)) 
	* (1.0 / Constants::fine_structure_constant) ;
    }
  }



}
//=======================================================================================================//
void 	MaxwellEquations::build_integrated_quantities_description (
						 std::vector< std::string > &legend, 
						 std::vector< std::string > &description)
{
  legend.resize(0);


  if (plot_solution("PhotonEnergy"))
  {
    unsigned int n = solution.size(); 
    legend.resize(n);
    for (unsigned int i = 0; i < n; i++)
    {


      ostringstream temp;
      temp << i ;
      legend[i] = temp.str();
    }



    description.resize(1);
    description[0] = "Photon energy [eV]";

  }



}
//=======================================================================================================//
void MaxwellEquations::build_nodal_results(const std::set<std::string>& variables,
						std::vector<double>& results, std::vector<std::string>& legend)
{

 

  const set<string>::const_iterator varend(variables.end());

  const MeshBase& mesh1 = system->get_mesh();


  const bool output = variables.find("OpticalModes") != varend ;
 
  if (output)
  {

    MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
    const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();

    unsigned int number_of_points = 0;
   
    const unsigned int num_functions = solution.size();
    unsigned int temp = 0;
  
    for ( ; nd != nd_el ; ++nd)  number_of_points++;

    legend.resize(num_functions);
    results.resize(number_of_points*num_functions);
   

    for (unsigned int i = 0; i < num_functions; i++)
    {
      std::ostringstream i_str;
      i_str << "state_number_" << i ; //The states are numbered starting from 0 
      legend[i] = i_str.str();

      std::vector<double> prob_data(number_of_points, 0.0);


      prepare_field_mod_squared(i,  prob_data);

 
      
      for (unsigned int i1 = 0; i1<number_of_points; i1++)
	results[ num_functions * i1 + i] = prob_data[i1];

      
    }
    
    
    
    

  }

 
  

}



//=======================================================================================================//

void MaxwellEquations::prepare_field_mod_squared(const unsigned int mode_number, std::vector<double>& data)
{
 DofMap& dof_map = system->get_dof_map();

 MeshBase::const_node_iterator       nd     = mesh->active_nodes_begin();
 const MeshBase::const_node_iterator nd_el  = mesh->active_nodes_end();
 
 
 unsigned int number_of_points = 0;
 
 for ( ; nd != nd_el ; ++nd)  number_of_points++;

 data.resize(number_of_points, 0.0);

 vector< vector<Complex> >  field_data;
 field_data.resize(number_of_points);

 for (unsigned int i = 0; i < number_of_points; i++) field_data[i].resize( number_of_field_components, Complex(0.0, 0.0) );

 MeshBase::const_element_iterator it = mesh->active_elements_begin();
 const MeshBase::const_element_iterator end =  mesh->active_elements_end();

 std::vector<unsigned int> dof_indices;

 //!calculation of probability function
 for ( ; it != end; ++it)
 {
   const Elem* elem = *it;
      

   for (short psi_index = 0; psi_index < number_of_field_components; psi_index++)
   {
     dof_map.dof_indices (elem, dof_indices, psi_index);
     for (unsigned int n = 0; n < elem->n_nodes(); n++)
     { 
       unsigned int  node_id =  elem->node(n);
       
       Complex value =  (solution[mode_number].eigen_vector[ dof_indices[n] ]);
       
       field_data[node_id][psi_index] = value;

	    
	      
     }
   }
 }

 //done

 double t1;
 for (unsigned int i = 0; i < number_of_points; i++) 
 {
   for (unsigned int j = 0; j < number_of_field_components; j++)
   {

     t1 = std::abs(field_data[i][j]);
     data[i] += t1 * t1; 
   }
    
 }



}




//=======================================================================================================//
PhysicalModel*  MaxwellEquations::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{
  MaxwellPhysicalModel* model = dynamic_cast<MaxwellPhysicalModel*> ( PhysicalModelInterface::create("maxwell", options) );
 
  if (model == NULL)
    throw ModelErrorException("MaxwellEquations: cannot create MaxwellPhysicalMode");
 
  return(model);

}

//=======================================================================================================//

void MaxwellEquations::do_init()
{


  FEMEigenvalueProblem::do_init();


  SimulationEnvironment& si = get_environment();   

  _device = &( si.get_device() );
  


  const ModelOptions& mod_opt = get_options();

  opt.work_units = mod_opt.get_option("length_units", 1e-6);


  opt.scalar_approximation = mod_opt.get_option("scalar_approximation", false);

  if (opt.scalar_approximation)
    number_of_field_components = 1;
  else
    number_of_field_components = 3;


  es = &(get_equation_systems());
  
  mesh = &(es->get_mesh());
  
  system_name = get_equation_system_name ( );

  es->add_system<LinearImplicitSystem> (system_name);

  system = &( es->get_system<LinearImplicitSystem>( system_name ) );

  dim = mesh->mesh_dimension();

 
  system->add_variable("Ax",FIRST);
  if (!opt.scalar_approximation)
 
  {
    system->add_variable("Ay",FIRST);
    system->add_variable("Az",FIRST);
  }



 

  
   
  system->add_matrix("Ham_real"); //add matrix for a real part of the Hamiltonian

  Ham_real = & (system->get_matrix("Ham_real"));


  system->add_matrix("Ham_imag"); //add matrix for a real part of the Hamiltonian

  Ham_imag = & (system->get_matrix("Ham_imag"));



  system->add_matrix("S_real"); //add  matrix for real part of S matrix

  S_real = &( system->get_matrix("S_real") );

  
  system->add_matrix("S_imag"); //add matrix for imaginary part of S matrix

  S_imag = &( system->get_matrix("S_imag") );

  system->init();




   



}

//=======================================================================================================//
void MaxwellEquations::parse_options()
{

  FEMEigenvalueProblem::parse_options();

  const ModelOptions& mod_opt = get_options();




  double spectrum_shift_in_eV   = mod_opt.get_option("spectrum_shift", 0.0);
 
  opt.spectrum_shift    = TiberCad::pow_2(  spectrum_shift_in_eV/( Constants::Hartree) * (opt.work_units / Constants::bohr_radius)
				      * Constants::fine_structure_constant) ;

  solver_opt.solve_ev_problem_twice = mod_opt.get_option("solve_ev_problem_twice",false);


  solver_opt.preconditioner =  mod_opt.get_option("preconditioner","cholesky");


}

//=======================================================================================================//
void MaxwellEquations::do_solve()
{



  parse_options();

  if (solver_opt.Dirichlet_bc_everywhere)
    apply_diriclet_bc_at_all_boundaries();
  else
    create_dirichlet_dofs();
 

  make_constraints(); //creates a copy of them

  make_nodes_periodic();
  
  apply_periodic_bc();


  make_new_dofs();

  solve_eigen_value_problem( solver_opt.number_of_eigenstates);

}


//========================================================================================================//

void MaxwellEquations::read_SLEPC_solution(unsigned int number_of_ev)
{

  unsigned int number_of_converged_solutions;
 
  number_of_converged_solutions = EigenSolver::number_of_converged_eigenvalues();

  //--------------------------------------------------------------------
  //read eigenvalues

  vector<MaxwellEquations::eigen_value>   ev(number_of_converged_solutions);

  for (unsigned ind = 0; ind < number_of_converged_solutions; ind++)
    {
      
      //cerr << EigenSolver::get_eigenvalue(ind) << "\n";

      ev[ind].k_squared =  EigenSolver::get_eigenvalue(ind) + opt.spectrum_shift;
      

      ev[ind].global_number = ind;
        
    }

  
  
  
  //---------------------------------------------------------------------
 
   sort( ev.begin(), ev.end(), compare_eigenvalue);


  

   //let us find the ground state 
   bool finish = false;
   unsigned int ground_state_index = 0;

   for (unsigned int i = 0; (i < number_of_converged_solutions && (!finish) ); i++)
   {
     if (ev[i].k_squared > opt.spectrum_shift)
     {
       ground_state_index = i;
       finish = true;
     }
   }

  


   cerr << "ground_state_index  " << ground_state_index << "\n";

 



   unsigned int solution_size;

   if (number_of_converged_solutions - ground_state_index< number_of_ev)
     solution_size = number_of_converged_solutions - ground_state_index ;
   else
     solution_size = number_of_ev;



   // solution_size = number_of_converged_solutions;
   // ground_state_index = 0;



   {
     MaxwellEquations::eigen_problem_solution temp1;
     temp1.k_squared = 0;
     temp1.eigen_vector.resize(number_of_all_dofs, Complex(0.0, 0.0));

     solution.clear();
     solution.resize(solution_size, temp1);


    
   }

   map<unsigned int, unsigned int>  global_to_sol_index;
   map<unsigned int, unsigned int> :: iterator it;



   for (unsigned int i =  ground_state_index; i < ground_state_index +  solution_size ; i++)
   {
     global_to_sol_index.insert( make_pair(ev[i].global_number, i - ground_state_index )  );
     solution[i - ground_state_index].k_squared = ev[i].k_squared;
     
   }
 


   //--------------------------------------------------------------------
   //read eigenvectors
 
  

   //read solutions - only independent dofs

   //----------------------------------------------------------------------
   for (unsigned int ind = 0; ind < number_of_converged_solutions; ind++)
   {
      
    vector<Complex> temp;
      
    //EigenSolver::get_eigen_vector( ind, temp);
    unsigned int ind1 = 0;
    
    EigenSolver::get_eigen_vector( ind, temp);
   
    it = global_to_sol_index.find(ind);
     
    if (  (it  !=  global_to_sol_index.end()) )
    {

    

      unsigned int solution_number = it->second;
	 
  


 
      //-----------------------------------------------------------------------------
      //put independent dofs in the eigenvectors that may contain also non independent dofs
      for (unsigned j = 0; j < number_of_all_dofs; j++)
      {
	if (new_dofs[j].independent)
	{
	   
	  solution[solution_number].eigen_vector[j] = temp[new_dofs[j].new_number];
	    
	    
	}
      }


      //put constrained dofs
	

   
	
      for (unsigned int j = 0; j < number_of_all_dofs; j++)
      {

    	  
	std::map<unsigned int, DofConstraintRow> :: iterator it;
	  
	it = my_dof_constraints.find(j);
	    

	if (it != my_dof_constraints.end() )
	{
	  
	  DofConstraintRow constr_row = it->second;
	  
	  std::map<unsigned int, Real>::iterator  c =  constr_row.begin();

	  
	  for ( ; c != constr_row.end() ; ++c )
	  {
	    solution[solution_number].eigen_vector[j] += ( c->second ) * solution[solution_number].eigen_vector[(c->first)];
	  }
	    

	     
	}
	   
	  
      }
	    


    }
      
      //------------------------------------------------------------------------

   }

  
 
   //normalization
   for (unsigned int i = 0; i < solution_size; i++)
   {
     const double norm = eigenstate_norm(i);
     
    
     const unsigned int n1 =  solution[i].eigen_vector.size();

      
     for (unsigned int j = 0; j < n1; j++)
       solution[i].eigen_vector[j] /= Complex(norm, 0.0);
   }
  
  



   


}

//=======================================================================//
void MaxwellEquations::calculate_Hamiltonian_and_S(void)
{
  Ham_real->zero();
  Ham_imag->zero();

  S_real->zero();
  S_imag->zero(); 

  vector<unsigned int> fieldvar(number_of_field_components);
 
  fieldvar[0] = system->variable_number("Ax");
  if (!opt.scalar_approximation)
  {
    fieldvar[1] = system->variable_number("Ay");
    fieldvar[2] = system->variable_number("Az");
  }

  DofMap& dof_map = system->get_dof_map();
  
  FEType fe_type = dof_map.variable_type(fieldvar[0]);

  Scaling& scaling = get_scaling();

  scaling.set_length_scaling(opt.work_units);

  scaling.set_calc_mesh_units(_device->get_mesh_units());

  AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );
    
 
  QGauss qrule (dim, SECOND);

  fe -> attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<Point>& q_point = fe->get_xyz();
  
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();


  std::vector<unsigned int> dof_indices_component;
 
  std::vector<unsigned int> dof_indices;

  DenseMatrix<Number> ham_real;
  DenseMatrix<Number> ham_imag;
  DenseMatrix<Number> s_imag;
  DenseMatrix<Number> s_real;

  DenseSubMatrix<Number> ham_real_sub(ham_real);
  DenseSubMatrix<Number> s_imag_sub(s_imag);
  DenseSubMatrix<Number> s_real_sub(s_real);

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


  MaxwellPhysicalModel* model;

  Tensor2Sym eps_real;


  for ( ; el != end_el ; ++el) 
  {//el
    const Elem* elem = *el;

    const ID subdomain = elem->subdomain_id();

    const Material* mat = _device->get_material(subdomain);

    model = dynamic_cast<MaxwellPhysicalModel*> (  mat->get_model(get_id()) ); 
    

    model->get_dielectric_constant()->get_dielectric_real(eps_real);


  


    dof_map.dof_indices (elem, dof_indices); 
    const unsigned int n_dofs   = dof_indices.size();


    ham_real.resize(n_dofs, n_dofs);
    ham_imag.resize(n_dofs, n_dofs);
    s_imag.resize(n_dofs, n_dofs);
    s_real.resize(n_dofs, n_dofs);

    fe->reinit (elem);

    

    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {//qp

      for (short i = 0; i < number_of_field_components; i++)
      {
	 dof_map.dof_indices (elem, dof_indices_component, fieldvar[i]);
	 const unsigned int n_i_dofs = dof_indices_component.size();

	 for (short j = 0; j < number_of_field_components; j++)
	 {
	    ham_real_sub.reposition(fieldvar[i]*n_i_dofs, fieldvar[j]*n_i_dofs, n_i_dofs, n_i_dofs);

	    s_real_sub.reposition(fieldvar[i]*n_i_dofs,   fieldvar[j]*n_i_dofs, n_i_dofs, n_i_dofs);

	    s_imag_sub.reposition(fieldvar[i]*n_i_dofs,   fieldvar[j]*n_i_dofs, n_i_dofs, n_i_dofs);

	   for (unsigned int p1=0; p1<n_i_dofs; p1++)
	   {
	     for (unsigned int p2=0; p2<n_i_dofs; p2++)
	     {
	       double value = 0;
	       complex<double> s_value = Complex(0.0, 0.0);

	       double eps_real_value;
	       if ( i > j)
		 eps_real_value = eps_real( i+1, j+1);
	       else
		 eps_real_value = eps_real( j+1, i+1);
 
	        s_value += JxW[qp] * phi[p1][qp] * phi[p2][qp] * eps_real_value;
		
	      
	       for (short n = 0; n < dim; n++)
		 value +=  JxW[qp] * dphi[p1][qp](n) * dphi[p2][qp](n) * delta_Kronecker(i,j);

	       
	       ham_real_sub(p1,p2) += value;
	     
 
	       s_real_sub(p1,p2) += s_value.real();
	       s_imag_sub(p1,p2) += s_value.imag();
	       
	     }

	   }

	 }
      }




    }



    ham_real.add( - opt.spectrum_shift, s_real);
    ham_imag.add( - opt.spectrum_shift, s_imag);
    
    vector<unsigned int> dof_indices_tmp;

   

    dof_indices_tmp = dof_indices;
    dof_map.constrain_element_matrix(ham_real, dof_indices_tmp);
    Ham_real->add_matrix(ham_real,dof_indices_tmp);


    dof_indices_tmp = dof_indices;
    dof_map.constrain_element_matrix(ham_imag, dof_indices_tmp);
    Ham_imag->add_matrix(ham_imag,dof_indices_tmp);


    dof_indices_tmp = dof_indices;
    dof_map.constrain_element_matrix(s_real, dof_indices_tmp);
    S_real->add_matrix(s_real,dof_indices_tmp);

    dof_indices_tmp = dof_indices;
    dof_map.constrain_element_matrix(s_imag, dof_indices_tmp);
    S_imag->add_matrix(s_imag,dof_indices_tmp);

    

  }
  

 


  copy_H_matrix_to_solver( );
  copy_S_matrix_to_solver( );

 

  // dof_map.print_dof_constraints();
 
}

//---------------------------------------------------------------------//



//---------------------------------------------------------------------// 
void MaxwellEquations::copy_S_matrix_to_solver(void)
{
  int size_matrix = S_real->n();
  

  EigenSolver::init_S_matrix(number_of_new_dofs);

  
  PetscMatrix<Number>* S_real_matrix = static_cast<PetscMatrix<Number>* >(S_real);

  S_real_matrix->close();


  PetscMatrix<Number>* S_imag_matrix = static_cast<PetscMatrix<Number>* >(S_imag);

  S_imag_matrix->close();


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
      
      ierr = MatGetRow(S_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
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

      ierr = MatGetRow(S_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
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

     
     
      
      
      ierr = MatRestoreRow(S_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = MatRestoreRow(S_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  
	  

    } 
  }

  EigenSolver::preallocate_S_matrix(number_of_new_dofs,  non_zeros_number);
  

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

	ierr = MatGetRow(S_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
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

	ierr = MatGetRow(S_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
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
     
	EigenSolver::insert_S_row( new_dofs[row].new_number, column_vector, row_values);
	
	ierr = MatRestoreRow(S_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	CHKERRABORT(libMesh::COMM_WORLD,ierr);

	ierr = MatRestoreRow(S_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  
	  

      } 
    }
  //------------------------------------------------------------------------------


  EigenSolver::finalize_S_assembly();


}

//---------------------------------------------------------------------------------//

double MaxwellEquations::get_new_spectrum_shift( )
{


  double st_shift_value ;

  read_SLEPC_solution(1);
  
  assert(solution.size() == 1);


  st_shift_value = (solution[0].k_squared - opt.spectrum_shift);
    
 
  st_shift_value -= 0.01;
 
  return st_shift_value;

}
//======================================================================//
double  MaxwellEquations::eigenstate_norm(unsigned int state_number)
{
  double  result;
  
  const vector< Complex > &  eigen_vector =  solution[state_number].eigen_vector;
  

  DofMap& dof_map = system->get_dof_map();

  FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

  Scaling& scaling = get_scaling();

  scaling.set_length_scaling(opt.work_units);

  scaling.set_calc_mesh_units(_device->get_mesh_units());


  AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  QGauss qrule (dim, SECOND);

  fe -> attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  std::vector<unsigned int> dof_indices;


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


 
 


  Complex temp(0.0, 0.0);
  Complex eigen_f_value1, eigen_f_value2;

  for ( ; el != end_el ; ++el) 
  {//el
      
    const Elem* elem = *el;
    fe->reinit (elem);
    
    for (short psi_index = 0; psi_index < number_of_field_components; psi_index++)
    {
      dof_map.dof_indices (elem, dof_indices, psi_index);
      const unsigned int n_psi_dofs = dof_indices.size();

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {//qp	      
	for (unsigned int p1=0; p1<n_psi_dofs; p1++)
	{
	  eigen_f_value1 = eigen_vector[dof_indices[p1]];
	  for (unsigned int p2=0; p2<n_psi_dofs; p2++)
	  {
	    eigen_f_value2 = eigen_vector[dof_indices[p2]];
	    temp += ( JxW[qp] * phi[p1][qp] * eigen_f_value1 *  phi[p2][qp] * conj(eigen_f_value2) );
	  }
	}
	      
      }


	
	  
    }
      

  }



  result = sqrt( abs(temp)  );


  return(result);

} 
