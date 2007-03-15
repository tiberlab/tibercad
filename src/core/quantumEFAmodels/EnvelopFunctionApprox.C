#include "SimulationEnvironment.h"
#include "EnvelopFunctionApprox.h"
#include "ModelOptions.h"
#include "EFAbulkModel.h"
#include "Material.h"
#include "Boundary.h"
#include <gnuplot_io.h>

using namespace std;
using namespace Constants;

Device*   EnvelopFunctionApprox:: _device;

//---------------------------------------------------------------------------------//


void EnvelopFunctionApprox::build_nodal_results(const std::set<std::string>& variables,
						std::vector<double>& results, std::vector<std::string>& legend)
{
  const set<string>::const_iterator varend(variables.end());

  const Mesh& mesh1 = system->get_mesh();

  MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
  const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();

  unsigned int number_of_points = 0;

  for ( ; nd != nd_el ; ++nd)  number_of_points++;
  
  if (variables.find("eigen_functions") != varend)
  {

    legend.resize(solution.size());
    results.resize(number_of_points * solution.size());

    unsigned int num_var = solution.size();
    for (unsigned int i = 0; i < num_var; i++)
    {
      std::ostringstream i_str;
      i_str << "state_number_" << i ; //The states are numbered starting from 0 
      legend[i] = i_str.str();

      std::vector<double> prob_data;

      prepare_probability_function(i,  prob_data);

      for (unsigned int i1 = 0; i1<number_of_points; i1++)
	results[num_var * i1 + i] = prob_data[i1];


    }
  }

}

//---------------------------------------------------------------------------------//


 
void EnvelopFunctionApprox::get_integrated_quantities(const std::set<std::string>& names,
						      std::vector<double>& values)
{

} 

//===================================================//

const std::vector< EnvelopFunctionApprox::eigen_propblem_solution >& EnvelopFunctionApprox::get_solution() const
{
  return(solution);
}
  

//====================================================//
PhysicalModel* EnvelopFunctionApprox::create_physical_model(const ModelOptions& options) const throw (ModelErrorException)
{


  EFAbulkModel* model = dynamic_cast<EFAbulkModel*> ( PhysicalModelInterface::create("EFAmodel", options) );
 
  if (model == NULL)
    throw ModelErrorException("EnvelopFunctionApprox: cannot create EFAbulkModel");

  

  return(model);

}
//=====================================================//
BoundaryProperties* EnvelopFunctionApprox::create_boundary_model(const ModelOptions& options) const  throw (ModelErrorException)
{
  
  return NULL;
}


//====================================================//
double EnvelopFunctionApprox::get_band_edge() const
{
   

 
  
  
  DofMap& dof_map = system->get_dof_map();

  FEType fe_type = dof_map.variable_type(0); 

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  QGauss qrule (dim, FIFTH);

  fe -> attach_quadrature_rule (&qrule);

  const std::vector<Point>& q_point = fe->get_xyz();

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


  

  double cond_band_edge;
 
  double valence_band_edge;
  
  vector <double> electric_potential(1, 0.0);

  
  vector<double> band_edges;
 

 
 
  for (; el != end_el ; ++el )
    {
      const Elem* elem = *el;
     
      fe->reinit (elem);

      if (opt.consider_potential)
	{
	  poisson_equation -> get_electric_potential(elem, q_point, electric_potential);
	  
	}
       
      poisson_equation -> get_band_edges ( *el, band_edges );

      short n1 = electric_potential.size();

      if ( elem == *(mesh->active_elements_begin()) )
	{
	  cond_band_edge = band_edges[0] - electric_potential[0];
	  valence_band_edge = band_edges[1] - electric_potential[0];
	}

      for (unsigned int i = 0; i < n1; i++) 
	{     
	
	  if (cond_band_edge > band_edges[0] - electric_potential[i]) 
	    cond_band_edge = band_edges[0] - electric_potential[i];

	  if (valence_band_edge < band_edges[0] - electric_potential[i] ) 
	    valence_band_edge = band_edges[1] - electric_potential[i];
	}
  
    }


  
  
  if (opt.particle == "el")
    return(cond_band_edge);
  else  
    return(valence_band_edge);
  

}

//===================================================//
EnvelopFunctionApprox::EnvelopFunctionApprox()
{
  poisson_equation = NULL;

  strain = NULL;
}
//===================================================//
void EnvelopFunctionApprox::parse_options()
{
 
  const ModelOptions& mod_opt = get_options();



  opt.particle                = mod_opt.get_option("particle", "el");

  opt.periodicity[0]          = mod_opt.get_option("x-periodicity", false);
  opt.periodicity[1]          = mod_opt.get_option("y-periodicity", false);
  opt.periodicity[2]          = mod_opt.get_option("z-periodicity", false);

  opt.log_output              = mod_opt.get_option("log_output", false);

  opt.solver                  = mod_opt.get_option("solver","arnoldi");
  opt.eigen_solver_tolerance  = mod_opt.get_option("eigen_solver_tolerance",1e-9); 
  opt.max_iteration_number    = mod_opt.get_option("max_iteration_number",10000);

  opt.solve_ev_problem_twice  = mod_opt.get_option("solve_ev_problem_twice",true);
  opt.Dirichlet_bc_everywhere = mod_opt.get_option("Dirichlet_bc_everywhere",false);
  

  opt.mpi_command_line        = mod_opt.get_option("mpi_command_line","");
  opt.solver_command_line     = mod_opt.get_option("solver_command_line", 
						   " -eps_largest_magnitude -st_type sinvert -st_ksp_rtol 1e-10 -st_ksp_type bcgs" );

  opt.output_type             = mod_opt.get_option("output_type","GMV");

  opt.eigen_number_increase_factor = mod_opt.get_option("eigen_number_increase_factor",1.2);
  

  opt.number_of_eigenstates   = mod_opt.get_option("number_of_eigenstates", 6);


  opt.spectrum_shift          = mod_opt.get_option("spectrum_shift", 0.0);




  //-------------------------------------------------------------------------------------------//
  //strain model
  std::string  strain_model_name = mod_opt.get_option("strain_model_name","no_strain");

  if (strain_model_name != "no_strain" )
  {
    opt.consider_strain = true;
    
    strain = dynamic_cast< Macrostrain* > ( find_simulation( strain_model_name )   ); 
    
    if (strain == NULL)
      throw InitFailedException("Unknown strain model" + strain_model_name);
    
 
  }
  else
  {
    opt.consider_strain = false;    
  }


  //-------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------//
  //poisson
  std::string  poisson_model_name = mod_opt.get_option("poisson_model_name","no_poisson");
  if ( poisson_model_name != "no_poisson" )
  {
    opt.consider_potential = true;
    
    poisson_equation  = dynamic_cast< DriftDiffusion* > (find_simulation ( poisson_model_name ));
    
    if (poisson_equation == NULL)
      throw InitFailedException( "Unknown poisson model " + poisson_model_name);

  }
  else
  {
    opt.consider_potential = false;
  }
  //--------------------------------------------------------------------------------------------//
  //as default, we  estimate spectrum shift only in electric potential is defined
  opt.estimate_spectrum_shift =  opt.consider_potential;

  opt.estimate_spectrum_shift = mod_opt.get_option("estimate_spectrum_shift",  opt.estimate_spectrum_shift);

  if ( !opt.consider_potential && opt.estimate_spectrum_shift) throw InitFailedException( "EnvelopeFunctionApprox: cannot estimate spectrum shift without electric potential");

  //--------------------------------------------------------------------------------------------//

  if (opt.estimate_spectrum_shift) opt.spectrum_shift = get_band_edge();
  //--------------------------------------------------------------------------------------------//
  //k-vector
  std::vector<double> k_vec(3, 0.0);
  mod_opt.get_option("k_vector",k_vec);
  if (k_vec.size() == 3)
  {
   
    for (short i = 0; i < 3; i++) k_vector[i] = k_vec[i]; 
    
   
  }
  else
    throw InitFailedException( "EnvelopeFunctionApprox: k_vector size must be equal to 3 instead of " + k_vec.size()); 

  //--------------------------------------------------------------------------------//
  std::string  job_name = mod_opt.get_option("job","eigenstates");
  if (job_name == "eigenstates")
    opt.job = EIGENSTATES;
  else if (job_name == "density")
    opt.job = DENSITY;
  else
    throw InitFailedException( "EnvelopeFunctionApprox: Incorrect job " + job_name );  
  //---------------------------------------------------------------------------------//

  opt.Temperature = mod_opt.get_option("Temperature", 300.0);

  opt.initial_eigenstates_number = mod_opt.get_option("initial_eigenstates_number", 6);

}




//===================================================//


//===================================================//
void EnvelopFunctionApprox::do_init( )
{

 


  SimulationEnvironment& si = get_environment();   

  _device = &( si.get_device() );


  double mesh_units = get_environment().get_device().get_mesh_units();
 
  opt.length_scale = mesh_units/bohr_radius;

  es = &(get_equation_systems());

  mesh = &(es->get_mesh());
  
  system_name = get_equation_system_name ( );

  es->add_system<LinearImplicitSystem> (system_name);

  system = &( es->get_system<LinearImplicitSystem>( system_name ) );

  dim = mesh->mesh_dimension();

  //--------------------------------------------------------------------------------------------------------//
  opt.number_of_bands = calculate_number_of_bands( );
 
  
  //--------------------------------------------------------------------------------------------------------//
  //add variables
  psi_name.clear();
  for (short i = 0; i < opt.number_of_bands; i++)
    {
      std::ostringstream var_str;
      var_str << "psi" << i ;
      string name = var_str.str();
      psi_name.push_back(name);

      system->add_variable(name,FIRST);
    } 

  //add matrixes
  //---------------------------------------------------------------------------------------------------------//
   
  DofMap& dof_map = system->get_dof_map();
   
  system->add_matrix("Ham_real"); //add matrix for a real part of the Hamiltonian

  Ham_real = & (system->get_matrix("Ham_real"));

  system->add_matrix("Ham_imag");//add matrix for an imaginary part of the Hamiltonian

  Ham_imag = &(  system->get_matrix("Ham_imag") );

  system->add_matrix("S_real"); //add matrix for S matrix

  S_real = &( system->add_matrix("S_real") );
  
 

  //---------------------------------------------------------------------------------------------------------//
  //My Jacobian 

   my_Jacobian = 1.0;
   for (short i = 1; i <= dim; i++)
     my_Jacobian *= opt.length_scale;

   //--------------------------------------------------------------------------------------------------------//
   MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
   const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

   bool temp = true;

   for ( ; el != end_el ; ++el) 
    {
      const Elem* elem = *el;
      short n1 = elem->n_nodes();
      for (short i1 = 0; i1 < n1 ; i1++)
	{


	  const Point& p = elem->point(i1);
	  for (unsigned i = 0; i < 3; i++)
	  {

	    if (temp)
	      {
		min_coord[i] = p(i);
		max_coord[i] = p(i);
		temp = false;
	      }


	    if (min_coord[i] < p(i)) min_coord[i] = p(i);
	    if (max_coord[i] > p(i)) max_coord[i] = p(i);
	  
	  }
	  
	}

    }
 
   //-------------------------------------------------------------------------------------------------------//

   system->init();
  
 
 
   //------------------------------------------------------------------------------------------------------//
   
  
   

   //------------------------------------------------------------------------------------------------------//


}
//===========================================================//
void EnvelopFunctionApprox::do_solve()
{
  parse_options();

 
  if ( opt.job ==   EIGENSTATES )
    solve_eigen_value_problem( opt.number_of_eigenstates);
  else if ( opt.job == DENSITY )
    calculate_convergent_density(opt.Temperature);
    
  



}
//===========================================================//
void EnvelopFunctionApprox::calculate_Hamiltonian_and_S(void)
{




  Ham_real->zero();
  Ham_imag->zero();
  S_real->zero();

  //material list
  //assemble_material_list();

 vector<unsigned int> psivar(opt.number_of_bands);
 //get numbers of variables
 for (unsigned int i = 0; i < opt.number_of_bands; i++)
   {
     psivar[i] = system->variable_number(psi_name[i]);
   }

 DofMap& dof_map = system->get_dof_map();

 FEType fe_type = dof_map.variable_type(psivar[0]); //all the variable have the same FE representation

 AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);
  
  // Tell the finite element object to use our quadrature rule.

  fe -> attach_quadrature_rule (&qrule);

 // Here we define some references to cell-specific data that
 // will be used to assemble the linear system.
 //
 // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  //------------------------------------------------------------
 std::vector<unsigned int> dof_indices_component;
 
 std::vector<unsigned int> dof_indices;

  //-------------------------------------------------------------
  //matrixes to built the system

  DenseMatrix<Number> ham_real; 
  DenseMatrix<Number> ham_imag;
  DenseMatrix<Number> s_real;


  DenseSubMatrix<Number> ham_real_sub(ham_real);
  DenseSubMatrix<Number> ham_imag_sub(ham_imag);
  DenseSubMatrix<Number> s_real_sub(s_real);


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();
 
  Tensor2Sym strain_crystal_system(0);
  double electric_potential = 0;

  EFAbulkHamiltonian* element_hamiltonian;

  unsigned int el_number = 0;



  for ( ; el != end_el ; ++el) 
    {//el
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      const ID subdomain = elem->subdomain_id();
      const Material* mat = _device->get_material(subdomain);

      element_hamiltonian = (  dynamic_cast<EFAbulkModel*> (  mat ->get_model(get_id()) )  )->get_Hamiltonian_model(); 

      element_hamiltonian->set_k_vector(k_vector);
      element_hamiltonian->calculate_Hamiltonian_k_par();


      dof_map.dof_indices (elem, dof_indices); 
      const unsigned int n_dofs   = dof_indices.size();
      fe->reinit (elem);

      ham_real.resize(n_dofs, n_dofs);
      ham_imag.resize(n_dofs, n_dofs);
      s_real.resize(n_dofs, n_dofs);

      // complex<double> operator_sign = Complex(0.0, -1.0);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{//qp
	  //--------------------------------------------------------------------------------
	  /*
	    We assume that strain and electric potential may be different for different quadrature points
	    It is done for a sake of a multiscale generalization
	  */
	  if (opt.consider_strain) 
	    {
             
	      strain_crystal_system = strain->get_strain_crystal(elem , q_point[qp]);
	    }
 

	  if (opt.consider_potential)
	    {

	      electric_potential = poisson_equation -> get_electric_potential(elem, q_point[qp]);
	     
	    }





	  element_hamiltonian->apply_strain_and_potential(strain_crystal_system, electric_potential);
	  
	  //------------------------------------------------------------------------------------------


	  std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&  
	    model_Ham = ( element_hamiltonian->get_Hamiltonian() );
	  

	

	  for (unsigned int band1 = 0; band1 < opt.number_of_bands; band1++)
	    {//band1
	      dof_map.dof_indices (elem, dof_indices_component, psivar[band1]);
	      const unsigned int n_psi_dofs = dof_indices_component.size();
 
	       for (unsigned int band2 = 0; band2 < opt.number_of_bands; band2++)
		 {//band2
		  
		   //Hamiltonian
		  

		   ham_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
		   ham_imag_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
		   for (unsigned int p1=0; p1<n_psi_dofs; p1++)
		     {
		       for (unsigned int p2=0; p2<n_psi_dofs; p2++)
			 {		      
			   complex<double> value = (0.0, 0.0);
			   //constant
			   value += JxW[qp] * phi[p1][qp] * phi[p2][qp] * model_Ham[band1][band2].constant ;
			  
			  
			   //linear left
			
			   for (short i = 0; i < dim; i++)
			     value -= JxW[qp]* dphi[p1][qp](i) * phi[p2][qp] * model_Ham[band1][band2].linear_left[i]
			       * Complex(0.0, -1.0) /opt.length_scale;
 
                           //linear right

			   for (short i = 0; i < dim; i++)
			     value += JxW[qp]* dphi[p2][qp](i) * phi[p1][qp] * model_Ham[band1][band2].linear_right[i] 
			       * Complex(0.0, -1.0) /opt.length_scale;
			
			   //quadratic
			   
			   for (short i = 0; i < dim; i++)
			     for (short j = 0; j < dim; j++)
			       {
				 value -= JxW[qp] * dphi[p1][qp](i) * dphi[p2][qp](j)*model_Ham[band1][band2].quad[i][j]
				 * Complex(0.0,-1.0) * Complex(0.0, -1.0) /opt.length_scale / opt.length_scale;

			
				  
			       }
			   

			    value *= my_Jacobian;


			  
			    ham_real_sub(p1,p2) += value.real(); 
			    ham_imag_sub(p1,p2) += value.imag(); 

			 }
		     }
		   
		  

		  
		   //S-matrix
		   if (band1 == band2)
		     {
		       s_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
		       for (unsigned int p1=0; p1<n_psi_dofs; p1++)
			 {
			   for (unsigned int p2=0; p2<n_psi_dofs; p2++)
			     {
			       s_real_sub(p1,p2) += JxW[qp] * phi[p1][qp] * phi[p2][qp] * my_Jacobian;
			     }
			 }
		      
		     }
		   //--------------------------------------------------------------------------//

		 }
	    }
	  


	}

  

      ham_real.add( - opt.spectrum_shift/Hartree, s_real);//apply spectrum shift.


/*
      ham_real.scale(Hartree);
      ham_imag.scale(Hartree);
*/
      vector<unsigned int> dof_indices_tmp;

      dof_indices_tmp = dof_indices;

      dof_map.constrain_element_matrix(s_real, dof_indices_tmp);


      S_real->add_matrix(s_real,dof_indices_tmp);

      dof_indices_tmp = dof_indices;

      dof_map.constrain_element_matrix(ham_real, dof_indices_tmp);
      Ham_real->add_matrix(ham_real,dof_indices_tmp);

      dof_indices_tmp = dof_indices;
     
      dof_map.constrain_element_matrix(ham_imag, dof_indices_tmp);
      Ham_imag->add_matrix(ham_imag,dof_indices_tmp);

      el_number++;

    }

 
//this is only to test

   Ham_real->print_matlab("ham_r_matlab.m");
   Ham_imag->print_matlab("ham_i_matlab.m");
   S_real->print_matlab("s.m");



  save_S_matrix("S.out");
  save_H_matrix("H.out");

  //  dof_map.print_dof_constraints();
     
}
//============================================================//

void EnvelopFunctionApprox::save_S_matrix(const std::string & fname)
{
 
  unsigned int out_int;
  unsigned long long out_long_long;

  std::ofstream out (fname.c_str());
  assert (out.good());
   
  const short int_size = sizeof(int);

  const short double_size = sizeof(double);
 // unsigned int MAT_FILE_COOKIE = 1211216;
  
  int magic_number = 1211216;
   
  out_int = *(reinterpret_cast<unsigned int*> (& magic_number) );  endian_swap(out_int);
  out.write( reinterpret_cast<char *>(  &  out_int ), int_size);


 

  int size_matrix = S_real->m() ;


  out_int = *(reinterpret_cast<unsigned int*> (& number_of_new_dofs) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);

  out_int = *(reinterpret_cast<unsigned int*> (& number_of_new_dofs) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
 
  PetscMatrix<Number>* p_matrix = static_cast<PetscMatrix<Number>* >(S_real);

  p_matrix->close();
  //-------------------------------------------------------------------------------------
  //how many non-zero elements do we have?

  int Number_of_elements = 0;

  

  for (int row = 0 ; row < size_matrix; row++)
    {
      if (new_dofs[row].independent)
	{
	  int ierr = 0;
	  const  PetscScalar *petsc_row_vals;
	  const  PetscInt *petsc_cols;
	  int n_cols = 0;
	  
	  ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

       
	  for (int i1 = 0; i1 < n_cols; i1++)
	    if ((new_dofs[petsc_cols[i1]].independent) &&  (petsc_row_vals[i1] != 0.0))  Number_of_elements++;
	    
	
      

	  ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	}
    }

  // std::cout << "We have got " << Number_of_elements << " non-zero elements in S matrix\n"; 

  out_int = *(reinterpret_cast<unsigned int*> (& Number_of_elements) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
  //-----------------------------------------------------------------------------------
  //write number of columns in each row
  for (int row = 0 ; row < size_matrix; row++)
    { 
      if (new_dofs[row].independent)
	{
	  int ierr = 0;
	  const  PetscScalar *petsc_row_vals;
	  const  PetscInt *petsc_cols;
	  int n_cols = 0;
	  int n_new_cols = 0;

	  ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  for (int i1 = 0; i1 < n_cols; i1++) 
	    if ((new_dofs[petsc_cols[i1]].independent) &&  (petsc_row_vals[i1] != 0.0)) n_new_cols++ ;

	  out_int = *(reinterpret_cast<unsigned int*> (& n_new_cols) );  endian_swap(out_int);
	  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
	  
	  

	  ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}
    }
  //===============================================================//
  //write number of columns in each row
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
	  
	  for (int col = 0; col < n_cols; col++)
	    {
	      if ((new_dofs[petsc_cols[col]].independent) &&  (petsc_row_vals[col] != 0.0))
		{
		  int col_number = new_dofs[petsc_cols[col]].new_number;
		  out_int = *(reinterpret_cast<unsigned int*> (& col_number) );  endian_swap(out_int);
		  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
		}
	    }
	  
	  ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  
	}
    }

  //write data
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

	  for (int col = 0; col < n_cols; col++)
	    {
	      if ((new_dofs[petsc_cols[col]].independent) &&  (petsc_row_vals[col] != 0.0))
		{
	      
		  double value = petsc_row_vals[col];
		  double zero = 0.0;
	  
		  out_long_long = *(reinterpret_cast<unsigned long long*> (& value) );  endian_swap(out_long_long);
		  out.write(  reinterpret_cast<char *>( & out_long_long ), double_size);
		  
		  out_long_long = *(reinterpret_cast<unsigned long long*> (& zero) );  endian_swap(out_long_long);
		  out.write(  reinterpret_cast<char *>( & out_long_long ), double_size);
		}
	    }
	  
	  ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	}
    }

}


//============================================================//

void EnvelopFunctionApprox::save_H_matrix(const std::string & fname)
{

  unsigned int out_int;
  unsigned long long out_long_long;

  std::ofstream out (fname.c_str());
  assert (out.good());
   
  const short int_size = sizeof(int);

  const short double_size = sizeof(double);
 // unsigned int MAT_FILE_COOKIE = 1211216;
  
  int magic_number = 1211216;
  
  out_int = *(reinterpret_cast<unsigned int*> (& magic_number) );  endian_swap(out_int);
  out.write( reinterpret_cast<char *>(& out_int), int_size);

  int size_matrix = Ham_real->n();

  

  out_int = *(reinterpret_cast<unsigned int*> (& number_of_new_dofs) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);

  out_int = *(reinterpret_cast<unsigned int*> (& number_of_new_dofs) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);

  
  
  PetscMatrix<Number>* H_real_matrix = static_cast<PetscMatrix<Number>* >(Ham_real);

  H_real_matrix->close();


  PetscMatrix<Number>* H_imag_matrix = static_cast<PetscMatrix<Number>* >(Ham_imag);

  H_imag_matrix->close();


  //-------------------------------------------------------------------------------------//
  //how many non-zero elements do we have?

  int Number_of_elements = 0;

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
	  //  set<int>::iterator com_col_it, com_col_end;
	  insert_iterator<set<int> >  com_ins(complex_column,complex_column.begin() );
	  
	  ierr = MatGetRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  real_column.clear();
	  
	  for (int i = 0; i < n_cols_real; i++) 
	    if  ((new_dofs[petsc_cols_real[i]].independent) && (petsc_row_vals_real[i] != 0.0)) real_column.insert(petsc_cols_real[i]);
	  
	  ierr = MatGetRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  
	  imag_column.clear();
	  for (int i = 0; i < n_cols_real; i++) 
	    if  ( (new_dofs[petsc_cols_imag[i]].independent) &&  (petsc_row_vals_imag[i] !=0.0) ) imag_column.insert(petsc_cols_imag[i]);
	  
	  set_union(real_column.begin(), real_column.end(), imag_column.begin(), imag_column.end(), com_ins);
	  
	  Number_of_elements += complex_column.size();
	  
	  ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}

    }
  if (opt.log_output) {
    std::cout << "We have got " << Number_of_elements << " non-zero elements in the Hamiltonian matrix\n"; 
    std::cout.flush();
  }

  out_int = *(reinterpret_cast<unsigned int*> (& Number_of_elements) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
  //----------------------------------------------------------------------------------------------------//
  //write number of columns in each row
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
	  //  set<int>::iterator com_col_it, com_col_end;
	  insert_iterator<set<int> >  com_ins(complex_column,complex_column.begin() );
	  
	  ierr = MatGetRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  real_column.clear();

	  for (int i = 0; i < n_cols_real; i++) 
	    if (new_dofs[petsc_cols_real[i]].independent  && ( petsc_row_vals_real[i] != 0.0 ) ) real_column.insert(petsc_cols_real[i]);
	  
	  ierr = MatGetRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  
	  imag_column.clear();
	  for (int i = 0; i < n_cols_real; i++) 
	    if (new_dofs[petsc_cols_imag[i]].independent && (petsc_row_vals_imag[i] != 0.0) ) imag_column.insert(petsc_cols_imag[i]);
	  
	  set_union(real_column.begin(), real_column.end(), imag_column.begin(), imag_column.end(), com_ins);

	  int n1 = complex_column.size();

	  out_int = *(reinterpret_cast<unsigned int*> (& n1) );  endian_swap(out_int);
	  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
	  
	  
	  ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}

    }    
  //-----------------------------------------------------------------------------------------------------//
  //write number of columns in each row
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
	  insert_iterator<set<int> >  com_ins(complex_column,complex_column.begin() );
	  
	  ierr = MatGetRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  real_column.clear();

	  for (int i = 0; i < n_cols_real; i++) 
	    if (new_dofs[petsc_cols_real[i]].independent && (petsc_row_vals_real[i] != 0.0) )
	      real_column.insert(new_dofs[petsc_cols_real[i]].new_number);

	  ierr = MatGetRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  imag_column.clear();

	  for (int i = 0; i < n_cols_real; i++) 
	    if (new_dofs[petsc_cols_imag[i]].independent && (petsc_row_vals_imag[i] != 0.0) ) 
	      imag_column.insert(new_dofs[petsc_cols_imag[i]].new_number);

	  set_union(real_column.begin(), real_column.end(), imag_column.begin(), imag_column.end(), com_ins);

   
	  for (com_col_it = complex_column.begin(); com_col_it != complex_column.end(); com_col_it++)
	    {
	      int n1 = *com_col_it;
	      out_int = *(reinterpret_cast<unsigned int*> (& n1) );  endian_swap(out_int);
	      out.write(  reinterpret_cast<char *>( & out_int ), int_size);
	    }
     

    
	  ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
      

	} 
    }
  //-------------------------------------------------------------------------------
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
	  
   
	  for (com_col_it = complex_column.begin(); com_col_it != complex_column.end(); com_col_it++)
	    {
	      int n1 = *com_col_it;
	      double value;

	      //real part------	  
	      position = real_values.find(n1);
	      if (position != real_values.end()) 
		value = position->second;
	      else 
		value = 0.0;

	

	      out_long_long = *(reinterpret_cast<unsigned long long*> (& value) );  endian_swap(out_long_long);
	      out.write(  reinterpret_cast<char *>( & out_long_long ), double_size);
	  
	      //----------------
	      //imag part 
	      position = imag_values.find(n1);
	      if (position != imag_values.end()) 
		value = position->second;
	      else 
		value = 0.0;
	      
	      out_long_long = *(reinterpret_cast<unsigned long long*> (& value) );  endian_swap(out_long_long);
	      out.write(  reinterpret_cast<char *>( & out_long_long ), double_size);
	      
	      //----------------
	    }
     

    
	  ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  

	} 
    }
  //------------------------------------------------------------------------------
}
//===============================================================================//

void EnvelopFunctionApprox::solve_eigen_value_problem(unsigned int ev_number, double st_shift_value)
{
  if (opt.Dirichlet_bc_everywhere)
    apply_diriclet_bc_at_all_boundaries();
   else
      create_dirichlet_dofs();
 

  make_constraints(); //creates a copy of them

  
  make_nodes_periodic();
  
  apply_periodic_bc();

  make_new_dofs();

  calculate_Hamiltonian_and_S(); //calculate Hamiltonian and S matrix
 

  if (opt.solve_ev_problem_twice)
    {

      st_shift_value = 0.0;

      std::ostringstream  command_line_0;

      command_line_0 <<     opt.mpi_command_line   ;

      command_line_0 <<  "  eigen_solver  -f1 H.out   -f2 S.out "; 
      
      command_line_0 <<  "   -eps_nev     " << 1;
    
      command_line_0 <<  "   -eps_type    " << opt.solver <<"  ";
      
      command_line_0 <<  opt.solver_command_line << "  ";

      command_line_0 <<   " -eps_tol  1e-5  -eps_max_it 3000 ";
 
      command_line_0 <<  "    \n";

      if (opt.log_output) 
      {
	cout << command_line_0.str()<<"\n";
	cout.flush();
      }

      std::system( (command_line_0.str()).c_str());

      read_SLEPC_solution(1);

      assert(solution.size() == 1);

      

      st_shift_value = (solution[0].eigen_energy - opt.spectrum_shift)/Hartree;

      if (opt.particle == "el")
	 st_shift_value -= 0.01/Hartree;
      else
	st_shift_value += 0.01/Hartree;
    }

  std::ostringstream  command_line;

  command_line << opt.mpi_command_line   ;
  
  command_line <<  "  eigen_solver  -f1 H.out   -f2 S.out "; 
  
  command_line <<  "   -eps_nev     " << ev_number;
  
  command_line <<  "   -eps_type    " << opt.solver ;

  command_line << " -st_shift  " <<  st_shift_value ;

  command_line <<  opt.solver_command_line << "  ";

  command_line <<  "    \n";



  if (opt.log_output) 
    {
      cout << command_line.str()<<"\n";
      cout.flush();
    }



  std::system( (command_line.str()).c_str());

  read_SLEPC_solution(ev_number);

 

  
  

}

//=============================================================//
void EnvelopFunctionApprox::read_SLEPC_solution(unsigned int number_of_ev )
{//
  /*
  1) Read all eigenvalues
  2) Sort the eigenvalues and select those we want 
  3) Read eigenvectors that correspond to the eigenvalues we want
  4) normalize eigenfunctions
  5) calculate fermi energy for each state
  */


  const short int_size = sizeof(int);

  const short double_size = sizeof(double);

  string fname_eigvals = "eigvals_SLEPC.out";


  std::ifstream file_eigvals ( fname_eigvals.c_str() );

  assert (file_eigvals.good());

  char buffer[int_size];
  char buffer_double[double_size];
  unsigned int dummy;
  
  unsigned long long fict;

  file_eigvals.read(buffer, int_size);
  //--------------------------------------------------------------------
  //how many solutions do we have from SLEPC?
  unsigned int number_of_converged_solutions;
  file_eigvals.read(buffer, int_size); 
  number_of_converged_solutions =  *(reinterpret_cast<unsigned int*> ( buffer));  endian_swap(number_of_converged_solutions);

#ifdef DEBUG
  cout << " Number of converged solutions  " << number_of_converged_solutions << "\n";
#endif

  unsigned int solution_size;
  if (number_of_converged_solutions < number_of_ev)
    solution_size = number_of_converged_solutions;
  else
    solution_size = number_of_ev;



  {
    EnvelopFunctionApprox::eigen_propblem_solution temp1;
    temp1.eigen_energy = 0;
    temp1.Fermi_energy = 0;
    temp1.eigen_vector.resize(number_of_all_dofs, Complex(0.0, 0.0));
    solution.resize(solution_size, temp1);
  }
 

  vector<EnvelopFunctionApprox::eigen_energy>   ev(number_of_converged_solutions);
  
  

  //--------------------------------------------------------------------
  //read eigenvalues
  for (unsigned ind = 0; ind < number_of_converged_solutions; ind++)
    {
      
      
      file_eigvals.read(buffer_double, double_size);
      fict = *( reinterpret_cast<unsigned long long*>( buffer_double) ); endian_swap(fict);


      

     
	
      ev[ind].energy = *(  reinterpret_cast<double*>( &fict ) ) * Hartree + opt.spectrum_shift;
      ev[ind].global_number = ind;

     

      //read dummy imaginary part
      file_eigvals.read(buffer_double, double_size);
      
      
    }
  
  //--------------------------------------------------------------------

  
  //sorting of the solutions
  if (opt.particle == "el") sort( ev.begin(), ev.end(), compare_eigen_energy_electrons1);

  if (opt.particle == "hl") sort( ev.begin(), ev.end(), compare_eigen_energy_holes1 );
   


  map<unsigned int, unsigned int>  global_to_sol_index;
  map<unsigned int, unsigned int> :: iterator it;
  for (unsigned int i = 0; i < solution_size; i++)
    {
      global_to_sol_index.insert( make_pair(ev[i].global_number, i)  );
      solution[i].eigen_energy = ev[i].energy;
    } 
  //--------------------------------------------------------------------
  //read eigenvectors
  string fname_eigvects = "eigvects_SLEPC.out";


  std::ifstream file_eigvects ( fname_eigvects.c_str() );

  assert (file_eigvects.good());

  //read solutions - only independent dofs
  //----------------------------------------------------------------------
  for (unsigned int ind = 0; ind < number_of_converged_solutions; ind++)
    {

     
 
	


      file_eigvects.read(buffer, int_size);
      file_eigvects.read(buffer, int_size);
      unsigned int vector_size =  *(reinterpret_cast<unsigned int*> ( buffer));  endian_swap(vector_size);
      assert( vector_size == number_of_new_dofs);
    
    

      vector<Complex> temp(vector_size); //only independent dofs
     
      for (unsigned j = 0; j < vector_size; j++)
	{
          double re, im;
	  file_eigvects.read(buffer_double, double_size);
	  fict = *( reinterpret_cast<unsigned long long*>( buffer_double) ); endian_swap(fict);
	  re   = *(  reinterpret_cast<double*>( &fict ) );

	  file_eigvects.read(buffer_double, double_size);
	  fict = *( reinterpret_cast<unsigned long long*>( buffer_double) ); endian_swap(fict);
	  im   = *(  reinterpret_cast<double*>( &fict ) );
	  temp[j] = Complex(re,im);
	 

	 
	}

      it = global_to_sol_index.find(ind);
      if (  it  !=  global_to_sol_index.end() )
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
	  for (unsigned j = 0; j < number_of_all_dofs; j++)
	    {

	      std::map<unsigned int, DofConstraintRow> :: iterator it;
	      it = my_dof_constraints.find(j);
	      if (it != my_dof_constraints.end() )
		{
		  DofConstraintRow constr_row = it->second;
		  std::map<unsigned int, Real>::iterator  c =  constr_row.begin();
		
		  for ( ; c != constr_row.end() ; ++c )
		    {
		      solution[solution_number].eigen_vector[j] +=  solution[solution_number].eigen_vector[(c->first)]
			* (c->second);
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
	solution[i].eigen_vector[j] /= norm;
    }

  //Fermi energy calculation

  if (poisson_equation != NULL)
    for (unsigned int i = 0; i < solution_size; i++)
      {
	solution[i].Fermi_energy = calculate_fermi_averaged( i);
	
      }
 

 


  
}
//=============================================================//


//=======================================================================//

void EnvelopFunctionApprox::output_eigen_function(unsigned int state_number,  const std::string& filename)
{
  //===========================================

  DofMap& dof_map = system->get_dof_map();

  const Mesh& mesh1 = system->get_mesh();
 
  MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
  const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();
  

  unsigned int number_of_points = 0;

  for ( ; nd != nd_el ; ++nd)  number_of_points++;
    


  //vector of names
  vector<string> output_names(psi_name.size() * 2);
 

  unsigned int number_output_data = output_names.size() ;

  for (short i = 0; i < number_output_data/2; i++)
    {
      output_names[i*2    ] = psi_name[i] + "_r";
      output_names[i*2 + 1] = psi_name[i] + "_i";
    }
  

  
 
  //vector 
  unsigned int  point_index = 0;
  
  unsigned int output_size = mesh->n_nodes()  * opt.number_of_bands * 2;

  vector<double>  psi_data(output_size);
 
 
  MeshBase::const_element_iterator it = mesh1.active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh1.active_local_elements_end();

  std::vector<unsigned int> dof_indices;

 

  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;
      

      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  for (unsigned int n = 0; n < elem->n_nodes(); n++)
	    { 
	      unsigned int  node_id =  elem->node(n);
	      Complex value =  (solution[state_number].eigen_vector[ dof_indices[n] ]);
	      
	      psi_data[psi_index*2 + node_id*number_output_data] = value.real();
	      
	      psi_data[psi_index*2 + 1  + node_id*number_output_data] = value.imag();
	      
	      
	      
	    }
	}
    }


    


  //std :: cout << filename << "\n";


  if (opt.output_type == "GMV")     
    GMVIO(mesh1).write_nodal_data(filename, psi_data, output_names);

  if (opt.output_type == "tecplot") 
    TecplotIO(mesh1,false).write_nodal_data(filename,psi_data,output_names);

  if (opt.output_type == "gnuplot" || opt.output_type == "GNUplot") 
    GnuPlotIO(mesh1).write_nodal_data(filename, psi_data, output_names);
  

}

//=======================================================================//

void  EnvelopFunctionApprox::prepare_probability_function(const unsigned int state_number, std::vector<double>& prob_data)
{

  DofMap& dof_map = system->get_dof_map();

  const Mesh& mesh1 = system->get_mesh();

  MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
  const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();

  unsigned int number_of_points = 0;

  for ( ; nd != nd_el ; ++nd)  number_of_points++;

 
  
  prob_data.resize(number_of_points, 0.0);

  MeshBase::const_element_iterator it = mesh1.active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh1.active_local_elements_end();

  std::vector<unsigned int> dof_indices;


  //!calculation of probability function
  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;
      

      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  for (unsigned int n = 0; n < elem->n_nodes(); n++)
	    { 
	      unsigned int  node_id =  elem->node(n);
	      Complex value =  (solution[state_number].eigen_vector[ dof_indices[n] ]);
	     
	      prob_data[node_id] += std::abs(value) * std::abs(value);
	      
	      
	      
	    }
	}
    }
  //done

 
}


void EnvelopFunctionApprox::output_probability_function(unsigned int state_number,  const std::string& filename)
{
  //===========================================

  DofMap& dof_map = system->get_dof_map();

  const Mesh& mesh1 = system->get_mesh();

  MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
  const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();

  unsigned int number_of_points = 0;

  for ( ; nd != nd_el ; ++nd)  number_of_points++;
    


  //vector of names
  vector<string> output_names(1);
  output_names[0] = "|psi|^2";

  
  unsigned int  point_index = 0;
  
  

  vector< vector<Complex> >  psi_data;
  psi_data.resize(number_of_points);
  for (unsigned int i = 0; i < number_of_points; i++) psi_data[i].resize( opt.number_of_bands, Complex(0.0, 0.0) ); 
 
  vector<double>  probability_data(number_of_points, 0.0);
 
  MeshBase::const_element_iterator it = mesh1.active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh1.active_local_elements_end();

  std::vector<unsigned int> dof_indices;


  //!calculation of psi

  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;
      

      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  for (unsigned int n = 0; n < elem->n_nodes(); n++)
	    { 
	      unsigned int  node_id =  elem->node(n);
	      Complex value =  (solution[state_number].eigen_vector[ dof_indices[n] ]);
	      
	      psi_data[node_id][psi_index] = value;
	      
	    
	      
	      
	      
	    }
	}
    }


  //calculation of |psi|^2
  double t1;
  for (unsigned int i = 0; i < number_of_points; i++)
    for (unsigned int j = 0; j < opt.number_of_bands; j++)
      {
	t1 = std::abs(psi_data[i][j]);
	probability_data[i] += t1 * t1; 
      }


  //std :: cout << filename << "\n";


  if (opt.output_type == "GMV")     
    GMVIO(mesh1).write_nodal_data(filename, probability_data, output_names);

  if (opt.output_type == "tecplot") 
    TecplotIO(mesh1,false).write_nodal_data(filename, probability_data,output_names);

  if (opt.output_type == "gnuplot" || opt.output_type == "GNUplot") 
    GnuPlotIO(mesh1).write_nodal_data(filename, probability_data,output_names);

}



//=======================================================================//


//=======================================================================//



void EnvelopFunctionApprox::apply_diriclet_bc_at_all_boundaries()
{
  MeshBase::const_element_iterator it = mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh->active_local_elements_end();

  dirichlet_dofs.clear();
  DofMap& dof_map = system->get_dof_map();
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
	    side_is_external = true;
	  else
	    {
	      std::vector< const Elem * > active_family;
	      if ( el1->has_children() )
		{
		  el1->active_family_tree (active_family);
	 
		  if (active_family.size() == 0)
		    side_is_external = true;
		  //TODO
		  /*
		    has to be corrected because it may contain active child that does not belong to boundary
		  */
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
			for (short band = 0 ; band <  opt.number_of_bands; band++)
			  {
			    dof_map.dof_indices (elem, dof_indices,band); 
			    dirichlet_dofs.insert(dof_indices[nd]);
			  }
		    }
		      
		}
	      else
		{//1D
		  for (short band = 0 ; band <  opt.number_of_bands; band++)
		    {
		      dof_map.dof_indices (elem, dof_indices,band);
		      dirichlet_dofs.insert(dof_indices[i]);
		    }
		}
	    }
	      
	      
	      
	      
	}
	   
	   
    }
    


}
//======================================================================//
void  EnvelopFunctionApprox::create_dirichlet_dofs( )
{
  
  SimulationEnvironment& se = get_environment(); 

  DofMap& dof_map = system->get_dof_map();

  MeshBase::const_element_iterator it = mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh->active_local_elements_end();

  dirichlet_dofs.clear();


 
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

	  if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) )
	    if ( bd->get_name() == "Dirichlet" || bd->get_name() == "dirichlet" )
	    {
	      for (short band = 0 ; band <  opt.number_of_bands; band++)
	      {
		dof_map.dof_indices (elem, dof_indices,band); 
		dirichlet_dofs.insert(dof_indices[n]);
	      }
	     
	    }

	  
	}
      
    }
  

}

//=======================================================================//


//=======================================================================//
void EnvelopFunctionApprox::make_constraints(void)
{
 
  DofMap& dof_map = system->get_dof_map();
  
 
  //----------------------------------------------------------------------------//
  //I recalculate my copy of the dof constraints because I need them explicitely!
  my_dof_constraints.clear();

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
//=======================================================================//

void EnvelopFunctionApprox::make_new_dofs( )
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
//=======================================================================//
bool EnvelopFunctionApprox::compare_eigen_energy_electrons(eigen_propblem_solution state1, eigen_propblem_solution state2)
{
  return(state1.eigen_energy < state2.eigen_energy );
}


//=======================================================================//
bool EnvelopFunctionApprox::compare_eigen_energy_holes(eigen_propblem_solution state1, eigen_propblem_solution state2)
{
  return(state1.eigen_energy > state2.eigen_energy );
}

//=======================================================================//

bool EnvelopFunctionApprox::compare_eigen_energy_holes1(eigen_energy state1, eigen_energy state2)
{
  return(state1.energy> state2.energy );
}

//=======================================================================//
bool EnvelopFunctionApprox::compare_eigen_energy_electrons1(eigen_energy state1, eigen_energy state2)
{
  return(state1.energy<  state2.energy );
}


//=======================================================================//
EnvelopFunctionApprox:: ~EnvelopFunctionApprox(void)
{
  es->delete_system(system_name);
}

//=======================================================================//
void EnvelopFunctionApprox::apply_periodic_bc()
{

  

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Periodic bc. Assembly",false);


  
 
  // The dimension that we are running
  //const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving


  vector<unsigned int> uvar(opt.number_of_bands);

  for (unsigned int i = 0; i < opt.number_of_bands; i++) 
    {
      uvar[i] = system->variable_number(psi_name[i]);
    }

  unsigned int system_number=system->number();
  
  DofMap& dof_map = system->get_dof_map();
  
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
      if  (opt.periodicity[i]) //Check if the periodic b.c. are applied along the direction i
	
	{
	 
	  std :: vector <const Node*>& vec =  nodes_periodic[i];

	  for (unsigned int n = 0; n < vec.size(); n++) // Loop over all the nodes
	    {
	      const Node* node1 = vec[n];
	      
	    
		 
	      for (unsigned int var_index = 0 ; var_index  < opt.number_of_bands ;  var_index ++)
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
		      MeshBase::const_element_iterator el3  = mesh->level_elements_begin(refinement_level);
		      MeshBase::const_element_iterator end_el3 = mesh->level_elements_end(refinement_level);
		      
		      Elem*  elem1;
		      for ( ; ( (el3 != end_el3) ) ; ++el3)  
			{
			  Elem* elem = *el3;
			  if (element_on_boundary(elem))
			    {
			      if (elem->contains_point(point2))
				{
				  elem1 = elem;
				  
				  break;
				  
				}
			    }
			}
		      
		      //children of the  most coarse element 
		      while ( !( elem1->active() ) )
			{
			  
			  for (unsigned int i=0 ; i < elem1->n_children() ; i++)
			    {
			      Elem* 	child = elem1->child(i);
			      if (element_on_boundary(child))
				{
				  if (child->contains_point(point2))
				    {
				      elem1 = child;
				      break;
				    }
				}
			    }
			}
			
		      

		      
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

void EnvelopFunctionApprox::make_nodes_periodic()
{
  const double pos_tol = 1e-10;
 
  nodes_periodic.clear();

 
  

  for (unsigned dir = 0; dir <=dim-1; dir++)
    {//directions
      std::vector < const Node*> temp_vec;
      temp_vec.clear();
    
      if (opt.periodicity[dir]) 
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


//-----------------------------------------------------------------------------//
double  EnvelopFunctionApprox::eigenstate_norm(unsigned int state_number)
{
  double  result;
  
  const vector< Complex > &  eigen_vector =  solution[state_number].eigen_vector;
  

  DofMap& dof_map = system->get_dof_map();

  FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  QGauss qrule (dim, SECOND);

  fe -> attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  std::vector<unsigned int> dof_indices;


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


  //My Jacobian 

  double length_scale = opt.length_scale;

  my_Jacobian = 1.0;
  for (short i = 1; i <= dim; i++)
    my_Jacobian *= length_scale;


  Complex temp(0.0, 0.0);
  Complex eigen_f_value1, eigen_f_value2;

  for ( ; el != end_el ; ++el) 
    {//el
      
      const Elem* elem = *el;
      fe->reinit (elem);

      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
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



  result = sqrt( abs(temp) * my_Jacobian );


  return(result);

} 



//==========================================================//


double EnvelopFunctionApprox::calculate_fermi_averaged(unsigned int i)
{

  Complex  result(0.0,0.0);


  //-----------------------------------------------------//

 



  const vector< Complex >&   eigen_vector =  solution[i].eigen_vector;
  

 
 
  //----------------------------------------------------//

  
  const Mesh* mesh = &(es->get_mesh());


  unsigned int dim = mesh->mesh_dimension();
  


  system = &( es->get_system<LinearImplicitSystem>(system_name));

  DofMap& dof_map = system->get_dof_map();
  


  

  //My Jacobian 

  double length_scale = opt.length_scale;

  my_Jacobian = 1.0;
  for (short i = 1; i <= dim; i++)
    my_Jacobian *= length_scale;

  
   system->init();
   

   FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

   AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

   // A 5th order Gauss quadrature rule for numerical integration.
   QGauss qrule (dim, SECOND);

   // Tell the finite element object to use our quadrature rule.
   fe -> attach_quadrature_rule (&qrule);

   // The element Jacobian * quadrature weight at each integration point.   
   const std::vector<Real>& JxW = fe->get_JxW();

   // properties at the quadrature points.
   const std::vector<Point>& q_point = fe->get_xyz();
   
   // The element shape functions evaluated at the quadrature points.
   const std::vector<std::vector<Real> >& phi = fe->get_phi();
  

   //------------------------------------------------------------
   std::vector<unsigned int> dof_indices_component;
 
   std::vector<unsigned int> dof_indices;

   //-------------------------------------------------------------
   
  //----------------------------------------------------//


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  
  Complex eigen_f_value1;
  Complex eigen_f_value2;


  DriftDiffusion::Solution dd_solution;
  double chem_pot_value_eV;
  
  for ( ; el != end_el ; ++el) 
    {//el

      const Elem* elem = *el;
      fe->reinit (elem);


      Point center = elem->centroid();

      poisson_equation->get_solution(elem, center, dd_solution);

      if (opt.particle == "el")
	chem_pot_value_eV = dd_solution.fermi_e;
      else
	chem_pot_value_eV = dd_solution.fermi_h;

     
      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
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
		      result += ( JxW[qp] * phi[p1][qp] * eigen_f_value1 *  phi[p2][qp] * conj(eigen_f_value2) ) * chem_pot_value_eV;
		    
		    }
		}

	    }


	
	  
	}
      

    }

 


  
  result *= my_Jacobian;


 
  return(result.real());
}

//---------------------------------------------------------------------------//

void EnvelopFunctionApprox::calculate_density(double Temperature)
{
  DofMap& dof_map = system->get_dof_map();

  const Mesh& mesh1 = system->get_mesh();

  
  _density.clear();



  MeshBase::const_element_iterator       nd     = mesh1.active_elements_begin();
  const MeshBase::const_element_iterator nd_el  = mesh1.active_elements_end();

 

  vector<double> density_of_state;

  unsigned int number_of_eigenfunctions = solution.size();

  for (unsigned int i = 0; i < number_of_eigenfunctions; i++)
    {
      
      

      const double Fermi_energy = solution[i].Fermi_energy;

      const double Energy = solution[i].eigen_energy;
      
      double prob_factor = Fermi_statistics_probability(Energy, Fermi_energy, Temperature); //Thermal probability

     
      density_of_state = calculate_cell_prob_function(i);

      MeshBase::const_element_iterator       it     = mesh1.active_elements_begin();
      const MeshBase::const_element_iterator it_end  = mesh1.active_elements_end();
     
      unsigned int el_number = 0;
      
      for( ;it !=  it_end ;++it)
      {
	const Elem* el = *it;
	
	double temp = density_of_state[el_number] * prob_factor;
	
	if (it == it_end)
	  _density.insert(pair<const Elem*, double> (el, temp));
	else
	  _density[el] += temp;
	
	el_number++;
      }
      
     

      
    }


  
  
}

//===========================================================//
vector<double>  EnvelopFunctionApprox::calculate_cell_prob_function(unsigned int state_number)
{

  const vector< Complex > &  eigen_vector =  solution[state_number].eigen_vector;

  DofMap& dof_map = system->get_dof_map();

  const Mesh& mesh = system->get_mesh();

  FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  QGauss qrule (dim, SECOND);

  fe -> attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  std::vector<unsigned int> dof_indices;
  

  MeshBase::const_element_iterator       el1     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el1 = mesh.active_elements_end();

  //--------------------------------------------------------------//
  //define results
  unsigned int active_el_number = 0;

  for ( ; el1 != end_el1; ++el1) active_el_number++;

  vector<double> result(active_el_number, 0);
  //--------------------------------------------------------------//


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();


 
  Complex eigen_f_value1, eigen_f_value2;
  unsigned int el_number = 0;

  for ( ; el != end_el ; ++el) 
    {//el
      
      double el_volume = 0;

      const Elem* elem = *el;
      fe->reinit (elem);

      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
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

		      Complex temp =  eigen_f_value1 * conj(eigen_f_value2);

		      result[el_number] +=  JxW[qp] * phi[p1][qp] *  phi[p2][qp] * temp.real() ;

		    }
		}
	      
	    }


	  
	  
	}
      
      //-------------------------------------------------//
      //has to be removed for the new libmesh
	  
      for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
	el_volume += JxW[qp];
      
      //untill here
      //-------------------------------------------------//

      /*
	for the new libmesh
	el_volume = elem->volume
      */

      


      result[el_number] /= el_volume;
      el_number++;

    }
  

  

 


  return(result);

}




//===========================================================//



vector<double> EnvelopFunctionApprox::calculate_prob_function(unsigned int state_number)
{
//===========================================

  

  DofMap& dof_map = system->get_dof_map();

  const Mesh& mesh1 = system->get_mesh();

  MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
  const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();

  unsigned int number_of_points = 0;

  for ( ; nd != nd_el ; ++nd)  number_of_points++;
    
  unsigned int  point_index = 0;
  

  vector< vector<Complex> >  psi_data;
  psi_data.resize(number_of_points);
  for (unsigned int i = 0; i < number_of_points; i++) psi_data[i].resize( opt.number_of_bands, Complex(0.0, 0.0) ); 
 
  vector<double>  probability_data(number_of_points, 0.0);
 
  MeshBase::const_element_iterator it = mesh1.active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh1.active_local_elements_end();

  std::vector<unsigned int> dof_indices;


  //!calculation of psi

  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;
      

      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  for (unsigned int n = 0; n < elem->n_nodes(); n++)
	    { 
	      unsigned int  node_id =  elem->node(n);
	      Complex value =  (solution[state_number].eigen_vector[ dof_indices[n] ]);
	      
	      psi_data[node_id][psi_index] = value;
	      
	    
	      
	      
	      
	    }
	}
    }


  //calculation of |psi|^2
  double t1;
  for (unsigned int i = 0; i < number_of_points; i++)
    for (unsigned int j = 0; j < opt.number_of_bands; j++)
      {
	t1 = std::abs(psi_data[i][j]);
	probability_data[i] += t1 * t1; 
      }


  
  //------

  return(probability_data);

  //-----

}


//=================================================================//

double EnvelopFunctionApprox::get_integrated_probability()
{
  double result = 0;
  unsigned int number_of_eigs = solution.size();
  for (unsigned int i = 0 ; i <  number_of_eigs; i++)
    { 
    
      result += Fermi_statistics_probability(solution[i].eigen_energy, solution[i].Fermi_energy,opt.Temperature);
    }
 
  return(result);
}


//=================================================================//
void EnvelopFunctionApprox::calculate_convergent_density(double T)
{

 

  unsigned int number_of_states = opt.initial_eigenstates_number;

  solve_eigen_value_problem(number_of_states);

  unsigned int n1 = solution.size();

  double last_state_density = Fermi_statistics_probability(solution[n1-1].eigen_energy, 
							   solution[n1-1].Fermi_energy,
							   T);

  double total_density =  get_integrated_probability();
 
  bool converged =( last_state_density/total_density  < opt.relative_density_tolerance ) ; 

  if  (opt.convergent_density && (!converged))
    {
      number_of_states = (unsigned int) (number_of_states * opt.eigen_number_increase_factor) + 1;
      
      opt.solve_ev_problem_twice = false;

      double st_shift_value = (solution[0].eigen_energy - opt.spectrum_shift)/Hartree;

      solve_eigen_value_problem(number_of_states, st_shift_value);

      double total_density1 = get_integrated_probability();

      if ( abs(total_density1 - total_density)/total_density < opt.relative_density_tolerance )  converged = true;

      total_density = total_density1;
    }

  calculate_density(T );
 
}

//==============================================================================//  

unsigned int EnvelopFunctionApprox::get_number_of_active_cells()
{
  unsigned int result = 0;
  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  for ( ; el != end_el ; ++el)
    result++;


  return(result);


  
}

//===============================================================================//
short EnvelopFunctionApprox::calculate_number_of_bands(void) const
{
  
  short result = 0;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  EFAbulkHamiltonian* element_hamiltonian;
  
  const Elem* elem = *el;
  const ID subdomain = elem->subdomain_id();
  const Material* mat = _device->get_material(subdomain);


  element_hamiltonian = (  dynamic_cast<EFAbulkModel*> ( mat ->get_model(get_id()) )    ) -> get_Hamiltonian_model() ; 

  std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&  
	    model_Ham = ( element_hamiltonian->get_Hamiltonian() );

  result = model_Ham.size();

  if (result == 0)  throw InitFailedException("EnvelopFunctionApprox: Hamiltonian is empty");

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    const ID subdomain = elem->subdomain_id();
    const Material* mat = _device->get_material(subdomain);


    element_hamiltonian = ( dynamic_cast<EFAbulkModel*> (  mat ->get_model(get_id()) ) )->get_Hamiltonian_model() ;
 
    const std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&  
      model_Ham = ( element_hamiltonian->get_Hamiltonian() );


    short result1 = model_Ham.size();

    if (result1 != result) throw InitFailedException("EnvelopFunctionApprox: Hamiltonians of different materials have different number of bands");


  }

  return(result);
 
}
//========================================================================================//
