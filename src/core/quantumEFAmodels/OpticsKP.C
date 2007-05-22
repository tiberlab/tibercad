
#include "OpticsKP.h"
#include <stdlib.h>

#include "SimulationInterface.h"
#include "KPbulkHamiltonian.h"
#include "EnvelopFunctionApprox.h"
#include "EFAbulkModel.h"
#include "SimulationEnvironment.h"
#include "Material.h"
using namespace std; 
using namespace Constants;

Device*   OpticsKP:: _device;

OpticsKP::~OpticsKP()
{
  // es->delete_system(system_name);
}
//===============================================//
PhysicalModel* OpticsKP::create_physical_model(const ModelOptions& options) const throw (ModelErrorException)
{
  ModelOptions kp8x8options = options;

  kp8x8options["model"] = "kp";
  kp8x8options["kp_model"] = "8x8";

  EFAbulkModel* model = dynamic_cast<EFAbulkModel*> ( PhysicalModelInterface::create("EFAmodel", kp8x8options ) );

  if (model == NULL)
    throw ModelErrorException("OpticsKP: cannot create EFAbulkModel");
 
  return(model);
  
}

//===============================================//
BoundaryProperties* OpticsKP::create_boundary_model(const ModelOptions& options) const  throw (ModelErrorException)
{
  
  return NULL;
}

//===============================================//
OpticsKP::OpticsKP()

{

  initial_state_model = NULL;
  final_state_model = NULL;


}

//==============================================//

void OpticsKP::parse_options()
{
  const ModelOptions& mod_opt = get_options();
  
  _initial_eigen_state_numbers.clear();
  const std::vector<EnvelopFunctionApprox::eigen_propblem_solution>& in_solution = initial_state_model->get_solution();

  {
    std::vector<unsigned int> temp;
    mod_opt.get_option("initial_eigenstates", temp);
    if (temp.size() == 2)
      if (temp[0] <= temp[1])
	if (temp[0] >= 0 && temp.size() <=in_solution.size())
	{
	  _initial_eigen_state_numbers.resize(temp[1] - temp[0] + 1 ); 
	  unsigned int j = 0;
	  for (unsigned i = temp[0]; i <= temp[1]; i++)
	  {
	    _initial_eigen_state_numbers[j] = i;
	    j++;
	  }
	}
	
    if (_initial_eigen_state_numbers.size() == 0)
      throw InitFailedException("OpticsKP: initial states are not acceptable\n");
  }

  _final_eigen_state_numbers.clear();
  const std::vector<EnvelopFunctionApprox::eigen_propblem_solution>& fin_solution = final_state_model->get_solution();

  {
    std::vector<unsigned int> temp;
    mod_opt.get_option("final_eigenstates", temp);
    if (temp.size() == 2)
      if (temp[0] <= temp[1])
	if (temp[0] >= 0 && temp.size() <=fin_solution.size())
	{
	  _final_eigen_state_numbers.resize(temp[1] - temp[0] + 1 ); 
	  unsigned int j = 0;
	  for (unsigned i = temp[0]; i <= temp[1]; i++)
	  {
	    _final_eigen_state_numbers[j] = i;
	    j++;
	  }
	}
	
    if (_final_eigen_state_numbers.size() == 0)
      throw InitFailedException("OpticsKP: final states are not acceptable\n");
  }


}

//==============================================//
void OpticsKP::do_init()
{

  SimulationEnvironment& si = get_environment();   

  _device = &( si.get_device() );

  

  //initial state----------------
  const ModelOptions& mod_opt = get_options();
  if  (mod_opt.find_option("initial_state_model"))
  {
    std::string quantum_model;
    quantum_model = mod_opt.get_option("initial_state_model" , "");
    initial_state_model = dynamic_cast<EnvelopFunctionApprox*> ( find_simulation ( quantum_model));
    if (initial_state_model == NULL)
      throw InitFailedException("OpticsKP: initial_state_model " + quantum_model + " does not exist\n");
  }
  else
    throw InitFailedException("OpticsKP: initial_state_model must be defined\n");
  //--------------------------

  //final state----------------
  if  (mod_opt.find_option("final_state_model"))
  {
    std::string quantum_model;
    quantum_model = mod_opt.get_option("final_state_model" , "");
    final_state_model = dynamic_cast<EnvelopFunctionApprox*> ( find_simulation ( quantum_model));
    if (final_state_model == NULL)
      throw InitFailedException("OpticsKP: final_state_model " + quantum_model + " does not exist\n");
  }
  else
    throw InitFailedException("OpticsKP: final_state_model must be defined\n");
  //---------------------------
  system_name = get_equation_system_name();
  es = &(get_equation_systems());

  //-------------------------------
 //---------------------------------------------------------------------------------------------------------//
  //My Jacobian 
  double mesh_units = get_environment().get_device().get_mesh_units();
  length_scale = mesh_units/bohr_radius;

  unsigned int dim = (es->get_mesh()).mesh_dimension();

  my_Jacobian = 1.0;
  for (short i = 1; i <= dim; i++)
    my_Jacobian *= length_scale;

   //--------------------------------------------------------------------------------------------------------//
  
}

//==============================================//
void OpticsKP::do_solve()
{

  parse_options();

  calculate_matrix();

  //temporary_solution------------------------------
  unsigned int n1 =  _initial_eigen_state_numbers.size();
  unsigned int n2 =  _final_eigen_state_numbers.size();
  std::vector< std::vector <std::vector <  Complex  >  >  >   P_matrix;
  get_P_matrix_elements(P_matrix);

 

  std::ostringstream os3;
  os3 << "output/optics.out";
	  
  std::ofstream out_optics(  (os3.str()).c_str()   );
	  
  for (int i = 0; i < n1; i++) 
    for (int j = 0; j < n2; j++) 
    {
      out_optics << i << "   " << j <<  "  " ;
      for (int p = 0; p < 3; p++) 
      {
	cerr << p << "  " << i  <<"   " << j   <<  P_matrix[p][i][j] << "\n";

	out_optics << P_matrix[p][i][j].real() << "   " <<  P_matrix[p][i][j].imag() <<"         " ;
      }
      out_optics << "\n";
    }

  //------------------------------------------------
}


//=====================================================================================================
void OpticsKP::get_P_matrix_elements (std::vector< std::vector <std::vector <  Complex  >  >  > &  P_matrix) 
{
  unsigned int n_i =  _initial_eigen_state_numbers.size();
  unsigned int n_f =  _final_eigen_state_numbers.size();


  P_matrix.clear();
  P_matrix.resize(3);
  for (unsigned i = 0; i < 3; i++)
    {
      P_matrix[i].resize(n_i);
      for (unsigned j = 0; j < n_i; j++)   P_matrix[i][j].resize(n_f);
    } 

  for (unsigned int i1 = 0; i1 < n_i; i1++)
    for (unsigned int i2 = 0; i2 < n_f; i2++)
      {
	vector<Complex> matr_elements =  calculate_matrix_element(i1, i2);
	for (unsigned i = 0; i < 3; i++)  P_matrix[i][i1][i2] = matr_elements[i];
      }



}

//=========================================================================//
void OpticsKP::calculate_matrix(void)
{
  
  
  const Mesh* mesh = &(es->get_mesh());


  unsigned int dim = mesh->mesh_dimension();

  es->add_system<LinearImplicitSystem> (system_name);

  system = &( es->get_system<LinearImplicitSystem>(system_name));


  //---------------------------------------------------------------------------------------
  //add variables

  vector<string> psi_name;

  psi_name.clear();

  for (short i = 0; i < 8; i++)
  {
    std::ostringstream var_str;
    var_str << "psi" << i ;
    string name = var_str.str();
    psi_name.push_back(name);

    system->add_variable(name,FIRST);
  }



  vector<unsigned int> psivar(8);
  for (unsigned int i = 0; i < 8; i++)
  {
    psivar[i] = system->variable_number(psi_name[i]);
  }


  //add matrixes

  


  DofMap& dof_map = system->get_dof_map();



  system->add_matrix("Px_real"); 
  Px_matr_real = & (system->get_matrix("Px_real"));

  system->add_matrix("Py_real"); 
  Py_matr_real = & (system->get_matrix("Py_real"));

  system->add_matrix("Pz_real"); 
  Pz_matr_real = & (system->get_matrix("Pz_real"));

  system->add_matrix("Px_imag"); 
  Px_matr_imag = & (system->get_matrix("Px_imag"));

  system->add_matrix("Py_imag"); 
  Py_matr_imag = & (system->get_matrix("Py_imag"));

  system->add_matrix("Pz_imag"); 
  Pz_matr_imag = & (system->get_matrix("Pz_imag"));


 

  system->init();




  FEType fe_type = dof_map.variable_type(psivar[0]); //all the variable have the same FE representation

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe -> attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe->get_JxW();

  // properties at the quadrature points.
  const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  

  // The element shape function gradients evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();


  //------------------------------------------------------------
  std::vector<unsigned int> dof_indices_component;
 
  std::vector<unsigned int> dof_indices;

  //-------------------------------------------------------------

  DenseMatrix<Number> Px_real; 
  DenseMatrix<Number> Px_imag;

  DenseMatrix<Number> Py_real; 
  DenseMatrix<Number> Py_imag;

  DenseMatrix<Number> Pz_real; 
  DenseMatrix<Number> Pz_imag;
 

  DenseSubMatrix<Number> Px_real_sub(Px_real);
  DenseSubMatrix<Number> Px_imag_sub(Px_imag);

  DenseSubMatrix<Number> Py_real_sub(Py_real);
  DenseSubMatrix<Number> Py_imag_sub(Py_imag);

  DenseSubMatrix<Number> Pz_real_sub(Pz_real);
  DenseSubMatrix<Number> Pz_imag_sub(Pz_imag);


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();
  //                                                              //
  // we do not apply neither electric potential nor strain        //
  // because in our model they do not affect optical properties   //
  //                                                              //

  EFAbulkHamiltonian* element_hamiltonian;
  KPbulkHamiltonian* element_kp_hamiltonian;

  unsigned int el_number = 0;

  double temp;
  for ( ; el != end_el ; ++el) 
    {//el
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      const ID subdomain = elem->subdomain_id();
      const Material* mat = _device->get_material(subdomain);

      element_hamiltonian = (  dynamic_cast<EFAbulkModel*> (  mat ->get_model(get_id()) )  )->get_Hamiltonian_model(); 

      element_kp_hamiltonian = dynamic_cast<KPbulkHamiltonian*>  (element_hamiltonian);
   
      

      dof_map.dof_indices (elem, dof_indices); 
      const unsigned int n_dofs   = dof_indices.size();
      fe->reinit (elem);

      Px_real.resize(n_dofs, n_dofs );
      Px_imag.resize(n_dofs, n_dofs );

      Py_real.resize(n_dofs, n_dofs );
      Py_imag.resize(n_dofs, n_dofs );

      Pz_real.resize(n_dofs, n_dofs );
      Pz_imag.resize(n_dofs, n_dofs );

 

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{//qp
	  const vector < vector <vector <EFAbulkHamiltonian::MatrixElement> > >&  
	    P = element_kp_hamiltonian->get_optical_operator() ;
	  


	  for (unsigned int band1 = 0; band1 < 8; band1++)
	    {//band1
	      dof_map.dof_indices (elem, dof_indices_component, psivar[band1]);
	      
	      const unsigned int n_psi_dofs = dof_indices_component.size();
	      for (unsigned int band2 = 0; band2 < 8; band2++)
		 {//band2
		   Px_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
		   Px_imag_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);

		   Py_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
		   Py_imag_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);

		   Pz_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
		   Pz_imag_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);

		   //Integration----------------------------------------------------------------
		   for (unsigned int p1=0; p1<n_psi_dofs; p1++)
		     {//p1
		       for (unsigned int p2=0; p2<n_psi_dofs; p2++)
			 {//p2
			  
			  
			   vector<Complex> value(3, Complex(0.0,0.0));

			   //----constant part-----------------------------------
			   temp = JxW[qp] * phi[p1][qp] * phi[p2][qp];

			   for (short pol = 0; pol < 3; pol++)
			     value[pol] += temp * P[pol][band1][band2].constant ;

			 

			   //----            ------------------------------------
			  
			   //--------linear left-----------------------------------------
			   for (short i = 0; i < dim; i++)
			     {
			       temp = JxW[qp]* dphi[p1][qp](i) * phi[p2][qp];
			       for (short pol = 0; pol < 3; pol++)
				 value[pol] -= temp * P[pol][band1][band2].linear_left[i]*Complex(0.0, -1.0) /length_scale; 
			     } 
			   //--------           ----------------------------------------


			   //--------linear right---------------------------------------
			   for (short i = 0; i < dim; i++)
			     {
			       temp = JxW[qp]* dphi[p2][qp](i) * phi[p1][qp];
			       for (short pol = 0; pol < 3; pol++)
				value[pol] += temp * P[pol][band1][band2].linear_right[i] * Complex(0.0, -1.0) /length_scale;
			     }
                           //-----------------------------------------------------------
			     
			   for (short pol = 0; pol < 3; pol++) value[pol] *= my_Jacobian;

			   Px_real_sub(p1,p2) += value[0].real();
			   Py_real_sub(p1,p2) += value[1].real();
			   Pz_real_sub(p1,p2) += value[2].real();

			   Px_imag_sub(p1,p2) += value[0].imag();
			   Py_imag_sub(p1,p2) += value[1].imag();
			   Pz_imag_sub(p1,p2) += value[2].imag();
			   

			 }
		  

		     }


		 }
	      
	    }

	}
      Px_matr_real->add_matrix(Px_real,dof_indices );
      Py_matr_real->add_matrix(Py_real,dof_indices );
      Pz_matr_real->add_matrix(Pz_real,dof_indices );
      
      Px_matr_imag->add_matrix(Px_imag,dof_indices );
      Py_matr_imag->add_matrix(Py_imag,dof_indices );
      Pz_matr_imag->add_matrix(Pz_imag,dof_indices );
	  
      
      el_number++;
    
    }
  
}


//----------------------------------------------------------------------------------------//

std::vector<Complex> OpticsKP::calculate_matrix_element(unsigned int i, unsigned int j)
{
  vector<Complex> result(3,Complex(0.0, 0.0));

  vector<  PetscMatrix<Number>* > P_real_p(3);
  vector<  PetscMatrix<Number>* > P_imag_p(3);
  
  const Mesh& mesh = system->get_mesh();

  
 
  
  
  DofMap& dof_map = system->get_dof_map();
  
  unsigned int number_of_nodes  = dof_map.n_dofs()/dof_map.n_variables(); //only nodes that belong to active elements

 

  P_real_p[0]  = static_cast< PetscMatrix<Number>* > (Px_matr_real); P_real_p[0]->close();  
  P_real_p[1]  = static_cast< PetscMatrix<Number>* > (Py_matr_real); P_real_p[1]->close();
  P_real_p[2]  = static_cast< PetscMatrix<Number>* > (Pz_matr_real); P_real_p[2]->close();
  P_imag_p[0]  = static_cast< PetscMatrix<Number>* > (Px_matr_imag); P_imag_p[0]->close();
  P_imag_p[1]  = static_cast< PetscMatrix<Number>* > (Py_matr_imag); P_imag_p[1]->close();
  P_imag_p[2]  = static_cast< PetscMatrix<Number>* > (Pz_matr_imag); P_imag_p[2]->close();

  const vector< Complex >   eigen_vector_i =  (initial_state_model->get_solution())[i].eigen_vector;
  const vector< Complex >   eigen_vector_f =  (final_state_model->get_solution())[j].eigen_vector;


 
  int size_matrix = Px_matr_real->n();



  //!number of bands in initial state
  short    num_bands_initial = initial_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_in = initial_state_model->get_kp_bands();

  assert( kp_bands_map_in.size() > 0);

  //number of bands in final state
  short    num_bands_final   = final_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_fi = final_state_model->get_kp_bands();

  assert( kp_bands_map_fi.size() > 0);

  map<short, short>::const_iterator  band_it;
  
 

  for (short pol = 0; pol < 3; pol++)
    {//polarization
      for (int row = 0 ; row < size_matrix; row++)
	{//rows of P matrix
	  
	 // short band_number1 = row%8;
	  short band_number1 = row/number_of_nodes;
	 
	  band_it = kp_bands_map_in.find( band_number1 );
	  

	  if (band_it != kp_bands_map_in.end())
	    {//band exists in kp model of the initial state
	     
	      unsigned int dof_in_initial_eigenvector = band_it->second * number_of_nodes + row%number_of_nodes;
	      int ierr = 0;
	      const  PetscScalar *petsc_row_vals_real;
	      const  PetscInt *petsc_cols_real;
	      int n_cols_real = 0;

	      const  PetscScalar *petsc_row_vals_imag;
	      const  PetscInt *petsc_cols_imag;
	      int n_cols_imag = 0;
	    

	      map<int,double> real_values, imag_values;
	      map<int, double>::iterator  position;
	 

	      set<int> real_column, imag_column, complex_column;
	      set<int>::iterator com_col_it;

	      insert_iterator<set<int> >  com_ins(complex_column,complex_column.begin() );
	      
	      
	      ierr = MatGetRow(P_real_p[pol]->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	      real_column.clear();
	    
	      for (int i = 0; i < n_cols_real; i++) 
		if (petsc_row_vals_real[i] != 0.0)  
		  {
		    real_column.insert(petsc_cols_real[i]);
		 
		    real_values.insert(make_pair(petsc_cols_real[i],petsc_row_vals_real[i] ));
		  }
	    
	      
	      ierr = MatGetRow(P_imag_p[pol]->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	      imag_column.clear();

	     
	      for (int i = 0; i < n_cols_real; i++) 
		if (petsc_row_vals_imag[i] != 0.0)   
		  {	    
		    imag_column.insert(petsc_cols_imag[i]);
		    imag_values.insert(make_pair(petsc_cols_imag[i],petsc_row_vals_imag[i] ));
		  }
	    

	      set_union(real_column.begin(), real_column.end(), imag_column.begin(), imag_column.end(), com_ins);
	      
	      for (com_col_it = complex_column.begin(); com_col_it != complex_column.end(); com_col_it++)
		{
		  int n1 = *com_col_it;
		
		  short band_number2  = n1/number_of_nodes;

		  band_it = kp_bands_map_fi.find(band_number2 );

		  if (band_it != kp_bands_map_fi.end())
		    {//band exists in kp model of the final state
		      double value_real;
		      double value_imag;
	    
		      //real part------	  
		      position = real_values.find(n1);
		      if (position != real_values.end()) 
			value_real = position->second;
		      else 
			value_real = 0.0;
		        
		      //imag part 
		      position = imag_values.find(n1);
		      if (position != imag_values.end()) 
			value_imag = position->second;
		      else 
			value_imag = 0.0;
		
		      Complex value_complex = Complex(value_real, value_imag);
		      
		     
		      unsigned int dof_in_final_eigenvector = band_it->second * number_of_nodes + n1%number_of_nodes;

		      
		      result[pol] += value_complex * conj(eigen_vector_i[dof_in_initial_eigenvector]) * eigen_vector_f[dof_in_final_eigenvector];


		    }
		}


	      ierr = MatRestoreRow(P_real_p[pol]->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	    
	      ierr = MatRestoreRow(P_imag_p[pol]->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	    }
	}
    }
  
  return(result);
}





//----------------------------------------------------------------------------------------//


