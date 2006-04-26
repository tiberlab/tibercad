
#include "EnvelopFunctionApprox.h"
using namespace std;
EnvelopFunctionApprox:: EnvelopFunctionApprox(options& opt1, Mesh& mesh, MeshData& mesh_data_in)
{
  //Initoalization
  opt = opt1;

  es = new EquationSystems(mesh);

  es->add_system<LinearImplicitSystem> ("Schroedinger");

  dim = mesh.mesh_dimension();
  
  //---------------------------------------------------------------------------------------
  //add variables
  psi_name.clear();
  for (short i = 0; i < opt.number_of_bands; i++)
    {
      std::ostringstream var_str;
      var_str << "psi" << i ;
      string name = var_str.str();
      psi_name.push_back(name);

      es->get_system("Schroedinger").add_variable(name,FIRST);
    } 

  //add matrixes

  //-----------------------------------------------------------------------------------------------------//
  LinearImplicitSystem& ls = static_cast<LinearImplicitSystem&>(es->get_system("Schroedinger"));

  ls.add_matrix("Ham_real"); //add matrix for a real part of the Hamiltonian

  Ham_real = & (ls.get_matrix("Ham_real"));

  ls.add_matrix("Ham_imag");//add matrix for an imaginary part of the Hamiltonian

  Ham_imag = &(  ls .get_matrix("Ham_imag") );

  ls.add_matrix("S_real"); //add matrix for S matrix

  S_real = &( ls.add_matrix("S_real") );

  //---------------------------------------------------------------------------------------------------------//
  //My Jacobian 

   my_Jacobian = 1.0;
   for (short i = 1; i <= dim; i++)
     my_Jacobian *= opt.length_scale;

   //--------------------------------------------------------------------------------------------------------//
  

   es->init();

   //-------------------------------------------------------------------------------------------------------//

   //material list

   meshdata = &mesh_data_in;
   assemble_material_list();

   //------------------------------------------------------------------------------------------------------//
   
}

//============================================================//
void EnvelopFunctionApprox::set_material_parameters(std::vector<EFAbulkHamiltonian*>&  bulkHamiltonian_in)
{
  bulkHamiltonian = bulkHamiltonian_in;
}



//===========================================================//
void EnvelopFunctionApprox::calculate_Hamiltonian_and_S(void)
{



 LinearImplicitSystem& system = es->get_system<LinearImplicitSystem>("Schroedinger");

 

 Mesh& mesh = system.get_mesh(); 

 vector<unsigned int> psivar(opt.number_of_bands);
 //get numbers of variables
 for (unsigned int i = 0; i < opt.number_of_bands; i++)
   {
     psivar[i] = system.variable_number(psi_name[i]);
   }

 DofMap& dof_map = system.get_dof_map();

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


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  unsigned int el_number = 0;
  for ( ; el != end_el ; ++el) 
    {//el
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;
      const unsigned int mat = material_of_elem[el_number];

      std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&  model_Ham = ( (bulkHamiltonian[mat])->get_Hamiltonian() );

      dof_map.dof_indices (elem, dof_indices); 
      const unsigned int n_dofs   = dof_indices.size();
      fe->reinit (elem);

      ham_real.resize(n_dofs, n_dofs);
      ham_imag.resize(n_dofs, n_dofs);
      s_real.resize(n_dofs, n_dofs);

      complex<double> operator_sign;

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{//qp
	  for (unsigned int band1 = 0; band1 < opt.number_of_bands; band1++)
	    {//band1
	      dof_map.dof_indices (elem, dof_indices_component, psivar[band1]);
	      const unsigned int n_psi_dofs = dof_indices_component.size();
 
	       for (unsigned int band2 = 0; band2 < opt.number_of_bands; band2++)
		 {//band2
		   
		   if (band1 < band2)
		     operator_sign = Complex(0.0, -1.0);
		   else
		     operator_sign = Complex(0.0, 1.0);

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
			       * operator_sign /opt.length_scale;

                           //linear right

			   for (short i = 0; i < dim; i++)
			     value += JxW[qp]* dphi[p2][qp](i) * phi[p1][qp] * model_Ham[band1][band2].linear_right[i] 
			       * operator_sign /opt.length_scale;
			
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



  Ham_real->print_matlab("ham_r_matlab.m");
  Ham_imag->print_matlab("ham_i_matlab.m");
  S_real->print_matlab("s.m");


  dof_map.print_dof_constraints();
     
}
//============================================================//
void EnvelopFunctionApprox::assemble_material_list(void)
{
Mesh& mesh = es->get_mesh();
 

  const unsigned int N_elem = mesh.n_active_elem();

 

  material_of_elem.resize( N_elem );

 
  
  MeshBase::const_element_iterator el  = mesh.active_elements_begin();
  MeshBase::const_element_iterator end_el = mesh.active_elements_end();

    

  unsigned int el_number = 0;
  
  for ( ; el != end_el ; ++el) 
    {
      const Elem* elem = *el;
      unsigned int mat;
      //-------------------------------------
      
      mat  = ( (unsigned int) (*meshdata)(elem->top_parent(),0) )  - 1;
      material_of_elem[el_number] = mat;
      Point p = elem->centroid();
   
      el_number++;
      
    }

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


 

  int size_matrix = S_real->n();
  out_int = *(reinterpret_cast<unsigned int*> (& size_matrix) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);

  out_int = *(reinterpret_cast<unsigned int*> (& size_matrix) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
 
  PetscMatrix<Number>* p_matrix = static_cast<PetscMatrix<Number>* >(S_real);

  p_matrix->close();
  //-------------------------------------------------------------------------------------
  //how many non-zero elements do we have?

  int Number_of_elements = 0;

  

  for (int row = 0 ; row < size_matrix; row++)
    {
      int ierr = 0;
      const  PetscScalar *petsc_row_vals;
      const  PetscInt *petsc_cols;
      int n_cols = 0;

      ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      Number_of_elements += n_cols;
      

      ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

    }

  std::cout << "We have got " << Number_of_elements << " non-zero elements in S matrix\n"; 

  out_int = *(reinterpret_cast<unsigned int*> (& Number_of_elements) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
  //-----------------------------------------------------------------------------------
  //write number of columns in each row
  for (int row = 0 ; row < size_matrix; row++)
    {
      int ierr = 0;
      const  PetscScalar *petsc_row_vals;
      const  PetscInt *petsc_cols;
      int n_cols = 0;

      ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      out_int = *(reinterpret_cast<unsigned int*> (& n_cols) );  endian_swap(out_int);
      out.write(  reinterpret_cast<char *>( & out_int ), int_size);

      ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

    }
  //===============================================================//
  //write number of columns in each row
  for (int row = 0 ; row < size_matrix; row++)
    {
      int ierr = 0;
      const  PetscScalar *petsc_row_vals;
      const PetscInt *petsc_cols;
      int n_cols = 0;

      ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      for (short col = 0; col < n_cols; col++)
	{
	  int col_number = petsc_cols[col];
	  out_int = *(reinterpret_cast<unsigned int*> (& col_number) );  endian_swap(out_int);
	  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
	}

      ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

    }

  //write data
  for (int row = 0 ; row < size_matrix; row++)
    {
      int ierr = 0;
      const  PetscScalar *petsc_row_vals;
      const PetscInt *petsc_cols;
      int n_cols = 0;

      ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      for (short col = 0; col < n_cols; col++)
	{
	  double value = petsc_row_vals[col];
	  double zero = 0.0;
	  
	  out_long_long = *(reinterpret_cast<unsigned long long*> (& value) );  endian_swap(out_long_long);
	  out.write(  reinterpret_cast<char *>( & out_long_long ), double_size);

	  out_long_long = *(reinterpret_cast<unsigned long long*> (& zero) );  endian_swap(out_long_long);
	  out.write(  reinterpret_cast<char *>( & out_long_long ), double_size);
	  
	}

      ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

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

  out_int = *(reinterpret_cast<unsigned int*> (& size_matrix) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);

  out_int = *(reinterpret_cast<unsigned int*> (& size_matrix) );  endian_swap(out_int);
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

      for (int i = 0; i < n_cols_real; i++) real_column.insert(petsc_cols_real[i]);

      ierr = MatGetRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      imag_column.clear();
      for (int i = 0; i < n_cols_real; i++) imag_column.insert(petsc_cols_imag[i]);

      set_union(real_column.begin(), real_column.end(), imag_column.begin(), imag_column.end(), com_ins);

      Number_of_elements += complex_column.size();

      ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      

    }
  std::cout << "We have got " << Number_of_elements << " non-zero elements in the Hamiltonian matrix\n"; 

  out_int = *(reinterpret_cast<unsigned int*> (& Number_of_elements) );  endian_swap(out_int);
  out.write(  reinterpret_cast<char *>( & out_int ), int_size);
  //----------------------------------------------------------------------------------------------------//
  //write number of columns in each row
  for (int row = 0 ; row < size_matrix; row++)
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

      for (int i = 0; i < n_cols_real; i++) real_column.insert(petsc_cols_real[i]);

      ierr = MatGetRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      imag_column.clear();
      for (int i = 0; i < n_cols_real; i++) imag_column.insert(petsc_cols_imag[i]);

      set_union(real_column.begin(), real_column.end(), imag_column.begin(), imag_column.end(), com_ins);

      int n1 = complex_column.size();

      out_int = *(reinterpret_cast<unsigned int*> (& n1) );  endian_swap(out_int);
      out.write(  reinterpret_cast<char *>( & out_int ), int_size);

    
      ierr = MatRestoreRow(H_real_matrix->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = MatRestoreRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      

    }    
  //-----------------------------------------------------------------------------------------------------//
  //write number of columns in each row
  for (int row = 0 ; row < size_matrix; row++)
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

      for (int i = 0; i < n_cols_real; i++) real_column.insert(petsc_cols_real[i]);

      ierr = MatGetRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      imag_column.clear();
      for (int i = 0; i < n_cols_real; i++) imag_column.insert(petsc_cols_imag[i]);

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
  //-------------------------------------------------------------------------------
  //write data of columns in each row
  
  for (int row = 0 ; row < size_matrix; row++)
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
	  real_column.insert(petsc_cols_real[i]);
	  real_values.insert(make_pair(petsc_cols_real[i],petsc_row_vals_real[i] ));
	}

      ierr = MatGetRow(H_imag_matrix->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      imag_column.clear();
      imag_values.clear();
      for (int i = 0; i < n_cols_real; i++)
	{ 
	  imag_column.insert(petsc_cols_imag[i]);
	  imag_values.insert(make_pair(petsc_cols_imag[i],petsc_row_vals_imag[i] ));
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
  
  //------------------------------------------------------------------------------
}
//===============================================================================//

void EnvelopFunctionApprox::solve_eigen_value_problem(unsigned int ev_number)
{

  create_dirichlet_dofs();

  make_constraints();


  


  calculate_Hamiltonian_and_S(); //calculate Hamiltonian and S matrix

  //----------------------------------------------------------------
  //write matrixes
  save_S_matrix("S.out");
  save_H_matrix("H.out");
  //----------------------------------------------------------------

  std::ostringstream  command_line;

  command_line <<  "eigen_solver  -f1 H.out   -f2 S.out  -eps_gen_hermitian -eps_smallest_magnitude ";
  command_line <<  "   -eps_nev   " << ev_number;
  command_line <<  "   -eps_type  " << opt.solver;
  command_line <<  "   -eps_tol   " << opt.eigen_solver_tolerance;
  command_line <<  "   -eps_monitor \n";

  system( (command_line.str()).c_str());

  read_SLEPC_solution();

}

//=============================================================//
void EnvelopFunctionApprox::read_SLEPC_solution( )
{
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
  cout << " Number of converged solutions  " << number_of_converged_solutions << "\n";

  solution.resize(number_of_converged_solutions);
  //--------------------------------------------------------------------
  //read eigenvalues
  for (unsigned i = 0; i < number_of_converged_solutions; i++)
    {
      
      
      file_eigvals.read(buffer_double, double_size);
      fict = *( reinterpret_cast<unsigned long long*>( buffer_double) ); endian_swap(fict);
      solution[i].eigen_energy = *(  reinterpret_cast<double*>( &fict ) );
      cerr << solution[i].eigen_energy << "\n";

      file_eigvals.read(buffer_double, double_size);
      
      
    }
  
  //--------------------------------------------------------------------
  //read eigenvectors
  string fname_eigvects = "eigvects_SLEPC.out";


  std::ifstream file_eigvects ( fname_eigvects.c_str() );

  assert (file_eigvects.good());

  for (unsigned i = 0; i < number_of_converged_solutions; i++)
    {
      file_eigvects.read(buffer, int_size);
      file_eigvects.read(buffer, int_size);
      unsigned int vector_size =  *(reinterpret_cast<unsigned int*> ( buffer));  endian_swap(vector_size);
      (solution[i].eigen_vector).resize(vector_size);

      for (unsigned j = 0; j < vector_size; j++)
	{
          double re, im;
	  file_eigvects.read(buffer_double, double_size);
	  fict = *( reinterpret_cast<unsigned long long*>( buffer_double) ); endian_swap(fict);
	  re   = *(  reinterpret_cast<double*>( &fict ) );

	  file_eigvects.read(buffer_double, double_size);
	  fict = *( reinterpret_cast<unsigned long long*>( buffer_double) ); endian_swap(fict);
	  im   = *(  reinterpret_cast<double*>( &fict ) );

	  solution[i].eigen_vector[j] = Complex(re,im);
	}
    }
  
}
//=============================================================//


//=======================================================================//

void EnvelopFunctionApprox::output_eigen_functions(unsigned int state_number,  std::string& filename)
{
  //===========================================
  const Mesh& mesh = es->get_mesh();
  LinearImplicitSystem& system = es->get_system<LinearImplicitSystem>("Schroedinger");
  DofMap& dof_map = system.get_dof_map();

  MeshBase::const_node_iterator       nd     = mesh.active_nodes_begin();
  const MeshBase::const_node_iterator nd_el  = mesh.active_nodes_end();

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
  
  unsigned int output_size = ( solution[state_number].eigen_vector.size()) * 2;

  vector<double>  psi_data(output_size);
 
 
  MeshBase::const_element_iterator it = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh.active_local_elements_end();

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
	      Complex value =  (solution[state_number].eigen_vector[psi_index * number_of_points + dof_indices[n] ]);
	      
	      psi_data[psi_index*2 + node_id*number_output_data] = value.real();
	      
	      psi_data[psi_index*2 + 1  + node_id*number_output_data] = value.imag();
	      
	      
	      
	    }
	}
    }


    


  //std :: cout << filename << "\n";


  if (opt.output_type == "GMV")     GMVIO(mesh).write_nodal_data(filename, psi_data, output_names);

  if (opt.output_type == "tecplot") TecplotIO(mesh,false).write_nodal_data(filename,psi_data,output_names);

  

}

//=======================================================================//
void EnvelopFunctionApprox::assign_mesh_data(MeshData& mesh_data_in)
{
   meshdata = &mesh_data_in;
   assemble_material_list();
}

//=======================================================================//

void EnvelopFunctionApprox::define_diriclet_nodes(std::vector<unsigned int>&  dirichlet_nodes_input)
{
  dirichlet_nodes = dirichlet_nodes_input;


}

//======================================================================//
void  EnvelopFunctionApprox::create_dirichlet_dofs( )
{
  const Mesh& mesh = es->get_mesh();
  LinearImplicitSystem& system = es->get_system<LinearImplicitSystem>("Schroedinger");
  DofMap& dof_map = system.get_dof_map();

  MeshBase::const_element_iterator it = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh.active_local_elements_end();

  dirichlet_dofs.clear();
  // insert_iterator<set<int> >  dir_ins(dirichlet_dofs,dirichlet_dof.begin() ); 


  const std::vector<unsigned int> :: const_iterator  n_begin = dirichlet_nodes.begin();
  const std::vector<unsigned int> :: const_iterator  n_end   = dirichlet_nodes.end();
  std::vector<unsigned int> :: const_iterator n_it;


  std::vector<unsigned int> dof_indices;

  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
	{ 
	  unsigned int  node_id =  elem->node(n);
	  //does a node belong to a a dirichlet nodes set?
	  if (find(n_begin, n_end, node_id) != n_end)
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
void EnvelopFunctionApprox::make_constraints(void)
{
  const Mesh& mesh = es->get_mesh();
  LinearImplicitSystem& system = es->get_system<LinearImplicitSystem>("Schroedinger");
  DofMap& dof_map = system.get_dof_map();
  
  DofConstraintRow constraint;
  constraint.clear();
  
  std::set<unsigned int> :: const_iterator dof_begin = dirichlet_dofs.begin();
  std::set<unsigned int> :: const_iterator dof_end   = dirichlet_dofs.end();
  std::set<unsigned int> :: iterator dof_it;

  for (dof_it = dof_begin; dof_it !=  dof_end; dof_it++)
    {
      unsigned int n_dof = *dof_it;      
      dof_map.add_constraint_row (n_dof, constraint) ; 
    }

}
//=======================================================================//

//=======================================================================//
EnvelopFunctionApprox:: ~EnvelopFunctionApprox(void)
{
  delete(es);
}



