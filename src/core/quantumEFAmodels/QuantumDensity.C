#include "QuantumDensity.h"


using namespace std;

//============================================//

QuantumDensity::QuantumDensity()
{

}

//============================================//

void QuantumDensity::set_options( QuantumDensity::options& options   )
{
  opt =  options;
}


//============================================//

QuantumDensity:: QuantumDensity( EnvelopFunctionApprox* model )
{
  quantum_model = model;
  real_space_density_size = quantum_model->get_number_of_active_cells();
}

//============================================//

QuantumDensity:: QuantumDensity( EnvelopFunctionApprox* model, QuantumDensity::options& options  )
{
  quantum_model = model;
  opt = options;
  real_space_density_size = quantum_model->get_number_of_active_cells();

}  



//============================================//
QuantumDensity:: ~QuantumDensity()
{
  delete eq;
  delete kmesh;
  
}


//============================================//
std::vector<double>& QuantumDensity::get_density(void)
{
  return(real_space_density);
}

//============================================//
void QuantumDensity::build_k_grid()
{


  //build mesh
  kmesh = new Mesh(k_dim);

  ElemType type;

  if (k_dim == 1)
    {
      type = EDGE3;
    }

  if (k_dim == 2)
    {
      type = QUAD8;
    }

  if (k_dim == 3)
    {
      type = HEX27;
    }


  MeshTools::Generation::build_cube (*kmesh, 
				     num_nodes[0], num_nodes[1], num_nodes[2], 
				     kmin[0], kmax[0], 
				     kmin[1], kmax[1], 
				     kmin[2], kmax[2],
				     type);
  cerr << kmin[0] <<"   "<<  kmax[0] <<"   "<<  kmin[1] <<"   "<<  kmax[1] << "  "<<  kmin[2] <<"   "<<
    kmax[2] << "  " <<  num_nodes[0] <<"   "<<  num_nodes[1] << "   "  << num_nodes[2] << "\n";


  
  rotate_mesh(kmesh, transform_matrix);

  kmesh->print_info();


  

}
//


//============================================//
void QuantumDensity::calculate_density()
{

 


  //----------------------------------------------------------------------//
  //build equation system object

 

 

  const DofMap& dof_map = system->get_dof_map();
    
  FEType fe_type = dof_map.variable_type(0);
    
  AutoPtr<FEBase> fe (FEBase::build(k_dim, fe_type));
    
  QGauss qrule (k_dim, FIRST);
    
  fe->attach_quadrature_rule (&qrule);

  
  const std::vector<Real>& JxW = fe->get_JxW();
  
  const std::vector<Point>& q_point = fe->get_xyz();
  
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  MeshBase::element_iterator       it     = kmesh->active_elements_begin();
  const MeshBase::element_iterator it_el  = kmesh->active_elements_end();


  real_space_density.resize(real_space_density_size,0.0);


  std::vector<unsigned int> dof_indices;


  for ( ; it != it_el ; ++it) 
    {

      const Elem* elem = *it;
      fe->reinit (elem);

      dof_map.dof_indices (elem, dof_indices, 0);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{//qp
	  unsigned int num_nodes =  phi.size();
	  for (unsigned int i = 0; i <  num_nodes; i++)
	    {
	      const Node* nd = elem->get_node(i);
	      map <const Node*, vector <double> > :: iterator node_it;

	      node_it = k_point_density.find(nd);
	      if (node_it != k_point_density.end())
		{
		  vector<double>&  dens_at_k_node = node_it->second;

		  for (unsigned int j = 0; j <  real_space_density_size; j++)
		    {
		      real_space_density[j] += phi[i][qp] * dens_at_k_node[j]*JxW[qp];
		    }
		    

		}
	      else
		{
		  cerr << "WARNING! QuantumDensity   node is missing\n";
		}
	      

	    }

	}
      
    }


  //--------------------------------------------------------------------------//

  // is in <cmath> : M_PI
  //double pi = 4.0 * atan(1.0);

  double factor = 1;

  for (short i = 0; i < k_dim; i++)  factor /= (2.0 * M_PI);

  for (unsigned int j = 0; j <  real_space_density_size; j++)  real_space_density[j] *= factor;
   


  //--------------------------------------------------------------------------//

}
//================================================================//
void  QuantumDensity::define_k_space(Tensor1 k_vector, unsigned int n)
{

  double norm_k = norm(k_vector);
  kmin[0] = -norm_k/2.0;  kmax[0] = norm_k/2.0; num_nodes[0] = n;

  
  //TODO hasto be changed!!
  Tensor1 basis1 = k_vector/norm_k;

  if (basis1(1) == 1)
    transform_matrix = Tensor2Sym(1);
  else
    {
      Tensor1 basis2;
     
      basis2(1) = 0;
      basis2(2) = -basis1(3);
      basis2(3) =  basis1(2);

      basis2 = basis2/norm(basis2);

      Tensor1 basis3 = vectorProduct(basis1, basis2);

      
      for (short i = 1; i < 4; i++)
	{
	  transform_matrix(i,1) = basis1(i);
	  transform_matrix(i,2) = basis2(i);
	  transform_matrix(i,3) = basis3(i);
	}
      
    }

 

  k_dim = 1;
}

//================================================================//
void QuantumDensity::define_k_space(Tensor1 k_vector1,unsigned int n,  Tensor1 k_vector2, unsigned int m)
{

  double norm_k1 = norm(k_vector1);
  kmin[0] = -norm_k1/2.0;  kmax[0] = norm_k1/2.0; num_nodes[0] = n;


  
  double norm_k2 = norm(k_vector2);
  kmin[1] = -norm_k2/2.0;  kmax[1] = norm_k2/2.0; num_nodes[1] = m;

  
  Tensor1 basis1 = k_vector1/norm_k1;
  Tensor1 basis2 = k_vector2/norm_k2;
  Tensor1 basis3 = vectorProduct(basis1, basis2);


  for (short i = 1; i < 4; i++)
    {
      transform_matrix(i,1) = basis1(i);
      transform_matrix(i,2) = basis2(i);
      transform_matrix(i,3) = basis3(i);
    }

  k_dim = 2;

}

//================================================================//

void QuantumDensity::define_k_space(Tensor1 k_vector1, unsigned int n, Tensor1 k_vector2, 
				    unsigned int m, Tensor1 k_vector3, unsigned int k)
{

  double norm_k1 = norm(k_vector1);
  kmin[0] = -norm_k1/2.0;  kmax[0] = norm_k1/2.0; num_nodes[0] = n;


  
  double norm_k2 = norm(k_vector2);
  kmin[1] = -norm_k2/2.0;  kmax[1] = norm_k2/2.0; num_nodes[1] = m;


  
  double norm_k3 = norm(k_vector3);
  kmin[2] = -norm_k3/2.0;  kmax[2] = norm_k3/2.0; num_nodes[2] = k;


  for (short i = 1; i < 4; i++)
    {
      transform_matrix(i,1) = k_vector1(i)/norm_k1; 
      transform_matrix(i,2) = k_vector2(i)/norm_k2;
      transform_matrix(i,3) = k_vector3(i)/norm_k3;
    }

  k_dim = 3;

}

//================================================================//
void QuantumDensity::calculate_at_each_k_point()
{

 

  MeshBase::node_iterator       it     = kmesh->active_nodes_begin();
  const MeshBase::node_iterator it_el  = kmesh->active_nodes_end();



  for ( ; it != it_el ; ++it) 
    {
      
      const Node*  nd  = *it;

      map<const Node*, vector<double> >::iterator it1;

      it1 =  k_point_density.find(nd);

      if (it1 ==  k_point_density.end() )

	{
	  double k_vector[3];


	  k_vector[0] = (*nd)(0);
	  k_vector[1] = (*nd)(1);
	  k_vector[2] = (*nd)(2);


	  
	  
	  std::cout << "Kvector (kx,ky,kz)   " << k_vector[0] << "   " <<  k_vector[1] << "    " <<  k_vector[2] << "\n";
	  
	  quantum_model->apply_k_vector(k_vector);
	  quantum_model->solve_eigen_value_problem(10);
	  

	  vector<double>  dens = quantum_model->calculate_convergent_density(opt.Temperature);

          unsigned int size_of_dens = dens.size();

	  for (unsigned int i2 = 0; i2 < size_of_dens  ; i2++)
	    dens[i2] *= opt.degeneracy;
     
     
	  k_point_density.insert( pair<const Node*, vector<double> > (nd, dens) );

	  double rho = (quantum_model->get_integrated_probability(opt.Temperature)) * opt.degeneracy;
	  
	  k_point_charge.insert(pair<const Node*, double > (nd, rho));
	  

	  

	}

    }


  
 
}


//==============================================================//

void QuantumDensity::rotate_mesh(Mesh* mesh, Tensor2Gen& RotMatrix)
{

  Tensor1 vec1;

  for (unsigned int n=0; n < mesh->n_nodes(); n++)
    {
       const Point p = mesh->node(n);
  
       vec1(1) = p(0);
       vec1(2) = p(1);
       vec1(3) = p(2);

       vec1 = RotMatrix * vec1;

       mesh->node(n) = Point( vec1(1), vec1(2), vec1(3) );

     }

}




//=============================================================//
void QuantumDensity::calculate_convergent_density()
{

  build_k_grid();

  eq = new EquationSystems(*kmesh);

  eq->add_system<LinearImplicitSystem> ("k-integration");
  
  system = &(eq->get_system<LinearImplicitSystem> ("k-integration"));
  
  system->add_variable("u", SECOND);

  
  //!system vector that contains charge density in k space from previous iteration
  NumericVector<Number>& old_density = system->add_vector("old density");


  eq->init();

  

  calculate_at_each_k_point();

 

  calculate_density();



  prepare_system_solution();


  


  if (opt.k_domain_refinement) 
    {
      //----------------------------------------
      //refinement block
      //---------------------------------------
     
     

      old_density = * (system->solution);
      

      MeshRefinement mesh_refinement(*kmesh);

      
      double norm_of_error = opt.relative_accuracy;

      for ( ; (norm_of_error >=  opt.relative_accuracy) ;  ) 
	{//for



	  if (opt.uniform_refinement)
	    mesh_refinement.uniformly_refine(1);
	  else
	    {
	      
	      ErrorVector error;
	      
	      
	      KellyErrorEstimator error_estimator;
	      
	      Tensor2Gen RotM_inv =  transform_matrix.transpose() ;
	      
	      rotate_mesh(kmesh,  RotM_inv );
	      
	      error_estimator.estimate_error (*system,error);
	      
	      rotate_mesh(kmesh, transform_matrix);
	      
	      
	      mesh_refinement.flag_elements_by_error_fraction (error,opt.refine_fraction,0.0, 10);
	      
	      
	      mesh_refinement.refine_and_coarsen_elements();
	      
	      kmesh->print_info();
	      
	      std::cout << "k-mesh after refinement  " << "\n";
	     
	      kmesh->print_info();

	      eq->reinit();

	      calculate_at_each_k_point();

	      calculate_density();

	      prepare_system_solution();


	      old_density.add(-1.0, *(system->solution));
	      
	      old_density.close();

	      double x1 = old_density.linfty_norm();

	      system->solution->close();

	      double x2 = ( system->solution->linfty_norm() );
	      
	      norm_of_error = x1/x2;


	      old_density = * (system->solution);

	      
	      
	    }

	}

    
      



    }//end of refinement block


}


//=================================================================//

void QuantumDensity::prepare_system_solution()
{

 

  unsigned int num_nodes = kmesh->n_nodes();

  system->solution->init(num_nodes);

 
  const  DofMap& dof_map = system->get_dof_map();
    
  FEType fe_type = dof_map.variable_type(0);


    
  AutoPtr<FEBase> fe (FEBase::build(k_dim, fe_type));


  QGauss qrule (k_dim, FIRST);
    
  fe->attach_quadrature_rule (&qrule);

  std::vector<unsigned int> dof_indices;


  MeshBase::element_iterator       it     = kmesh->active_elements_begin();
  const MeshBase::element_iterator it_el  = kmesh->active_elements_end();


  cerr << "----------------------------------------\n";




  for ( ; it != it_el ; ++it) 
    {

   
      const Elem* elem = *it;

      fe->reinit (elem);

     

      dof_map.dof_indices (elem, dof_indices);


    

      unsigned int n_nodes = elem->n_nodes();

    
      
      for (unsigned int n = 0; n < n_nodes; n++)
	{
	  const	  Node*  node =  elem->get_node(n);


	  cerr << (*node)(0) << "   " << (*node)(1) << "    " << (*node)(2) << "    " <<  k_point_charge[node] << "\n";

	 

	  system->solution->set( dof_indices[n] ,  k_point_charge[node]);
	  


	 

	}  
      

    }


}

//=================================================================//
const Mesh& QuantumDensity::get_k_mesh() const
{
  return(*kmesh);

}

//=================================================================// 
std::vector<double>   QuantumDensity::get_density_in_k_space(void)  const
{
  
  vector<double> result ;

  eq->build_solution_vector  	( result) ;     


 
  
  return(result);
}
