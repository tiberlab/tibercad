#include "KspaceIntegration.h"
using namespace std;

KspaceIntegration::KspaceIntegration()
{
  eq = NULL;

  system = NULL;

}

//-------------------------------------------------------//
KspaceIntegration::~KspaceIntegration()
{
  delete eq;
}
//-------------------------------------------------------//


void KspaceIntegration::calculate_density()
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
 
   


 
  

  std::vector<unsigned int> dof_indices;


  for ( ; it != it_el ; ++it) //loop over k space elements
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
	map <const Node*, map< const Elem*, double> > :: iterator node_it;

	node_it = k_point_density.find(nd);
	if (node_it != k_point_density.end())
	{
	  map<const Elem*, double>&  dens_at_k_node = node_it->second;
	  
	  map<const Elem*, double>::iterator   dens_at_k_node_it = dens_at_k_node.begin();
	  map<const Elem*, double>::iterator   dens_at_k_node_end = dens_at_k_node.end();

	

	  for ( ; dens_at_k_node_it != dens_at_k_node_end ; ++dens_at_k_node_it) //loop over real space elements
	  {
	    const Elem* el = dens_at_k_node_it->first;
	    
		     

	    double temp =  phi[i][qp] * (dens_at_k_node_it->second)*JxW[qp];
	   

	    real_space_density[el] +=  temp;
	    
	    
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
  {
    double factor = 1.0;

    for (short i = 0; i < k_dim; i++)  factor /= (2.0 * M_PI);
    
    map<const Elem*, double>::iterator   it = real_space_density.begin();
    map<const Elem*, double>::iterator   end_it = real_space_density.end();


    for (; it != end_it; ++it)  (it->second) *= factor;
  }


  //--------------------------------------------------------------------------//

}


//-------------------------------------------------------------------------------//
void KspaceIntegration::calculate_convergent_density()
{
  

  build_k_grid();

  eq = new EquationSystems(*kmesh);

  eq->add_system<LinearImplicitSystem> ("k-integration");
  
  system = &(eq->get_system<LinearImplicitSystem> ("k-integration"));
  
  system->add_variable("u", integration_order);

  
  //!system vector that contains charge density in k space from previous iteration
  NumericVector<Number>& old_density = system->add_vector("old density");

  
  eq->init();

  
  
 
  calculate_at_each_k_point();
 
 
  real_space_density.clear();



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
	      
	      
	mesh_refinement.flag_elements_by_elem_fraction (error,opt.refine_fraction,0.0, 10);

	      
	mesh_refinement.refine_and_coarsen_elements();
	     

	eq->reinit();

	calculate_at_each_k_point();
	
	real_space_density.clear();
	calculate_density();

	prepare_system_solution();


	old_density.add(-1.0, *(system->solution));
	      
	old_density.close();

	double x1 = old_density.linfty_norm();

	system->solution->close();

	double x2 = ( system->solution->linfty_norm() );
	      
	norm_of_error = x1/x2;


	old_density = *(system->solution);

	std::cout << "k space grid has " << kmesh->n_nodes() << " nodes " << flush;
	std::cout <<  "quantum density error " << norm_of_error << endl << flush;

	      
      }

    }

    
      



  }//end of refinement block


}

//----------------------------------------------------------------------------------------//
void KspaceIntegration::prepare_system_solution()
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


 


  for ( ; it != it_el ; ++it) 
  {

   
    const Elem* elem = *it;

    fe->reinit (elem);

     

    dof_map.dof_indices (elem, dof_indices);
      

    

    unsigned int n_nodes = elem->n_nodes();

    
      
    for (unsigned int n = 0; n < n_nodes; n++)
    {
      const	  Node*  node =  elem->get_node(n);
      
#ifdef DEBUG
      cerr << (*node)(0) << "   " << (*node)(1) << "    " << (*node)(2) << "    " <<  k_point_charge[node] << "\n";
#endif
      

      system->solution->set( dof_indices[n] ,  k_point_charge[node]);
 	  


	 

    }  
      

  }


}
//------------------------------------------------------------------------------------//

void KspaceIntegration::parse_options( )
{

  const ModelOptions& mod_opt = get_options();
 


  opt.uniform_refinement      = mod_opt.get_option("uniform_refinement",false);

  opt.refine_fraction         = mod_opt.get_option("refine_fraction", 0.3);
  opt.maximum_ref_level       = mod_opt.get_option("maximum_ref_level", 8);
  opt.relative_accuracy       = mod_opt.get_option("relative_accuracy", 1e-2);

  opt.degeneracy                = mod_opt.get_option("degeneracy",1);
  opt.k_domain_refinement       = mod_opt.get_option("refine_k_space", false);



}
//------------------------------------------------------------------------------------//

void KspaceIntegration::do_solve( )
{


  parse_options();



  calculate_convergent_density();



}

//--------------------------------------------------------------------------------------//
void KspaceIntegration::do_init(void)
{
  Kspace::do_init();

}
