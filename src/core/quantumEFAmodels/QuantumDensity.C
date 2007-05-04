#include "QuantumDensity.h"
#include "SimulationEnvironment.h"
#include "Control.h"
#include "gnuplot_io.h"
using namespace std;

void QuantumDensity::do_plot (void)
{
  //---------------------------------------------------------------------------
  //standard output
  SimulationInterface::do_plot();
  //---------------------------------------------------------------------------
  //k-space output
  const Device& dev = get_environment().get_device();

  string suffix = get_control().get_filename_suffix();
  string outdir = get_control().get_output_dir();
  string format = get_control().get_output_format();

  string suff;
  if (format == "gmv")
    suff = ".gmv";
  else if (format == "ise")
    suff = ".plt";

  const std::set< std::string >& plotvariables = get_control().get_plotvariables();

  if (plotvariables.find("k-space_density") != plotvariables.end())
  {
    string filename(outdir + "/" + get_name() +
        "_k_space" + suffix + suff);

   

    vector<string> names(1,"density[atomic_units]");

    const vector<double> results = get_density_in_k_space();
  


    if (format == "gmv")
      GMVIO(get_k_mesh()).write_nodal_data(filename, results, names);
    else if (format == "gnuplot")
      GnuPlotIO(get_k_mesh()).write_nodal_data(filename, results, names);
    else if (format == "ise")
      TecplotIO(get_k_mesh()).write_nodal_data(filename, results, names);
    else
    {
      cout << "Output format not supported. Falling back to GMV." << endl;
      GMVIO(get_k_mesh()).write_nodal_data(filename, results, names);
    }

  }

  //----------------------------------------------------------------------------
}



void QuantumDensity::get_particle_density(const Elem* element, const std::vector<double>& quad_points, std::vector<double> density)
{

  unsigned int n = quad_points.size();
  density.resize(n,0.0);

  map< const Elem*, double >::iterator it;
  it = real_space_density.find(element);
  if (it != real_space_density.end())
  {//if the element elem is active we put the same density for all the point
    unsigned int n = quad_points.size();
   
    for (unsigned int i = 0; i < n; i++)
      density[i] = it->second;
  }
  else
  {//if the element is not included, we have to check his children or parents
    //-------------------------------------------------------------------------
    //1) may be it has a parent the belongs to the  real_space_density map
    const Elem* el1 = element->parent();
    bool out = false;
    bool found = false;
    
    while ( !out )
    {
      if ( el1 != NULL )
      {
	it = real_space_density.find(el1);
	if ( it != real_space_density.end() )
	{
	  unsigned int n = quad_points.size();


	  for (unsigned int i = 0; i < n; i++)
	    density[i] = it->second;
	  
	  out = true;
	  found = true;
	}
	else
	{
	  el1 = el1->parent();
	}
      }
      else
      {
	out = true;
      }
    }
   
    if (!found) //2) may be it has a child that belongs to the real_space_density map
    {
      std::vector< const Elem * > active_children;
      element -> active_family_tree ( active_children, true);
      unsigned int n = active_children.size();

      for (unsigned int i1 = 0; ( i1 < n  ); i1++)
      {
	it  =  real_space_density.find(active_children[i1]);
	if (it !=  real_space_density.end() )
	{
	  //we have to check if this child contains a quadrature points
	  const unsigned int n = quad_points.size();
	  
	  for (unsigned int i = 0; i < n; i++)
	    if (active_children[i1]->contains_point(quad_points[i]) )
	    {
	      density[i] = it->second;
     
	    }
	  
	}
	
      

      }
    }
  }
}

//==================================================================================//
void QuantumDensity::build_elemental_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend)
{
  if (variables.find("quantum_density") != variables.end() )
  {
    legend.resize(1, "particle_density[cm^-3]");

    //------------------------------------------------------------------------------------//
    //!from map to vector
    Mesh & mesh = get_environment().get_mesh();
    
    map<const Elem*, double>::iterator it = real_space_density.begin();

    MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
    
    results.resize(real_space_density_size);

    unsigned int el_number = 0;

    for (; el !=end_el ; ++el)
    {
      const Elem* elem = *el; 

      it = real_space_density.find(elem);

      if (it != real_space_density.end())
	results[el_number] = it->second;
      else
	cerr << "WARNING: not all elements have a calculated density\n";

      el_number++;
    }
    //-------------------------------------------------------------------------------------//

    //now I have to transform from atomic units to [cm^-3]    

    const double coeff =  1.0/ ( (Constants::bohr_radius) * (Constants::bohr_radius) * (Constants::bohr_radius) * 1.0e6 );

    for (unsigned int i = 0; i < real_space_density_size; i++)
    {
      results[i] *= coeff;
    }
  


    //-------------------------------------------------------------------------------------//

  }

}

//===================================================================================//

QuantumDensity::QuantumDensity()
{
  quantum_model = NULL;

  system = NULL;

  eq = NULL;


}

//============================================//
void QuantumDensity::do_init( )
{
  

  const ModelOptions& mod_opt = get_options();

  Kspace::do_init();

  //--kspace domain---------------------------//

  //---------quantum model---------------------------------------------------------------------//
  
  std::string quantum_simul_name;
  if (mod_opt.find_option("quantum_simulation"))
  {
    quantum_simul_name = mod_opt.get_option("quantum_simulation","");
    quantum_model = dynamic_cast<  EnvelopFunctionApprox* > ( find_simulation( quantum_simul_name  )   );
    if (quantum_model == NULL)
      throw  InitFailedException("QuantumDensity: quantum_simulation " + quantum_simul_name + " does not exist");
  }
  else
  {
    throw  InitFailedException("QuantumDensity: quantum_simulation  has to be specified");
  }

 
  //--------------------------------------------------------------------------------------------//


}


//============================================//
void QuantumDensity::do_solve ( )
{


 parse_options();



 calculate_convergent_density();



}
//============================================//
void QuantumDensity::parse_options( )
{
 const ModelOptions& mod_opt = get_options();
 
 opt.Temperature             = mod_opt.get_option("Temperature", opt.Temperature);
 opt.log_output              = mod_opt.get_option("log_output", false);
 opt.uniform_refinement      = mod_opt.get_option("uniform_refinement",false);

 opt.refine_fraction         = mod_opt.get_option("refine_fraction", 0.3);
 opt.maximum_ref_level       = mod_opt.get_option("maximum_ref_level", 8);
 opt.relative_accuracy       = mod_opt.get_option("relative_accuracy", 1e-2);

 opt.degeneracy                = mod_opt.get_option("degeneracy",1);
 opt.k_domain_refinement       = mod_opt.get_option("refine_k_space", false);
 opt.intial_eigenstates_number = mod_opt.get_option("intial_eigenstates_number", 6);

 {//-------------------------------------------------------------------------------------//
   std::string  job_name = mod_opt.get_option("job","density");
   if (job_name == "density")
     opt.job = DENSITY;
   else if (job_name == "dispersion")
     opt.job = DISPERSION;
   else
     throw InitFailedException( "QuantumDensity: Incorrect job " + job_name ); 
   //-------------------------------------------------------------------------------------//
 }
}



//============================================//
QuantumDensity:: ~QuantumDensity()
{
  delete eq;


  
}








//========================================================================//
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
 
   
  real_space_density.clear();

 
  

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

//================================================================//
void QuantumDensity::calculate_at_each_k_point()
{

 

  MeshBase::node_iterator       it     = kmesh->active_nodes_begin();
  const MeshBase::node_iterator it_el  = kmesh->active_nodes_end();



  for ( ; it != it_el ; ++it) 
    {
      
      const Node*  nd  = *it;

      map<const Node*, map< const Elem*, double> >::iterator it1;

      it1 =  k_point_density.find(nd);

      if (it1 ==  k_point_density.end() )

      {
	vector<double> k_vector(3, 0.0);


	k_vector[0] = (*nd)(0);
	k_vector[1] = (*nd)(1);
	k_vector[2] = (*nd)(2);
	  

	ModelOptions quantum_model_opts;
	  

	quantum_model_opts.set_option("k_vector",  k_vector);
	quantum_model_opts.set_option("initial_eigenstates_number",opt.intial_eigenstates_number ); 
	quantum_model_opts["job"] = "density";

 
	  
	quantum_model->set_options(quantum_model_opts);


	quantum_model->solve();

	map<const Elem*, double> dens = quantum_model->get_density();

	if (it == kmesh->active_nodes_begin()) real_space_density_size = dens.size();
	  
	map<const Elem*, double>::iterator density_it = dens.begin();
	map<const Elem*, double>::iterator density_end   = dens.end();
	


	for ( ; density_it != density_end  ; ++density_it )
	{
	  density_it->second *= opt.degeneracy;	 
	}
     
	k_point_density.insert( pair< const Node*, map<const Elem*, double> > (nd, dens) );

	double rho = (quantum_model->get_integrated_probability()) * opt.degeneracy;

	 
	  
	k_point_charge.insert(pair<const Node*, double > (nd, rho));
	  

	  

      }

    }


#ifdef DEBUG
  cerr << "Schroedinger equation at each point is solved\n";
#endif 
 
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


	      old_density = *(system->solution);

	      std::cout << "k space grid has " << kmesh->n_nodes() << " nodes " << flush;
	      std::cout <<  "quantum density error " << norm_of_error << endl << flush;

	      
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


//=================================================================// 
std::vector<double>   QuantumDensity::get_density_in_k_space(void)  const
{
  
  vector<double> result ;

  eq->build_solution_vector  	( result) ;     
  
  return(result);
}
