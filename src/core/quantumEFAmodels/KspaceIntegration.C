#include "KspaceIntegration.h"
#include "Control.h"
#include "gnuplot_io.h"
#include "SimulationEnvironment.h"
#include "GMVIO_cell.h"
#include "tecplot_IO_cell.h"
#include "VTKIO.h"

 
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
void KspaceIntegration::do_plot()
{
 //---------------------------------------------------------------------------
  //standard output
  SimulationInterface::do_plot();
  //---------------------------------------------------------------------------
  //k-space output
  if (kmesh != NULL)
  {

    const Device& dev = get_environment().get_device();

    string suffix = get_control().get_filename_suffix();
    string outdir = get_control().get_output_dir();
    string format = get_control().get_output_format();

    string suff;

    if (format == "gmv")
      suff = ".gmv";
    else if (format == "ise")
      suff = ".plt";
    else if (format == "vtk")
      suff = ".vtk";
    else 
      suff = ".gmv"; 

    const std::set< std::string >& plotvariables = get_control().get_plotvariables();

    if (plotvariables.find("k-space") != plotvariables.end())
    {
      string filename(outdir + "/" + get_name() +
		      "_k_space" + suffix + additional_name_suffix + suff);

   

      vector<string> names(1,"value[atomicUnits]");

   

      const vector<double> results = get_density_in_k_space();
  

      if (format == "gmv")
	GMVIO_cell(get_k_mesh()).write_ascii_cell_data(filename, results, names);
      else if (format == "ise")
	TecplotIO_cell(get_k_mesh()).write_cell_data(filename, results, names);
      else if (format == "vtk")
	VTKIO(get_k_mesh()).write_elemental_data(filename, results, names); 
      else
      {
	cout << "Output format not supported. Falling back to GMV." << endl;
	GMVIO_cell(get_k_mesh()).write_ascii_cell_data(filename, results, names);
      }
      
    }
  }
  //----------------------------------------------------------------------------

}
//-------------------------------------------------------//
void KspaceIntegration::calculate_density()
{

 

  const DofMap& dof_map = system->get_dof_map();
    
  FEType fe_type = dof_map.variable_type(0);

 
    
  AutoPtr<FEBase> fe (FEBase::build(k_dim, fe_type));

  QGauss qrule (k_dim, THIRD);
    
  fe->attach_quadrature_rule (&qrule);

 
  
  const std::vector<Real>& JxW = fe->get_JxW();
  
  const std::vector<Point>& q_point = fe->get_xyz();
  
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

 
 
  MeshBase::element_iterator       it     = kmesh->active_elements_begin();
  const MeshBase::element_iterator it_el  = kmesh->active_elements_end();
 
   

 
  double factor = 1.0;
  {
    for (short i = 0; i < k_dim; i++)  factor /= (2.0 * M_PI);
    
    factor *= get_degeneracy_factor() * opt.degeneracy;
  }

 

  std::vector<unsigned int> dof_indices;


 


  for ( ; it != it_el ; ++it) //loop over k space elements
  {
    const Elem* elem = *it;
    

    std::map<const Elem*, std::map <const Elem*, double> >::iterator it_k_elem;
    
    it_k_elem = kspace_local_density.find(elem);

    if (it_k_elem == kspace_local_density.end())
    {  

      map<const Elem*, double>  dens_at_k_elem;

      fe->reinit (elem);

      dof_map.dof_indices (elem, dof_indices, 0);
    
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {//qp
	
      

	map<const Elem*, double>  dens_at_k_point ;

	double integrated_quantity;

	calculate_for_k_point( q_point[qp], dens_at_k_point, integrated_quantity);
	

	kspace_integral[elem] += integrated_quantity * JxW[qp] * factor;


	  

	
	map<const Elem*, double>::iterator   dens_at_k_node_it = dens_at_k_point.begin();
	map<const Elem*, double>::iterator   dens_at_k_node_end = dens_at_k_point.end();
	
	

	for ( ; dens_at_k_node_it != dens_at_k_node_end ; ++dens_at_k_node_it) //loop over real space elements
	{
	  const Elem* el = dens_at_k_node_it->first;
	  
	  dens_at_k_elem[el] += (dens_at_k_node_it->second)*JxW[qp] * factor;
	}
	  
 
	     
	
	
      }

      kspace_local_density.insert(pair<const Elem*,map<const Elem*, double> >(elem,dens_at_k_elem));
    }
  }

  //--------------------------------------------------------------------------//
  {//real space density calculation
  
  
   

    MeshBase::element_iterator       it_k_space     = kmesh->active_elements_begin();
    const MeshBase::element_iterator it2  = kmesh->active_elements_end();
    

    for ( ; it_k_space != it2 ; ++it_k_space)
    {
      const Elem* el = *it_k_space;

      map<const Elem*, double>& k_cell_map = kspace_local_density[el];

      map<const Elem*, double>::iterator it_real_space =   k_cell_map.begin();
      map<const Elem*, double>::iterator it_3 =   k_cell_map.end();


      for ( ; it_real_space != it_3 ; ++it_real_space )
      {
	real_space_density[it_real_space->first] += it_real_space->second;
      }

    }



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

  
 

  
  eq->init();

  
  kspace_local_density.clear();
 
  kspace_integral.clear();

  volume.clear();
 
 
  real_space_density.clear();

  calculate_volumes();

  calculate_density();

 

  double x1;
  double x2;
 
  

  if (opt.k_domain_refinement) 
  {
    //----------------------------------------
    //refinement block
    //---------------------------------------
     
     

      

    MeshRefinement mesh_refinement(*kmesh);

      
    double norm_of_error = opt.relative_accuracy;

    
 
    for ( ; (norm_of_error >=  opt.relative_accuracy) ;  ) 
    {//for

      if (opt.uniform_refinement)
	mesh_refinement.uniformly_refine(1);
      else
      {
	      

	      
	ErrorVector error = ErrorVector(kmesh->n_elem(), kmesh);  
	
	      
	Tensor2Gen RotM_inv =  transform_matrix.transpose() ;
	      
	rotate_mesh(kmesh,  RotM_inv );
	
	estimate_error_for_refinement(error);
	
	rotate_mesh(kmesh, transform_matrix);
	      

        mesh_refinement.refine_fraction() = opt.refine_fraction;

	mesh_refinement.max_h_level() = 10;

	mesh_refinement.coarsen_fraction() = 0.0;
	
	      
	mesh_refinement.flag_elements_by_error_fraction (error);

	      
	mesh_refinement.refine_and_coarsen_elements();
	     

	if (opt.log_output)
	{
	  kmesh->print_info();
	  cout.flush();
	}

	eq->reinit();


	old_real_space_density = real_space_density;

	real_space_density.clear();

	calculate_volumes();

	calculate_density();


	norm_of_error = estimate_error();


	if (opt.log_output)
	{
	  std::cerr << "k space grid has " << kmesh->n_nodes() << " nodes " << flush;
	  std::cerr <<  "quantum density error " << norm_of_error << endl << flush;
	}
	      
      }

    }

    
      



  }//end of refinement block


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
  opt.log_output                = mod_opt.get_option("log_output",  false);


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
//--------------------------------------------------------------------------------------//
void KspaceIntegration::estimate_error_for_refinement(ErrorVector& error)
{
  


  std::fill (error.begin(), error.end(), 0.0);

 
  MeshBase::const_element_iterator       elem_it1  = kmesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end1 = kmesh->active_elements_end(); 
  
  for (; elem_it1 != elem_end1; ++elem_it1)
  {

   
    const Elem* el = *elem_it1;
    const unsigned int el_id = el->id();


    error[el_id] = abs(kspace_integral[el]); //test
    
  

  }


 

}

//--------------------------------------------------------------------------------------//

void KspaceIntegration::calculate_volumes(void)
{
 
  MeshBase::const_element_iterator       elem_it  = kmesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end = kmesh->active_elements_end(); 
  

 

  for (; elem_it != elem_end; ++elem_it)
  {
 
    const Elem* el = *elem_it;
   
    if (volume.find(el) == volume.end())
    {	         
      volume[el] = el->volume(); 
    }
  }
  
  
}
//-------------------------------------------------------------------------------
unsigned int KspaceIntegration::how_many_elements_to_do()
{
  unsigned int result = 0;

  MeshBase::const_element_iterator       elem_it  = kmesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end = kmesh->active_elements_end();
  for (; elem_it != elem_end; ++elem_it)
  {
    const Elem* el = *elem_it;
    if ( kspace_integral.find(el) ==  kspace_integral.end() ) result++;
  }

  return result;

}
//---------------------------------------------------------------------------------
double  KspaceIntegration::estimate_error(void)
{

  double result;
  double t1 = 0.0; double t2 = 0.0;


  map<const Elem*, double>::iterator   it;
  map<const Elem*, double>::iterator   it1 = real_space_density.begin();
  map<const Elem*, double>::iterator   it2 = real_space_density.end();

  for (it = it1; it != it2; ++it)
  {
    const Elem* el = it->first;

    double volume;
    
    if (el != NULL) 
      volume =  el->volume();
    else
      volume = 1.0; // bulk calculus 

    t1 += real_space_density[el] * real_space_density[el] * volume;

    t2 += (real_space_density[el] - old_real_space_density[el]) *  (real_space_density[el] - old_real_space_density[el]) * volume;


  }

  result = t2/t1;

  return(result);

}
//=================================================================//
std::vector<double>   KspaceIntegration::get_density_in_k_space(void)  const
{
  
  vector<double> result ;

  MeshBase::const_element_iterator       elem_it  = kmesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end = kmesh->active_elements_end(); 

  unsigned int n_active_elements = 0;

  for ( ; elem_it !=  elem_end ; ++elem_it )
  {
    n_active_elements++;
  }
    
    
 

  result.resize( n_active_elements );
  elem_it  = kmesh->active_elements_begin();
    
  unsigned int j = 0;
  for ( ; elem_it !=  elem_end; ++elem_it )
  {
    const Elem* el = *elem_it;
   

    map <const Elem*, double >::const_iterator it1 = kspace_integral.find(el);
    

    map <const Elem*, double >::const_iterator it2 = volume.find(el);


 

    result[j] = it1->second / it2->second;
    

    j++;
  }


  return(result);
 
}
