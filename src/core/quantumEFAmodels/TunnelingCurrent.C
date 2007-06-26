#include "QuantumDensity.h"
#include "SimulationEnvironment.h"
#include "TunnelingCurrent.h"
#include "Control.h"
#include "gnuplot_io.h"

#include "mesh_tools.h"

using namespace std;

 

extern "C"
{
  int call_hetero(double potential, double *kpar, double *transmission, double* Energy, int N);
}
//===================================================================================//

TunnelingCurrent::TunnelingCurrent()
{
  Vmesh = NULL;

  //  kmesh = NULL;
}


//============================================//
TunnelingCurrent:: ~TunnelingCurrent()

{
  delete Vmesh;
  

  //delete kmesh;
}





//============================================//

void TunnelingCurrent::do_plot (void)
{
  //k-space output
  const Device& dev = get_environment().get_device();

  string suffix = get_control().get_filename_suffix();
  string outdir = get_control().get_output_dir();
  string format = get_control().get_output_format();

  string suff;
 
  suff = ".gmv";
  
  const std::set< std::string >& plotvariables = get_control().get_plotvariables();

  if (plotvariables.find("tunneling_current") != plotvariables.end())
 {
    string filename(outdir + "/" + get_name() +  suffix + suff);

   

    vector<string> names(1,"density[atomic_units]");

    vector<double> results;
    unsigned int el_number = 0; 
    results.resize(Vmesh->n_elem());


    map < const Elem*, double > :: iterator it1 = real_space_density.begin();
    for (; it1 !=  real_space_density.end() ;++it1)
    {
      results[el_number] = it1->second;
      el_number++;
      cerr << it1->second << "\n";
    }


   
      GMVIO_cell(*Vmesh).write_ascii_cell_data(filename, results, names);



      

  }

  //----------------------------------------------------------------------------


}

void TunnelingCurrent::build_elemental_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend)
{
}



//============================================//
void TunnelingCurrent::build_V_grid()
{


  

  //build applied voltage mesh
  Vmesh = new Mesh(1);

  unsigned int num_nodes = opt.voltage_steps;  // 5  2  10;

  Real start;
  Real end;

  start = opt.voltage_start;
  end = opt.voltage_stop;
  



  ElemType type;



  MeshTools::Generation::build_line (*Vmesh, 
				     num_nodes, start, end, 
				     EDGE2);


  

 
  

}

void TunnelingCurrent::parse_options(void)
{

  const ModelOptions& mod_opt = get_options();

  KspaceIntegration::parse_options();

  
  opt.Efermi_left = mod_opt.get_option("Efermi_left", 0.0);
  opt.Efermi_right = mod_opt.get_option("Efermi_right", 0.0);
  opt.voltage_start = mod_opt.get_option("voltage_start", 1.0);
  opt.voltage_stop = mod_opt.get_option("voltage_stop", 3.0);
  opt.voltage_steps = mod_opt.get_option("voltage_steps", 25); 

  opt.energy_int_refinement = mod_opt.get_option("energy_integration_refinement", false);
  opt.init_nodes_for_energy_int =  mod_opt.get_option("initial_nodes_for_energy_integration", 3);
  opt.energy_int_tolerance = mod_opt.get_option("energy_integration_tolerance", 1e-2);
  
}

//========================================================================//
void TunnelingCurrent::do_solve()
{

  parse_options();

  real_space_density.clear();

  build_V_grid();

  build_k_grid();


  


  MeshBase::element_iterator       it_elem     = Vmesh->active_elements_begin();
  const MeshBase::element_iterator it_elem_end  = Vmesh->active_elements_end(); 
 

 

  for ( ;  it_elem != it_elem_end ; ++it_elem) 
  {






    unsigned int max_refinement = 0;

    applied_voltage_elem  = *it_elem;

    cerr << "------------------------------\n";
    cerr << " applied volatage " << applied_voltage_elem->centroid() << "\n";
    cerr << "------------------------------\n ";  


    // kmesh = new Mesh(*(Kspace::kmesh));


    {

      ostringstream temp;
      temp << "mesh_" << max_refinement << ".xda";
      kmesh->write(temp.str()); 

    }

    eq = new EquationSystems(*kmesh);

    eq->add_system<LinearImplicitSystem> ("k-integration");
  
    system = &(eq->get_system<LinearImplicitSystem> ("k-integration"));
  
    system->add_variable("u", integration_order);

  
    //!system vector that contains charge density in k space from previous iteration
    // NumericVector<Number>& old_density = system->add_vector("old density");

  
    eq->init();


    kmesh->print_info();


    k_point_density.clear();
    k_point_charge.clear();


    kspace_integral.clear();
    volume.clear();
   

    calculate_volumes();

    calculate_density();
 

    //prepare_system_solution();

  

    double x1;
    double x2;


    {
      map < const Elem*, double > :: iterator it1 = real_space_density.find(applied_voltage_elem);
      x2 = it1->second;
    }

    cerr << "  x2 = "<< x2 << "\n";

    if (KspaceIntegration::opt.k_domain_refinement) 
    {


      //----------------------------------------
      //refinement block
      //---------------------------------------
     
     

      MeshRefinement mesh_refinement(*kmesh);
      

      

      
      double norm_of_error = KspaceIntegration::opt.relative_accuracy;

      
      
      for ( ; (norm_of_error >=  KspaceIntegration::opt.relative_accuracy) ;  ) 
      {//for
	
	if (KspaceIntegration::opt.uniform_refinement)
	  mesh_refinement.uniformly_refine(1);
	else
	{
	  
	  ErrorVector error = ErrorVector(kmesh->n_elem(), kmesh);
	  
	  
	  // KellyErrorEstimator error_estimator;
	  
	  Tensor2Gen RotM_inv =  transform_matrix.transpose() ;
	  
	  rotate_mesh(kmesh,  RotM_inv );
	  

	  estimate_error_for_refinement(error);

	 
	 
	  //mesh_refinement.flag_elements_by_elem_fraction (error,KspaceIntegration::opt.refine_fraction,0.0, 10);
	      
	  mesh_refinement.flag_elements_by_error_fraction (error,KspaceIntegration::opt.refine_fraction,0.0, 10);
 
	  mesh_refinement.refine_and_coarsen_elements();


	  rotate_mesh(kmesh, transform_matrix);

	  
	  eq->reinit();
	  
	  kmesh->print_info();
	  
	  cerr <<   "  we have to do  " <<  how_many_elements_to_do() << "  elements \n";
	  
	  calculate_volumes();
	  calculate_density();

	  {
	    map < const Elem*, double > :: iterator it1 = real_space_density.find(applied_voltage_elem);
	    x1 = it1->second;
	  }

	  cerr << "x1 = " << x1 << "  x2 = "<< x2 << "\n";
	  
	
	  
	  norm_of_error = std::abs(x1/x2 - 1.0);
	  x2 = x1 ;

	 

	  std::cout << "k space grid has " << kmesh->n_nodes() << " nodes " << flush;
	  std::cout <<  "quantum density error " << norm_of_error << endl << flush;
	  
	  
	}

	max_refinement++;
      }


      k_space_output();


   

      mesh_refinement.uniformly_coarsen(MeshTools::n_levels(*kmesh));

      

     
    }//end of refinement block
    
    


   
   

    
    delete(eq);
  

  }//voltage loop
  

 
}

//================================================================//
void TunnelingCurrent::calculate_density()
{
  const DofMap& dof_map = system->get_dof_map();
    
  FEType fe_type = dof_map.variable_type(0);

 
    
  AutoPtr<FEBase> fe (FEBase::build(k_dim, fe_type));

  QGauss qrule (k_dim, THIRD);
    
  fe->attach_quadrature_rule (&qrule);

 
  
  const std::vector<Real>& JxW = fe->get_JxW();
  
  const std::vector<Point>& q_point = fe->get_xyz();
  
 
  
  QGauss qrule_low (k_dim, THIRD);
 

  QGauss* qrule_used;

 
  MeshBase::element_iterator       it     = kmesh->active_elements_begin();
  const MeshBase::element_iterator it_el  = kmesh->active_elements_end();
  std::vector<unsigned int> dof_indices;

  real_space_density[applied_voltage_elem] = 0;

  for ( ; it != it_el ; ++it) //loop over k space elements
  {

    const Elem* elem = *it;
   
    if ( kspace_integral.find(elem) ==  kspace_integral.end() )
    {
      
      if (elem->level() == 0)
      {

	qrule_used = &qrule;

      }
      else
      {
	qrule_used = &qrule;

      }

      fe->attach_quadrature_rule (qrule_used);

      fe->reinit (elem);

   

      dof_map.dof_indices (elem, dof_indices, 0);
    
  

      for (unsigned int qp=0; qp < qrule_used->n_points(); qp++)
      {//qp
      
	double f = calculte_at_k_point(q_point[qp]) * JxW[qp]; 


	{
	  double factor = 1.0;

	  for (short i = 0; i < k_dim; i++)  factor /= (2.0 * M_PI);

	  f *= factor;
	}
	
	real_space_density[applied_voltage_elem] += f ;

	kspace_integral[elem] += f;
      
      }
    }
    else
    {
      
      real_space_density[applied_voltage_elem] += kspace_integral[elem];

    }
    


  }



 


  

}


//================================================================//



//==============================================================//
double TunnelingCurrent::calculte_at_k_point(const Point& k)
{
  double factor = Constants::e * Constants::e /
    (  M_PI * (Constants::hbar)*( Constants::bohr_radius * Constants::bohr_radius )    ) /1e4;

  //factor = 1.0; //---------------test only----------------


  double transm;


  double  k_vector[3];
	   

  //cerr << (*nd) << "\n";

	  


  k_vector[0] = k(0);
  k_vector[1] = k(1);
  k_vector[2] = k(2);
	  
	 
  for (short i1 = 0; i1 < 3; i1++)
  {
    k_vector[i1] =  k_vector[i1]/(Constants::bohr_radius * 1e9)  ;
    if (abs(k_vector[i1]) < 1e-5)  k_vector[i1] = 1e-5;
  }
	  


  	   
  // k_vector[2] = 1e-5;

	    
  k_vector[0] = std::ceil(k_vector[0] *1e5)/1e5;
  k_vector[1] = std::ceil(k_vector[1] *1e5)/1e5;
  k_vector[2] = std::ceil(k_vector[2] *1e5)/1e5;




  Point  x = applied_voltage_elem->centroid();

  double applied_voltage = x(0);
 
	   
  transm =integrate_over_energy(k_vector, applied_voltage);

  
  transm *= factor;


  return(transm);

}




//==============================================================//
void TunnelingCurrent::do_init()
{
  KspaceIntegration::do_init();
}
//==============================================================//
void TunnelingCurrent::k_space_output(void)
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
  
  

  
  const std::set< std::string >& plotvariables = get_control().get_plotvariables();

  double voltage = (applied_voltage_elem->centroid())(0);
  ostringstream v_s;
  v_s << voltage;
  if (plotvariables.find("k-space_tunneling") != plotvariables.end())
  {

    kmesh->print_info();

    string filename(outdir + "/" + get_name() +
		    "_k_space" + suffix +"_voltage_" +v_s.str() + suff );

  

    std::vector<double> results;
    std::vector<std::string> names;
    names.resize(1, "current");


   


    MeshBase::const_element_iterator       elem_it  = kmesh->active_elements_begin();
    const MeshBase::const_element_iterator elem_end = kmesh->active_elements_end(); 

    unsigned int n_active_elements = 0;

    for ( ; elem_it !=  elem_end ; ++elem_it )
    {
      n_active_elements++;
    }
    
    
    

    results.resize( n_active_elements );
    elem_it  = kmesh->active_elements_begin();
    
    unsigned int j = 0;
    for ( ; elem_it !=  elem_end; ++elem_it )
    {
      const Elem* el = *elem_it;
      results[j] =  kspace_integral[el]/volume[el];
      j++;
    }
    

    if (format == "gmv")
      GMVIO_cell(*kmesh).write_ascii_cell_data(filename, results, names);
    else if (format == "ise")
      TecplotIO_cell(*kmesh).write_cell_data(filename, results, names);
    else
    {
      cout << "Output format not supported. Falling back to GMV." << endl;
      GMVIO_cell(*kmesh).write_ascii_cell_data(filename, results, names);
    }



 
  }
  
}
//==========================================================================//
void TunnelingCurrent::estimate_error_for_refinement(ErrorVector& error)
{
  //------------------------------------------------
  //volume of active alements

 

  //error.resize (kmesh->n_elem());

  std::fill (error.begin(), error.end(), 0.0);

 

  


  //--------------------------------------------------
  
  double mean_value = 0;
  {

    double temp = 0;
    double volume_total = 0;
    
    MeshBase::const_element_iterator       elem_it  = kmesh->active_elements_begin();
    const MeshBase::const_element_iterator elem_end = kmesh->active_elements_end(); 
    for (; elem_it != elem_end; ++elem_it)
    { 
      const Elem* el = *elem_it;

      volume_total += volume[el];
      temp += abs(kspace_integral[el]);
    }
   
    mean_value = temp/volume_total; 
  }


  MeshBase::const_element_iterator       elem_it1  = kmesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end1 = kmesh->active_elements_end(); 
  for (; elem_it1 != elem_end1; ++elem_it1)
  {

   
    // e is necessarily an active element on the local processor
    const Elem* el = *elem_it1;
    const unsigned int el_id = el->id();


    double el_mean_value = kspace_integral[el]/volume[el];

    unsigned int n = el->n_neighbors();
    vector<double> neighbor_values;

    for (unsigned int i = 0; i < n; i++)
    {
      const Elem* el_neighbor = el->neighbor(i);
      if (el_neighbor != NULL)
      {
	if (el_neighbor->active())
	{
	  neighbor_values.push_back(  kspace_integral[el_neighbor]/volume[el_neighbor] );
	}
	else
	{
	  std::vector< const Elem * > active_family;

	  el_neighbor->active_family_tree(active_family);

	  const unsigned int active_family_size =  active_family.size();

	  double t1 = 0;
	  double t2 = 0;
	 
	  for (unsigned int i1 = 0; i1 < active_family_size; i1++)
	  {
	    if (active_family[i1]->is_neighbor(el))
	    {
	 
	      t1 += volume[active_family[i1]];
	      t2 += kspace_integral[ active_family[i1] ] * volume[active_family[i1]];
	    }
	  }
	 
	 
	  neighbor_values.push_back(t2/t1);

	}
      }
      else
      {
	//!here the periodicity (symmetry) of k-space has to be treated somehow....

      }

    }

    double error_cell = 0;

    
   

    for (unsigned int i = 0; i < neighbor_values.size(); i++)
    {
      
      error_cell += abs(el_mean_value - neighbor_values[i])/mean_value;

      error_cell += abs(el_mean_value - neighbor_values[i])/mean_value;

      //  cerr <<  " " << neighbor_values[i] << "    " << mean_value  ;
     

    }

    error[el_id] = kspace_integral[el]; //test
    
    error[el_id] = error_cell;

  }



  // for (unsigned int t = 0; t < error.size(); t++)
  //  cerr << t << "  " << error[t] << "\n";


}
//===================================================================
unsigned int TunnelingCurrent::how_many_elements_to_do()
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
//===================================================================

void TunnelingCurrent::calculate_volumes(void)
{
 
  MeshBase::const_element_iterator       elem_it  = kmesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end = kmesh->active_elements_end(); 
  

 

  for (; elem_it != elem_end; ++elem_it)
  {
 
    const Elem* el = *elem_it;
   
    if (volume.find(el) == volume.end())
    {	      
     
    
      FEType fe_type (FIRST , LAGRANGE);
      
      AutoPtr<FEBase> fe (FEBase::build(k_dim,
                                   fe_type));	      
      QGauss qrule (k_dim, FIRST);

      fe -> attach_quadrature_rule (&qrule);
      
      const std::vector<Real>& JxW = fe->get_JxW();
    
      fe->reinit(el);

      double el_volume = 0.0;
      
      for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
	el_volume += JxW[qp];
    
      volume[el] = el_volume; 
    }
  }
  
  
}

//========================================================================//

double TunnelingCurrent::integrate_over_energy(double kpar[3], double electric_potential)
{


 
  
  double result = 0.0;
  double result_old = 0.0;

  const double Ec = 1.17;

  //------------------------------------------------------//
  //Equation system
  Mesh* Emesh = new Mesh(1);

  unsigned int num_nodes = opt.init_nodes_for_energy_int ;  

  double E1;
  double E2;


  opt.Efermi_right = opt.Efermi_left - electric_potential;
  
 
  E1 = Ec;

  
  E2 = opt.Efermi_left;
  



  MeshTools::Generation::build_line (*Emesh, 
				     num_nodes, E1, E2, 
				     EDGE2);


  //------------------------------------------------------//
 
 

  result = integrate_over_fix_energy( Emesh, kpar, electric_potential  );
  
  cerr << result << "   ";

  if (opt.energy_int_refinement)
  {

    bool do_integration;

    if (abs(result) < 1e-50)
      do_integration = false;
    else
      do_integration = true;
    

    MeshRefinement mesh_refinement(*Emesh);
    

    for ( ; do_integration ; )
    {

      

      mesh_refinement.uniformly_refine(1);

      mesh_refinement.refine_and_coarsen_elements();

      result = integrate_over_fix_energy( Emesh,  kpar, electric_potential);

     

      cerr << result << "   ";

      if ( (abs(result) < 1e-20) || (abs(result_old - result)/abs(result) < opt.energy_int_tolerance) )
      {
	do_integration = false;
      }
      else
      {
	result_old = result;

	do_integration = true;
	
      }
    }  

  
  }


  cerr << "\n";


  delete Emesh;



  return(result);

 


}






//==========================================================//
double TunnelingCurrent::integrate_over_fix_energy(const Mesh* Emesh, double kpar[3], double electric_potential)
{

  double result = 0;

  FEType fe_type (FIRST , LAGRANGE);

  AutoPtr<FEBase> fe (FEBase::build(1, fe_type));

  QGauss qrule (1,  NINTH);

 
    
  fe->attach_quadrature_rule (&qrule);

 
  
  const std::vector<Real>& JxW = fe->get_JxW();
  
  const std::vector<Point>& q_point = fe->get_xyz();



  //------------------------------------------------------//
  vector<double> energy_values;
  vector<double> transmission_values;
  {
    MeshBase::const_element_iterator       elem_it  = Emesh->active_elements_begin();
    const MeshBase::const_element_iterator elem_end = Emesh->active_elements_end(); 
    for (; elem_it != elem_end; ++elem_it)
    {
      
      const Elem* el = *elem_it;
      fe->reinit (el);
      for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
      {      
	energy_values.push_back(q_point[qp](0));
      }
    
    }
  }
  int energy_size = energy_values.size();
  transmission_values.resize(energy_size);

  {

    double energy_array[energy_size];
    double transmission_array[energy_size];
    for (unsigned i1 = 0 ; i1 < energy_size  ; i1++)
    {
      energy_array[i1] = energy_values[i1];
      transmission_array[i1] = 0.0;
    }

   
    int status = call_hetero(electric_potential, kpar, transmission_array, energy_array,  energy_size);
    
    for (unsigned int i1 = 0 ; i1 < energy_size  ; i1++)
      transmission_values[i1] = transmission_array[i1];

  }


 

  MeshBase::const_element_iterator       elem_it  = Emesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end = Emesh->active_elements_end(); 
  unsigned int point = 0;
  for (; elem_it != elem_end; ++elem_it)
  {
    
    const Elem* el = *elem_it;
    fe->reinit (el);
    for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
    {

      double Energy = q_point[qp](0);    

     
      result += transmission_values[point] * JxW[qp];

      point++;
    }
  }
  

  return result;
}


//===================================================================================//

