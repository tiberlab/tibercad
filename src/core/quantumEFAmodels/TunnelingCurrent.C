#include "QuantumDensity.h"
#include "SimulationEnvironment.h"
#include "TunnelingCurrent.h"
#include "Control.h"
#include "gnuplot_io.h"
#include "VTKIO.h"
#include "mesh_tools.h"

using namespace std;

 

extern "C"
{
#ifdef ENABLE_HETERO
  int call_hetero(double potential, double *kpar, double *transmission, double* Energy, int N);
#else
  int call_hetero(double potential, double *kpar, double *transmission, double* Energy, int N) {};
#endif
}
//===================================================================================//

TunnelingCurrent::TunnelingCurrent()
{
  Vmesh = NULL;

  Ves = NULL;

}


//===================================================================================//
TunnelingCurrent:: ~TunnelingCurrent()

{
  delete Vmesh;
  
  delete  Ves;
  
}



TunnelingCurrent*  TunnelingCurrent::create()
{
#ifndef ENABLE_HETERO
  throw InitFailedException("Cannot create TunnelinCurrent model: Hetero is disabled.");
#endif
  return (new TunnelingCurrent );
}


//===================================================================================//
void TunnelingCurrent::write_current()
{
 std::ofstream current_file;
 current_file.open(opt.filename.c_str());

 if (!current_file.good())
 {
   cerr << "TunnelingCurrent: output file problem\n";
   //error();

 }

 MeshBase::node_iterator       it     = Vmesh->nodes_begin();
 const MeshBase::node_iterator it_end  = Vmesh->nodes_end(); 

 for ( ;  it != it_end ; ++it)
 {
   const Node* nd = *it;
   current_file << setprecision(14) << current[nd] << "\n";

 } 
 
 current_file.close();

}

//=====================================================================================// 
void TunnelingCurrent::read_current()
{
  std::ifstream current_file;
  current_file.open(opt.filename.c_str());

  if (!current_file.good())
  {
    cerr << "TunnelingCurrent: input file problem\n";
    //error();

 }

 MeshBase::node_iterator       it     = Vmesh->nodes_begin();
 const MeshBase::node_iterator it_end  = Vmesh->nodes_end(); 

 string string_from_file;

 for ( ;  it != it_end ; ++it)
 {
   const Node* nd = *it;
   
   getline(current_file, string_from_file );
  
   istringstream input_string(string_from_file);
  
   input_string >> current[nd];
   

 } 
 
 current_file.close();


}
//======================================================================================//

void TunnelingCurrent::do_plot (void)
{
  //IV-output
  const Device& dev = get_environment().get_device();

  string suffix = get_control().get_filename_suffix();
  string outdir = get_control().get_output_dir();
  string format = get_control().get_output_format();

  string suff;
 
  suff = ".dat";
  
  const std::set< std::string >& plotvariables = get_control().get_plotvariables();

  if (plotvariables.find("tunneling_current") != plotvariables.end())
  {
    string filename(outdir + "/" + get_name() +  suffix + suff);

   

    vector<string> names(1,"density[atomic_units]");

    vector<double> results;
    unsigned int el_number = 0; 

    //results.resize(Vmesh->n_elem());

    results.resize(Vmesh->n_nodes());

    map < const Node*, double > :: iterator it1 = current.begin();
    for (; it1 !=  current.end() ;++it1)
    {
      results[el_number] = it1->second;
      el_number++;
      cerr  << (*it1->first)(0)  << "     " << it1->second << "\n";
    }


   
    GnuPlotIO(*Vmesh).write_nodal_data(filename, results, names);



      

  }

  //----------------------------------------------------------------------------


}

void TunnelingCurrent::build_elemental_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend)
{
}



//=================================================================================//
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



  Ves = new EquationSystems(*Vmesh);
  
  Ves->add_system<LinearImplicitSystem> ("current");
        
  Ves->get_system("current").add_variable("I", FIRST);
        
  Ves->init();


 
  

}
//=========================================================================//
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
  opt.energy_int_zero_limit = mod_opt.get_option("energy_integration_zero_limit", 1e-12);
  opt.energy_int_uniform_refinement = mod_opt.get_option("energy_integral_uniform_refinement", false); 

  string suffix = get_control().get_filename_suffix();
  string outdir = get_control().get_output_dir();
 
  opt.filename = outdir + "/" + get_name() + suffix;
  opt.filename = mod_opt.get_option("filename",opt.filename);

  opt.read_results_from_file = mod_opt.get_option("read_results_from_file", false);
  
  opt.write_results_to_file = mod_opt.get_option("write_results_to_file", true); 

  
}
//========================================================================//
double TunnelingCurrent::get_current(double v)
{

  double result = 0;

  vector<Point> v_point(1);

  v_point[0](0) = v; 


  LinearImplicitSystem& system = Ves->get_system<LinearImplicitSystem> ("current");
  
  const DofMap& dof_map = system.get_dof_map();
    
  


  MeshBase::element_iterator       it     = Vmesh->active_elements_begin();
  const MeshBase::element_iterator it_elem_end  = Vmesh->active_elements_end();

  bool not_found = true; 

  for ( ;  ( (it != it_elem_end) && not_found) ; ++it) 
  {
    const Elem* el = *it;
    if (el->contains_point(v_point[0]))
    {
      not_found = false;

      FEType fe_type (FIRST , LAGRANGE);

      AutoPtr<FEBase> fe (FEBase::build(1,    fe_type));
   
      vector<Point> points(1);
      FEInterface::inverse_map(1, fe_type, el, v_point, points);

      const std::vector<std::vector<Real> >& phi = fe->get_phi();

      fe->reinit(el, &points);
      
      for (unsigned int i = 0; i < el->n_nodes(); i++)
      {
	const Node* nd = el->get_node(i);

	result  += phi[i][0] * current[nd];      
      }

    }

  }

  if (not_found) cerr <<"Warning!!!!\n";

  return result;

}


//========================================================================//
void TunnelingCurrent::do_solve()
{

  parse_options();

 

  build_V_grid();

  build_k_grid();


  if ( opt.read_results_from_file )
  {
    read_current();
  }
  else
  {

    real_space_density.clear();


    MeshBase::node_iterator       it_elem     = Vmesh->nodes_begin();
    const MeshBase::node_iterator it_elem_end  = Vmesh->nodes_end(); 
 

 

    for ( ;  it_elem != it_elem_end ; ++it_elem) 
    {






      unsigned int max_refinement = 0;
      
      applied_voltage_node  = *it_elem;

      cerr << "------------------------------\n";
      {  
	cerr << " applied volatage " << (*applied_voltage_node)(0) << "\n";
      }
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


  

      kspace_integral.clear();

      volume.clear();
    

      calculate_volumes();

      calculate_density();
 

      //prepare_system_solution();

  

      double x1;
      double x2;


      {
	map < const Node*, double > :: iterator it1 = current.find(applied_voltage_node);
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

	 
	      
	    mesh_refinement.flag_elements_by_error_fraction (error,KspaceIntegration::opt.refine_fraction,0.0, 10);
 
	    mesh_refinement.refine_and_coarsen_elements();


	    rotate_mesh(kmesh, transform_matrix);

	  
	    eq->reinit();
	  
	    kmesh->print_info();
	  
	    cerr <<   "  we have to do  " <<  how_many_elements_to_do() << "  elements \n";
	  
	    calculate_volumes();
	    calculate_density();

	    {
	      map < const Node*, double > :: iterator it1 = current.find(applied_voltage_node);
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

      {//k-space output
	std::ostringstream o;
	o << (*applied_voltage_node)(0);
	additional_name_suffix = o.str();
	
	//	KspaceIntegration::do_plot();
      }

    
      delete(eq);
  

    }//voltage loop
  

    if (opt.write_results_to_file) write_current();


  }


 
 
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

  current[applied_voltage_node] = 0;

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
	

	f *= get_degeneracy_factor();


	current[applied_voltage_node] += f ;

	kspace_integral[elem] += f;
      
      }
    }
    else
    {
      
      current[applied_voltage_node] += kspace_integral[elem];

    }
    


  }



 


  

}


//================================================================//



//==============================================================//
double TunnelingCurrent::calculte_at_k_point(const Point& k)
{
  double factor = Constants::e * Constants::e /
    (   (Constants::hbar)*( Constants::bohr_radius * Constants::bohr_radius )    ) /1e4;

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




  

  double applied_voltage = (*applied_voltage_node)(0);
 
	   
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
  else if (format == "vtk")
    suff = ".vtk";
  

  
  const std::set< std::string >& plotvariables = get_control().get_plotvariables();

  double voltage = (*applied_voltage_node)(0);
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
    else if (format == "vtk")
      TiberVTKIO(*kmesh).write_elemental_data(filename, results, names);
    else
    {
      cout << "Output format not supported. Falling back to GMV." << endl;
      GMVIO_cell(*kmesh).write_ascii_cell_data(filename, results, names);
    }



 
  }
  
}
//==========================================================================//




 






//========================================================================//

double TunnelingCurrent::integrate_over_energy(double kpar[3], double electric_potential)
{


 
  
  double result = 0.0;
  double result_old = 0.0;

  const double Ec = 1.17;

  //------------------------------------------------------//
  //Equation system
  Emesh = new Mesh(1);

  unsigned int num_nodes = opt.init_nodes_for_energy_int ;  

  double E1;
  double E2;


  opt.Efermi_right = opt.Efermi_left - electric_potential;
  
 
  E1 = Ec;

  
  E2 = opt.Efermi_left;
  

  E2 = opt.Efermi_left + 40.0 * SimulationOptions::temperature * Constants::k_Boltzmann;

  MeshTools::Generation::build_line (*Emesh, 
				     num_nodes, E1, E2, 
				     EDGE2);


  //------------------------------------------------------//
 
  energy_integral.clear();

  //cerr << "------------------------------------------\n";


  result = integrate_over_fix_energy( Emesh, kpar, electric_potential  );
  
  // cerr << result << "   ";

  if (opt.energy_int_refinement)
  {

    bool do_integration;

    if (abs(result) < opt.energy_int_zero_limit)
      do_integration = false;
    else
      do_integration = true;
    

    MeshRefinement mesh_refinement(*Emesh);
    

    for ( ; do_integration ; )
    {

      
      if (opt.energy_int_uniform_refinement)
	mesh_refinement.uniformly_refine(1);
      else
      {

	ErrorVector error = ErrorVector(Emesh->n_elem(), Emesh);
      
	estimate_error_for_energy_refinement(error);

	mesh_refinement.refine_fraction() = 0.2;
	mesh_refinement.coarsen_fraction() = 0.0;
	mesh_refinement.max_h_level() = 20;

	mesh_refinement.flag_elements_by_error_fraction(error);         

      }

      mesh_refinement.refine_and_coarsen_elements();

      result = integrate_over_fix_energy( Emesh,  kpar, electric_potential);

      



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


 


  delete Emesh;



  return(result);

 


}






//==========================================================//
double TunnelingCurrent::integrate_over_fix_energy(const Mesh* Emesh, double kpar[3], double electric_potential)
{

  double result = 0;

  FEType fe_type (FIRST , LAGRANGE);

  AutoPtr<FEBase> fe (FEBase::build(1, fe_type));

  // QGauss qrule (1,  NINTH);


  QGauss qrule (1,  FIFTH);
 
    
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

      if (energy_integral.find(el) ==  energy_integral.end())
      {

	fe->reinit (el);

	for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
	{      
	  energy_values.push_back(q_point[qp](0));
	}

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
      transmission_values[i1] = transmission_array[i1] * 
	thermal_probability( opt.Efermi_left, energy_array[i1] );

   
  }


 

  MeshBase::const_element_iterator       elem_it  = Emesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end = Emesh->active_elements_end(); 
 
  unsigned int point = 0;

  for (; elem_it != elem_end; ++elem_it)
  {
    
    const Elem* el = *elem_it;
    

    if (energy_integral.find(el) ==  energy_integral.end())
    {

      fe->reinit (el);

      energy_integral[el] = 0.0;

      for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
      {
	
	double Energy = q_point[qp](0);    

     
	energy_integral[el] += transmission_values[point] * JxW[qp] ;

	point++;

      }
    }
  }


  elem_it  = Emesh->active_elements_begin();
  for (; elem_it != elem_end; ++elem_it)
  {
    
    const Elem* el = *elem_it;
    result += energy_integral[el];

  }

  

  return result;
}


//===================================================================================//

void TunnelingCurrent::estimate_error_for_energy_refinement(ErrorVector& error)
{
 std::fill (error.begin(), error.end(), 0.0);
  

  MeshBase::const_element_iterator       elem_it1  = Emesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end1 = Emesh->active_elements_end(); 
  for (; elem_it1 != elem_end1; ++elem_it1)
  {

   
    // e is necessarily an active element on the local processor
    const Elem* el = *elem_it1;
    const unsigned int el_id = el->id();


    error[el_id] = abs(energy_integral[el]); 
    
  

  }


}
