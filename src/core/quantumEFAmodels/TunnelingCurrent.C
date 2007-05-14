#include "QuantumDensity.h"
#include "SimulationEnvironment.h"
#include "TunnelingCurrent.h"
#include "Control.h"
#include "gnuplot_io.h"
using namespace std;

 

extern "C"
{
   double call_hetero(double potential, double *kpar);
}
//===================================================================================//

TunnelingCurrent::TunnelingCurrent()
{
  Vmesh = NULL;
}

//============================================//
//============================================//
TunnelingCurrent:: ~TunnelingCurrent()

{
  delete Vmesh;
  
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
    string filename(outdir + "/" + get_name() +
        "_k_space" + suffix + suff);

   

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
   //  else if (format == "gnuplot")
//       GnuPlotIO(*Vmesh).write_nodal_data(filename, results, names);
//     else if (format == "ise")
//       TecplotIO(*Vmesh).write_nodal_data(filename, results, names);
//     else
//     {
//       cout << "Output format not supported. Falling back to GMV." << endl;
//       GMVIO(*Vmesh).write_nodal_data(filename, results, names);
//     }

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


  //build mesh
//  kmesh = new Mesh(k_dim);


  //build applied voltage mesh
  Vmesh = new Mesh(1);

  unsigned int num_nodes = 1;  // 5  2  10;

  Real start;
  Real end;

  start = 1.00;
  end = 2.00;




  ElemType type;



  MeshTools::Generation::build_line (*Vmesh, 
				     num_nodes, start, end, 
				     EDGE2);


  

 
  

}

void TunnelingCurrent::parse_options(void)
{
  KspaceIntegration::parse_options();

}

//========================================================================//
void TunnelingCurrent::do_solve()
{

  parse_options();

  real_space_density.clear();

  build_V_grid();

  build_k_grid();

  eq = new EquationSystems(*kmesh);

  eq->add_system<LinearImplicitSystem> ("k-integration");
  
  system = &(eq->get_system<LinearImplicitSystem> ("k-integration"));
  
  system->add_variable("u", integration_order);

  
  //!system vector that contains charge density in k space from previous iteration
  NumericVector<Number>& old_density = system->add_vector("old density");

  
  eq->init();


  MeshBase::element_iterator       it_elem     = Vmesh->active_elements_begin();
  const MeshBase::element_iterator it_elem_end  = Vmesh->active_elements_end(); 
 

  unsigned int max_refinement = 0;

  for ( ;  it_elem != it_elem_end ; ++it_elem) 
  {


    applied_voltage_elem  = *it_elem;


    k_point_density.clear();
    k_point_charge.clear();

    calculate_at_each_k_point();
 
   
    calculate_density();
 

    prepare_system_solution();

    double x1;
    double x2;


    {
      map < const Elem*, double > :: iterator it1 = real_space_density.find(applied_voltage_elem);
      x2 = it1->second;
    }


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
	  
	  // kmesh->print_info();
	  
	  eq->reinit();
	  
	  calculate_at_each_k_point();
	  
	  calculate_density();

	  {
	    map < const Elem*, double > :: iterator it1 = real_space_density.find(applied_voltage_elem);
	    x1 = it1->second;
	  }

	  cerr << "x1 = " << x1 << "\n";
	  
	  prepare_system_solution();
	  
	  	 
	  
	  system->solution->close();
	  
	 
	  
	  norm_of_error = std::abs(x1/x2 - 1.0);
	  x2 = x1 ;

	  old_density = *(system->solution);

	  std::cout << "k space grid has " << kmesh->n_nodes() << " nodes " << flush;
	  std::cout <<  "quantum density error " << norm_of_error << endl << flush;
	  
	  
	}

	max_refinement++;
      }

      mesh_refinement.uniformly_coarsen(max_refinement);
    }//end of refinement block
    
    
    k_space_output();

  }//voltage loop
  

 
}


//================================================================//

void TunnelingCurrent::calculate_at_each_k_point()
{


  const double factor = -Constants::e * Constants::e /
	 (  M_PI * (Constants::hbar)*( Constants::bohr_radius * Constants::bohr_radius )    ) /1e4; 

  Real  applied_voltage;

  double transm;

  MeshBase::node_iterator       it     = kmesh->active_nodes_begin();
  const MeshBase::node_iterator it_el  = kmesh->active_nodes_end();

  

  for ( ; it != it_el ; ++it) 
    {
      
      const Node*  nd  = *it;

      map<const Node*, map< const Elem*, double> >::iterator it1;

    

      it1 =  k_point_density.find(nd);

      if ( (it1 == k_point_density.end()) )
	{
	  //  vector<double> k_vector(3, 0.0);


           MeshBase::element_iterator       it_elem     = Vmesh->active_elements_begin();
           const MeshBase::element_iterator it_elem_end  = Vmesh->active_elements_end(); 

	   double  k_vector[3];
	   

	   //cerr << (*nd) << "\n";

	  


	   k_vector[0] = (*nd)(0);
	   k_vector[1] = (*nd)(1);
	   k_vector[2] = (*nd)(2);
	  
	 
	   for (short i1 = 0; i1 < 3; i1++)
	   {
	     k_vector[i1] =  k_vector[i1]/(Constants::bohr_radius * 1e9)  ;
	     if (abs(k_vector[i1]) < 1e-5)  k_vector[i1] = 1e-5;
	   }
	  


	   
	   // k_vector[2] = 1e-5;

	    
	    k_vector[0] = std::ceil(k_vector[0] *1e6)/1e6;
	    k_vector[1] = std::ceil(k_vector[1] *1e6)/1e6;
	    k_vector[2] = std::ceil(k_vector[2] *1e6)/1e6;




	   // Point Elem::centroid 
	   Point  x = applied_voltage_elem->centroid();

	   applied_voltage = x(0);
   

	   //cerr << "call hetero " << " k = " << k_vector[0] <<"  " << k_vector[1] << "   " <<k_vector[2] << "\n";  


	   
	   transm =call_hetero(applied_voltage,k_vector);

	   

	   //cerr<<  endl<<endl;
	  

	   // transm = -transm *  Constants::e  / ( (2.0* M_PI)* (2.0* M_PI)*(2.0* M_PI)  
	   //                                          * (Constants::hbar / Constants::e )    ) * 1.0e18 / 1.0e4 ;
	   //
	   

	  

	   transm *= factor;

	   //transm  *= -1.0;

	   // e[C], hbar[eVs]  ! A/nm^2 => A/cm^2
       /*
	   cerr<<  endl<<endl << " *****************************"<<endl;
	   cerr << "  V, current(kpar) = " << applied_voltage  << "        " <<  transm  <<  endl;
	   cerr<<  endl<<endl << " *****************************"<<endl ;
       */
	   

	   std::map<const Elem*, double> transmission_map;

	   transmission_map.insert( pair< const Elem*, double> (applied_voltage_elem, transm) );

	  
       

	   k_point_density.insert( pair< const Node*, map<const Elem*, double> > (nd,transmission_map ) );



	   k_point_charge.insert(pair<const Node*, double > (nd, transm));
	  

	  

	}

    }


#ifdef DEBUG
  cerr << "Schroedinger equation at each point is solved\n";
#endif 
 
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
    string filename(outdir + "/" + get_name() +
		    "_k_space" + suffix +"_voltage_" +v_s.str() + suff );

  

    std::vector<double> results;
    std::vector<std::string> names;
    names.resize(1, "current");
    
    results.resize( kmesh->n_nodes() );
    MeshBase::node_iterator       it     = kmesh->active_nodes_begin();
    const MeshBase::node_iterator it_el  = kmesh->active_nodes_end();
    for ( ; it != it_el ; ++it) 
    {
      const Node* nd  = *it;
      
      double t = k_point_density[nd][applied_voltage_elem];
      
      results[nd->id()] = t;

    
    }

   

  

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
  
}
