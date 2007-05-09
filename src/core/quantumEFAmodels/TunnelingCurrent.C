#include "QuantumDensity.h"
#include "SimulationEnvironment.h"
#include "TunnelingCurrent.h"
#include "Control.h"
#include "gnuplot_io.h"
using namespace std;

extern "C" 
{ 

  double call_hetero_(double& potential, double *kpar);
};

//===================================================================================//

TunnelingCurrent::TunnelingCurrent()
{
  quantum_model = NULL;

  system = NULL;

  eq = NULL;

  kmesh = NULL;
}

//============================================//
//============================================//
TunnelingCurrent:: ~TunnelingCurrent()

{
  delete eq;

  delete kmesh;
  
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
    map < const Elem*, double > :: iterator it1 = transmission_map.begin();
    for (; it1 != transmission_map.end() ;++it1)
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

//  if (k_dim == 1)
//    {
      type = EDGE2;
//     }

//   if (k_dim == 2)
//     {
//       type = QUAD8;
//     }

//   if (k_dim == 3)
//     {
//       type = HEX27;
//     }


// void MeshTools::Generation::build_line (Mesh& *Vmesh,
//                                          const unsigned int nx,
//                                         const Real xmin, const Real xmax,
//                                        const ElemType type,
//                                      const bool gauss_lobatto_grid)


  MeshTools::Generation::build_line (*Vmesh, 
				     num_nodes, start, end, 
				     type);


  
//   rotate_mesh(kmesh, transform_matrix);

//   kmesh->print_info();

 
  

}



//========================================================================//
//================================================================//

void TunnelingCurrent::calculate_at_each_k_point()
{

  Real  applied_voltage;

  double transm;

   MeshBase::node_iterator       it     = kmesh->active_nodes_begin();
   const MeshBase::node_iterator it_el  = kmesh->active_nodes_end();

  
   build_V_grid();

  

   k_point_density.clear();

   cerr <<"eeeeeeeeeeeeee   " << kmesh->n_nodes() << endl;

   for ( ; it != it_el ; ++it) 
    {
      
      const Node*  nd  = *it;

      map<const Node*, map< const Elem*, double> >::iterator it1;

      cerr << "test  " << *nd << "\n";

      it1 =  k_point_density.find(nd);

      if ( (it1 == k_point_density.end()) )

	{
	//  vector<double> k_vector(3, 0.0);


           MeshBase::element_iterator       it_elem     = Vmesh->active_elements_begin();
           const MeshBase::element_iterator it_elem_end  = Vmesh->active_elements_end(); 

          double  k_vector[3];

	  k_vector[0] = (*nd)(0);
	  k_vector[1] = (*nd)(1);
	  k_vector[2] = (*nd)(2);

//Constants::bohr_radius

          k_vector[0] = Constants::bohr_radius * 1e9 *  k_vector[0] ;
          if (k_vector[0] == 0.0)  k_vector[0] = 1e-5;
          k_vector[1] = Constants::bohr_radius * 1e9 *  k_vector[1] ;
          if (k_vector[1] == 0.0)  k_vector[1] = 1e-5;
          k_vector[2] = Constants::bohr_radius * 1e9 *  k_vector[2] ;
          if (k_vector[2] == 0.0)  k_vector[2] = 1e-5;

          cerr << " k  in  a.u. ********   " << ( 11.571 / (Constants::bohr_radius * 1e9) ) << endl<< endl;
    


	 //  ModelOptions quantum_model_opts;
	  

// 	  quantum_model_opts.set_option("k_vector",  k_vector);
// 	  quantum_model_opts.set_option("initial_eigenstates_number",opt.intial_eigenstates_number ); 
// 	  quantum_model_opts["job"] = "density";

 
	  
// 	  quantum_model->set_options(quantum_model_opts);


// 	  quantum_model->solve();

//  call hetero external  for   calculation of T
//  put T_E(kpar) in map <Vmesh_elem, Transm> 

 

        


          for ( ;  it_elem != it_elem_end ; ++it_elem) 
          {

            const Elem*  lin_elem  = *it_elem;



            // Point Elem::centroid 
            Point  x = lin_elem->centroid();

            applied_voltage = x(0);
            //  transm = call_hetero(applied_voltage)
 //result = call_hetero_(potential, kpar, fermi_level);

            cerr << "call hetero " << " k = " << (*nd) ; 
            transm =call_hetero_(applied_voltage,k_vector);
            cerr<<  endl<<endl;
           //  cerr << " transm, V = " << transm << "  " <<  applied_voltage <<  endl;
//             cerr<<  endl<<endl;

// tunn_current_test = -T_integrated_on_E  * el_charge / ( (2.d0*Pi)**3.d0 * hbar ) *&  ! e[C], hbar[eVs]
//            1.d18 / 1.d4; ! A/nm^2 => A/cm^2

            transm = -transm *  Constants::e  / ( (2.0* M_PI)* (2.0* M_PI)*(2.0* M_PI)  
                                                  * (Constants::hbar / Constants::e )    ) * 1.0e18 / 1.0e4 ;


            // e[C], hbar[eVs]  ! A/nm^2 => A/cm^2
       
            cerr<<  endl<<endl << " *****************************"<<endl;
            cerr << "  V, current(kpar) = " << applied_voltage  << "  " <<  transm  <<  endl;
            cerr<<  endl<<endl << " *****************************"<<endl ;


            transmission_map.insert( pair< const Elem*, double> (lin_elem, transm) );


          }

          k_point_density.insert( pair< const Node*, map<const Elem*, double> > (nd,transmission_map ) );




	  

	  

	}

    }


#ifdef DEBUG
  cerr << "Schroedinger equation at each point is solved\n";
#endif 
 
}


//==============================================================//
