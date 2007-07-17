#include "OptRecombinSpectrum.h"
#include "OpticsKP.h"
#include "SimulationEnvironment.h"
#include "Control.h"
#include "gnuplot_io.h"
#include "GraceIO.h"
using namespace std;

void OptRecombinSpectrum::do_plot()
{
 // const Device& dev = get_environment().get_device();

  string suffix = get_control().get_filename_suffix();
  string outdir = get_control().get_output_dir();
  string format = get_control().get_output_format();

  string suff;
  if (format == "gmv")
    suff = ".gmv";
  else if (format == "ise")
    suff = ".plt";
  else if (format == "grace")
    suff = ".dat";

  const std::set< std::string >& plotvariables = get_control().get_plotvariables();
  
  if (plotvariables.find("optical_spectrum") != plotvariables.end())
  {
     string filename(outdir + "/" + get_name() +
        "_spectrum" + suffix + suff);

     vector<string> names(1,"power_density[W/eV]"); 

     vector<double> results;

     MeshBase::const_element_iterator       elem_it  = _energy_mesh->active_elements_begin();
     const MeshBase::const_element_iterator elem_end = _energy_mesh->active_elements_end();

     for(;elem_it != elem_end; ++elem_it)
     {
       const Elem* el = *elem_it;
       double value;
       map<const Elem*, double>::iterator res_it = real_space_density.find(el);
       if (res_it != real_space_density.end())
       {
         value = res_it->second; //[ 1/a.u._of_time]
         value /= Constants::atomic_time; //[1/second]
         value *= Constants::elementary_charge; //[J/(eV*second)]
       }
       else
       {
         cerr << "WARNING!!!!";
       }
       


       results.push_back(value);

     } 

     
     if (format == "gmv")
       GMVIO_cell(*_energy_mesh).write_ascii_cell_data(filename, results, names);
     else if (format == "ise")
       TecplotIO_cell(*_energy_mesh).write_cell_data(filename, results, names);
     else if (format == "grace")
       GraceIO(*_energy_mesh).write_elemental_data(filename, results, names);
     else
     {
       cout << "Output format not supported. Falling back to GMV." << endl;
       GMVIO_cell(*_energy_mesh).write_ascii_cell_data(filename, results, names);
     }


  }


}




OptRecombinSpectrum::OptRecombinSpectrum()
{
  _quantum_model_initial_state = NULL;
  _quantum_model_final_state = NULL ;

 
  _optical_model = NULL;

 
  _energy_mesh = NULL;

}


//===============================================================//

OptRecombinSpectrum::~OptRecombinSpectrum()
{
  delete(_energy_mesh);
}



//================================================================//
void OptRecombinSpectrum::calculate_for_k_point(const Point& k_point, 
                                                std::map<const Elem*, double>& density, 
                                                double& integrated_quantity)
{


  

  vector<double> k_vector(3, 0.0);
  

  cerr << k_point << "\n";

  k_vector[0] = k_point(0);
  k_vector[1] = k_point(1);
  k_vector[2] = k_point(2);


  ModelOptions  quantum_model_opts;


  quantum_model_opts.set_option("k_vector",  k_vector);
  quantum_model_opts["job"] = "eigenstates";

  _quantum_model_initial_state->set_options(quantum_model_opts);

  _quantum_model_final_state->set_options(quantum_model_opts);


  ModelOptions  optical_model_opts;

  optical_model_opts.set_option("k_vector",  k_vector);

  _optical_model->set_options(optical_model_opts);




  _quantum_model_initial_state->solve();

  _quantum_model_final_state->solve();

  _optical_model->solve();

  // get  spectrum data ("density") = map <const Elem*, double>&


  std::map<const Elem*, double>& spectrum = density;


  cerr << opt.polariz(1) << "\n";
  cerr << opt.polariz(2) << "\n"; 
  cerr << opt.polariz(3) << "\n";

  _optical_model->calculate_spectrum(*_energy_mesh, opt.Gamma, opt.polariz,  spectrum );


// integrated_quantity needed for  mesh  refinement, for the  moment is  neglected.
  integrated_quantity = 0.0;

}


//=======================================================================================//
void OptRecombinSpectrum::do_init( )
{

  Kspace::do_init();//--kspace domain--------------


  const ModelOptions& mod_spectrum = get_options();

  //---------------------- options OptSpectrum -----------------------
  //  take  the  name of  Optics  simulation used in  Opt Spectrum
 
  
  std::string optics_simul_name ; 
  if (mod_spectrum.find_option("optical_matr_elem_model"))
  {
    optics_simul_name = mod_spectrum.get_option("optical_matr_elem_model","");
    _optical_model   = dynamic_cast<  OpticsKP* > ( find_simulation(optics_simul_name )   );
    if ( _optical_model== NULL)
      throw  InitFailedException("Optical Spectrum:optical_matr_elem_model  " + 
                                 optics_simul_name + " does not exist");
  }
  else
  {
    throw  InitFailedException("Optical Spectrum: optical_matr_elem_model   has to be specified");
  }

 
  //--------------------------------------------------------------------------------------------//

  


//  const ModelOptions& mod_opt_optics =  _optical_model->see_options();

 //---------quantum models for initial and final states (from OpticsKP module)---------------------------------
  _quantum_model_initial_state    = _optical_model->get_initial_state_model();

  _quantum_model_final_state    = _optical_model->get_final_state_model();
 

 

  // std::string quantum_in_simul_name,  std::string quantum_fin_simul_name;   ;

 //  if  (mod_opt_optics.find_option("initial_state_model"))
//   {
//     std::string quantum_model;
//     quantum_model = mod_opt_optics.get_option("initial_state_model" , "");
//     _quantum_model_initial_state    = dynamic_cast<EnvelopFunctionApprox*> ( find_simulation ( quantum_model));
//     if (  _quantum_model_initial_state == NULL)
//       throw InitFailedException("Optical Spectrum: initial_state_model " + quantum_model + " does not exist\n");
//   }
//   else
//     throw InitFailedException("Optical Spectrum: initial_state_model must be defined\n");
//   //--------------------------

 //  //final state----------------
//   if  (mod_opt_optics.find_option("final_state_model"))
//   {
//     std::string quantum_model;
//     quantum_model = mod_opt_optics.get_option("final_state_model" , "");
//     _quantum_model_final_state    = dynamic_cast<EnvelopFunctionApprox*> ( find_simulation ( quantum_model));
//     if (_quantum_model_final_state == NULL)
//       throw InitFailedException("Optical Spectrum: final_state_model " + quantum_model + " does not exist\n");
//   }
//   else
//     throw InitFailedException("Optical Spectrum: final_state_model must be defined\n");
//   //---------------------------

  if (mod_spectrum.find_option("Emin"))
    opt.Emin = mod_spectrum.get_option("Emin", 0.0);
  else
     throw InitFailedException("Optical Spectrum: Emin must be defined\n");


  
  if (mod_spectrum.find_option("Emax"))
    opt.Emax = mod_spectrum.get_option("Emax", 0.0);
  else
    throw InitFailedException("Optical Spectrum: Emin must be defined\n");



  if (opt.Emax < opt.Emin)  throw InitFailedException("Optical Spectrum: Emax < Emin");

  if (mod_spectrum.find_option("dE"))
    opt.dE = mod_spectrum.get_option("dE", 0.0);
  else
    throw InitFailedException("Optical Spectrum: dE must be defined\n");


  if (opt.dE <= 0)  throw InitFailedException("Optical Spectrum: dE <= 0");
   
  
  
  unsigned int num_nodes = (int)((opt.Emax - opt.Emin)/opt.dE) + 1; 

  _energy_mesh = new Mesh(1);
  
  MeshTools::Generation::build_cube (*_energy_mesh, 
				     num_nodes, 0, 0, 
				     opt.Emin, opt.Emax, 
				     0, 0, 
				     0, 0,
				     EDGE2);




  


}


//============================================//
void OptRecombinSpectrum::parse_options( )
{

  KspaceIntegration::parse_options();

  const ModelOptions& mod_spectrum = get_options();

  opt.Gamma = mod_spectrum.get_option("broadening", 0.007);


  vector<double> polariz;
  mod_spectrum.get_option("polarization", polariz);
  if (polariz.size() == 3)
  {
    opt.polariz(1) = polariz[0];
    opt.polariz(2) = polariz[1];
    opt.polariz(3) = polariz[2];
    
    if (norm (opt.polariz) != 0)
      opt.polariz = opt.polariz/norm( opt.polariz );
    else
      InitFailedException("Optical Spectrum: polarization vector must be non zero");
  }
  else
    InitFailedException("Optical Spectrum: polarization vector must be defined\n");




//Tensor1& polariz

 



}
