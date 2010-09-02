// $Id$

#include "OptRecombinSpectrum.h"
#include "OpticsKP.h"
#include "SimulationEnvironment.h"
#include "DataOutput.h"

using namespace std;


void OptRecombinSpectrum::do_plot()
{

  KspaceIntegration::do_plot();

  if (plot_solution("optical_spectrum"))
  {
     string filename(get_name() + "_spectrum" +
         TiberCad::get_filename_suffix());

     string format = get_options().get_option("output_format", "grace");

     DataOutput data_output(*_energy_mesh, format);
     data_output.set_output_directory(get_output_directory());


     string dimension;
     double area_dim_factor = 1;
     if (k_dim == 1)
     {
       dimension = "/cm";
       area_dim_factor  = (Constants::bohr_radius * 1e2);
     }
     else if (k_dim == 2)
     {
       dimension = "/cm^2";
       area_dim_factor  = (Constants::bohr_radius * 1e2) * (Constants::bohr_radius * 1e2);
     }
     else if (k_dim == 3)
     {
       dimension = "/cm^3";
       area_dim_factor  = (Constants::bohr_radius * 1e2) * (Constants::bohr_radius * 1e2) * (Constants::bohr_radius * 1e2);
     }

     vector<string> names(1,"power_density[W/eV" + dimension + "]");

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
         value = res_it->second; //[ 1/a.u._of_time/((bohr_radius)^kdim)]
         value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
         value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
	 value /= area_dim_factor ; //[J/(eV*second)/((cm)^kdim)]



       }
#ifdef DEBUG
       else
         cerr << "WARNING!!!!";
#endif



       results.push_back(value);

     }



     data_output.write_cell_data(filename, results, names);



  }


}

//=====================================================================================//


OptRecombinSpectrum::OptRecombinSpectrum(const ModelOptions& options)
  : KspaceIntegration(options)
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


  _quantum_model_initial_state->solve(); //calculate eigenstates

  _quantum_model_final_state->solve(); //calculate eigenstates

  _optical_model->solve(); //calculate matrix elements of P operator




  std::map<const Elem*, double>& spectrum = density;


  _optical_model->calculate_spectrum(*_energy_mesh, opt.Gamma, opt.polariz,  spectrum );

  //for integrated quantity I take a sum of the map
  integrated_quantity = 0.0;

  std::map<const Elem*, double>::iterator it = spectrum.begin();
  std::map<const Elem*, double>::iterator it1 = spectrum.end();

  for (;it != it1; ++it)
    integrated_quantity +=abs(it->second);


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


  if (!_optical_model->is_initialized())
    _optical_model->init();




 //---------quantum models for initial and final states (from OpticsKP module)---------------------------------
  _quantum_model_initial_state    = _optical_model->get_initial_state_model();

  _quantum_model_final_state    = _optical_model->get_final_state_model();



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


 std::string  job_name = mod_spectrum.get_option("process","recombination");
  if (job_name == "recombination")
    job = RECOMBINATION;
  else if (job_name == "absorption")
    job = ABSORPTION;

  else
    throw InitFailedException( "OptSpectrum: Incorrect process: " + job_name );







}
