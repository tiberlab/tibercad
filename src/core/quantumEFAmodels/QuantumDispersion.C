#include "QuantumDispersion.h"
#include "SimulationEnvironment.h"
#include "Control.h"
#include "gnuplot_io.h"
#include "GraceIO.h"
using namespace std;

//---------------------------------------------------------------//

QuantumDispersion::QuantumDispersion()
{
  quantum_model = NULL;
}

//---------------------------------------------------------------//

QuantumDispersion::~QuantumDispersion()
{
 
}

//---------------------------------------------------------------//

void QuantumDispersion::do_init(void)
{
  const ModelOptions& mod_opt = get_options();

  Kspace::do_init();//--kspace domain--------------

  //---------quantum model---------------------------------------------------------------------//
  
  std::string quantum_simul_name;
  if (mod_opt.find_option("quantum_simulation"))
  {
    quantum_simul_name = mod_opt.get_option("quantum_simulation","");
    quantum_model = dynamic_cast<  EnvelopFunctionApprox* > ( find_simulation( quantum_simul_name  )   );
    if (quantum_model == NULL)
      throw  InitFailedException("QuantumDispersion: quantum_simulation " + quantum_simul_name + " does not exist");
  }
  else
  {
    throw  InitFailedException("QuantumDispersion: quantum_simulation  has to be specified");
  }

 
  //--------------------------------------------------------------------------------------------//
}

//---------------------------------------------------------------//
void 	QuantumDispersion::do_solve (void)
{

  parse_options();

  calculate_eigen_energy();

}

//---------------------------------------------------------------//
void 	QuantumDispersion::parse_options (void)
{

  const ModelOptions& mod_opt = get_options();
  opt.min_eigenvalue_number = mod_opt.get_option("min_eigenvalue_number", 0);
  opt.max_eigenvalue_number = mod_opt.get_option("max_eigenvalue_number", 10);
  opt.bulk_calculation = mod_opt.get_option("bulk_calculation", false);
  

}
//---------------------------------------------------------------//
  
void 	QuantumDispersion::do_plot (void)
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
   else if (format == "grace")
    suff = ".dat";

  const std::set< std::string >& plotvariables = get_control().get_plotvariables();

 if (plotvariables.find("k-space_dispersion") != plotvariables.end())
  {
    string filename(outdir + "/" + get_name() +
		    "_k_space" + suffix + suff);

    
   
    std::vector<double> results;
    std::vector<std::string> names;

    unsigned int number_of_eigs_to_store = opt.max_eigenvalue_number - opt.min_eigenvalue_number + 1;
    names.resize(number_of_eigs_to_store);

    unsigned int number_of_k_points = kmesh->n_nodes();
    results.resize( number_of_eigs_to_store * number_of_k_points );

   
   

    for (unsigned int i = 0; i < number_of_eigs_to_store ; i++)
    {
      std::ostringstream i_str;
      i_str << "state_number_" << i + opt.min_eigenvalue_number ; //The states are numbered starting from 0 
      names[i] = i_str.str();

      for (unsigned int j = 0; j < number_of_k_points ; j++)
	results[number_of_eigs_to_store * j + i] = eigen_energy[j][i];
    }
    


    
    

    if (format == "gmv")
      GMVIO(get_k_mesh()).write_nodal_data(filename, results, names);
    else if (format == "gnuplot")
      GnuPlotIO(get_k_mesh()).write_nodal_data(filename, results, names);
    else if (format == "ise")
      TecplotIO(get_k_mesh()).write_nodal_data(filename, results, names);
    else if (format == "grace")
    {
      //1D only
      Mesh* kmesh1D = new Mesh(get_k_mesh());
      Tensor2Gen inv_transform_matrix = transform_matrix.transpose();
      rotate_mesh (kmesh1D, inv_transform_matrix);

      GraceIO(*kmesh1D).write_nodal_data(filename, results, names);

      delete kmesh1D;
    }
    else
    {
      cout << "Output format not supported. Falling back to GMV." << endl;
      GMVIO(get_k_mesh()).write_nodal_data(filename, results, names);
    } 

  }
}

//--------------------------------------------------------------------------//
void QuantumDispersion::calculate_eigen_energy()
{

  build_k_grid();
  

  unsigned int number_of_k_points = kmesh->n_nodes();
  unsigned int number_of_eigs_to_store = opt.max_eigenvalue_number - opt.min_eigenvalue_number + 1;

  vector<double> temp (number_of_eigs_to_store, 0.0);

  eigen_energy.resize(number_of_k_points, temp);


  for (unsigned int i = 0; i < number_of_k_points; i++)
  {
   
   
    const Node  k_point = kmesh->node(i);

   
    vector<double> k_vector(3, 0.0);

    k_vector[0] = k_point(0);
    k_vector[1] = k_point(1);
    k_vector[2] = k_point(2);

    ModelOptions quantum_model_opts;

    quantum_model_opts.set_option("k_vector",  k_vector);

    quantum_model_opts.set_option("number_of_eigenstates",opt.max_eigenvalue_number + 1); 
  
    if (opt.bulk_calculation)
      quantum_model_opts["job"] = "bulk";
    else
      quantum_model_opts["job"] = "eigenstates";


    quantum_model->set_options(quantum_model_opts);


    quantum_model->solve();

    const std::vector<EnvelopFunctionApprox::eigen_propblem_solution>& solution = quantum_model->get_solution();

    
    for (unsigned int j = 0 ; j < number_of_eigs_to_store ; j++)
    {
      eigen_energy[i][j] = solution[j + opt.min_eigenvalue_number].eigen_energy;
    }

    
  }

}
