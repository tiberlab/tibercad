// $Id: OpticsTB.C 1313 2009-03-04 05:11:04Z alex $

#include "OpticsTB.h"
#include "EigenvalueProblem.h"
#include "Control.h"
#include "DataOutput.h"
#include "SimulationOptions.h"

#include "mesh.h"
#include "elem.h"
#include "mesh_generation.h"

using namespace Constants;

OpticsTB::OpticsTB()

{
  _initial_state_model = NULL;
  _final_state_model = NULL;
  _energy_mesh = NULL;

}


OpticsTB::~OpticsTB()
{

  delete(_energy_mesh);

}

//===============================================//
void OpticsTB::parse_options()
{

  const ModelOptions& mod_opt = get_options();

  {
    std::string quantum_model = mod_opt.get_option("initial_state_model" , "");

    _initial_state_model=dynamic_cast<EigenvalueProblem*>(find_simulation(quantum_model));

    if (_initial_state_model == NULL)
      throw InitFailedException("OpticsTB: initial_state_model is not defined\n");
        
    _initial_state_particle = mod_opt.get_option("initial_state_particle","el");
    
    _initial_eigen_state_numbers.clear();
    
    
    std::vector<int> temp;
    mod_opt.get_option("initial_eigenstates", temp);
    
    if (temp.size() == 0)
      throw InitFailedException("OpticsTB: initial states are not acceptable\n");
    
    _initial_eigen_state_numbers.resize(temp.size());

    for (unsigned int i = 0; i < temp.size(); i++)
    {
      _initial_eigen_state_numbers[i] = temp[i];
    }  
    
  }
  
  {
    std::string quantum_model = mod_opt.get_option("final_state_model" , "");
    
    _final_state_model=dynamic_cast<EigenvalueProblem*>(find_simulation(quantum_model));

    // checks that initial and final state model are the same
    if(_final_state_model!=_initial_state_model)
      throw InitFailedException("OpticsTB: initial and final models differ\n"); 
    
    //check that the final state model has been defined
    if (_final_state_model == NULL)
      throw InitFailedException("OpticsTB: final_state_model is not defined\n");
        
    _final_state_particle = mod_opt.get_option("final_state_particle","hl");
    
    _final_eigen_state_numbers.clear();
      
    std::vector<int> temp;
    mod_opt.get_option("final_eigenstates", temp);
    
    if (temp.size() == 0)
      throw InitFailedException("OpticsTB: initial states are not acceptable\n");
    
    _final_eigen_state_numbers.resize(temp.size());
    for (unsigned int i = 0; i < temp.size(); i++)
    {
      _final_eigen_state_numbers[i] = temp[i];
    }
    
  }
 
}


//===============================================//
void OpticsTB::do_init()
{

  std::cout<<"tb optics init"<<std::endl;
  
  this->parse_options();

  std::cout<<"tb optics initialized"<<std::endl;

}

//===============================================//

void OpticsTB::do_solve()
{
  //check that the states are really available
   
  
  int verbose = SimulationOptions::verbose();

  
  if (verbose > 0)
    std::cout << "calculation of the optical matrix elements..." << std::endl;

  check_states();

  this->calculate_P_matrix_elements();   


  if (verbose > 0)  std::cout << "done\n" << std::flush;

  std::vector<double> is_pop, fs_pop;
  std::vector<double> is_ene, fs_ene;

  _initial_state_model->get_populations(_initial_state_particle, is_pop);
  _final_state_model->get_populations(_final_state_particle, fs_pop);
  _initial_state_model->get_eigenvalues(_initial_state_particle, is_ene);
  _final_state_model->get_eigenvalues(_final_state_particle, fs_ene);

  unsigned int n1 =  _initial_eigen_state_numbers.size();
  unsigned int n2 =  _final_eigen_state_numbers.size();

  // write down matrix elements for debugging
  if (verbose > 0)
  {
 
    std::cout<<"i-f      Ei       Ef    |Px|^2       fi    ff"
	     <<std::endl;
    for (int i = 0; i < n1; i++)
    {
      for (int j = 0; j < n2; j++)
      {

	unsigned int is = _initial_eigen_state_numbers[i];
	unsigned int fs = _final_eigen_state_numbers[j];

	std::cout << is << "-" << fs << ": " << std::flush;
	std::cout << is_ene[is] << ":" << fs_ene[fs] << " " << std::flush;      

	std::cout << std::norm(Px_matrix[i][j]) << "   ";
	std::cout << is_pop[is] << ", " << fs_pop[fs] << std::endl;

      }
    }

    std::cout<<"i-f      Ei       Ef    |Py|^2       fi    ff"
	     <<std::endl; 
    for (int i = 0; i < n1; i++)
    {
      for (int j = 0; j < n2; j++)
      {

	unsigned int is = _initial_eigen_state_numbers[i];
	unsigned int fs = _final_eigen_state_numbers[j];

	std::cout << is << "-" << fs << ": " << std::flush;
	std::cout << is_ene[is] << ":" << fs_ene[fs] << " " << std::flush;      

	std::cout << std::norm(Py_matrix[i][j]) << "   ";
	std::cout << is_pop[is] << ", " << fs_pop[fs] << std::endl;

      }
    }

    std::cout<<"i-f      Ei       Ef    |Pz|^2       fi    ff"
	     <<std::endl;
    for (int i = 0; i < n1; i++)
    {
      for (int j = 0; j < n2; j++)
      {

	unsigned int is = _initial_eigen_state_numbers[i];
	unsigned int fs = _final_eigen_state_numbers[j];

	std::cout << is << "-" << fs << ": " << std::flush;
	std::cout << is_ene[is] << ":" << fs_ene[fs] << " " << std::flush;      

	std::cout << std::norm(Pz_matrix[i][j]) << "   ";
	std::cout << is_pop[is] << ", " << fs_pop[fs] << std::endl;

      }
    }


  }


}
//========================================================================================
void OpticsTB::get_states()
{
  //get states from TB model (but these states are not currently used)
  _i_states = _initial_state_model->get_eigen_solution();
  _f_states = _final_state_model->get_eigen_solution();
}

//========================================================================================
void OpticsTB::check_states()
{
  
  //get_states();

  // Find maximum state index of initial states
  unsigned int i; 
  unsigned int n_i =_initial_eigen_state_numbers.size();
  unsigned int n_i_st =_initial_state_model->get_num_states(_initial_state_particle); 
  bool resized = false;

  for(i=0; i<n_i; i++)
  {
    if(_initial_eigen_state_numbers[i] > n_i_st-1)
    {
      _initial_eigen_state_numbers.erase(_initial_eigen_state_numbers.begin()+i);
      resized = true;

    }
  }

  if(resized)
  {
    std::cout<< "initial states redefined: ";
    for(i=0; i<_initial_eigen_state_numbers.size(); i++)
    {
      std::cout<<_initial_eigen_state_numbers[i];
    }
    std::cout<<std::endl;
  }

  // Find maximum state index of final states
  unsigned int n_f = _final_eigen_state_numbers.size();
  unsigned int n_f_st =_final_state_model->get_num_states(_final_state_particle); 
  resized = false;

  for(i=0; i<n_f; i++)
  {
    if(_final_eigen_state_numbers[i] > n_f_st-1)
    {
      _final_eigen_state_numbers.erase(_final_eigen_state_numbers.begin()+i);
      resized = true;
    }
  }

  if(resized)
  {
    std::cout<< "fianal states redefined: ";
    for(i=0; i<_final_eigen_state_numbers.size(); i++)
    {
      std::cout<<_final_eigen_state_numbers[i];
    }
    std::cout<<std::endl;
  }


  
}
//========================================================================================
void OpticsTB::do_plot()
{
  // get  spectrum calculation  options from  opticsTB model  section
  // for  calculation of  spectrum for k = 0

  const ModelOptions& mod_spectrum = get_options();


  const std::set< std::string >& plotvariables = get_control().get_plotvariables();
  if (plotvariables.find("optical_spectrum_k_0") != plotvariables.end())
  {

    std::cerr << "Computing spectrum in gamma point" << std::endl;

    double Gamma = mod_spectrum.get_option("broadening", 0.007);
    double Emin, Emax, dE;


    if (mod_spectrum.find_option("Emin"))
      Emin = mod_spectrum.get_option("Emin", 0.0);
    else
      throw InitFailedException("OpticsTB: Emin must be defined\n");
    

    if (mod_spectrum.find_option("Emax"))
      Emax = mod_spectrum.get_option("Emax", 0.0);
    else
      throw InitFailedException("OpticsKP: Emax must be defined\n");
    
     
    if (Emax < Emin)  throw InitFailedException("OpticsTB: Emax < Emin");

    
    if (mod_spectrum.find_option("dE"))
      dE = mod_spectrum.get_option("dE", 0.0);
    else
      throw InitFailedException("OpticsTB: dE must be defined\n");
    
    
    if (dE <= 0)  throw InitFailedException("OpticsKP: dE <= 0");
    
    double spin_deg = 2.0;
    //if(_initial_state_model->is_relativistic()) spin_deg = 1.0;

    unsigned int num_nodes = (int)((Emax - Emin)/dE) + 1;
    
    
    // do  energy_mesh
    _energy_mesh = new Mesh(1);

    
    MeshTools::Generation::build_cube (*_energy_mesh,
				       num_nodes, 0, 0,
				       Emin, Emax,
				       0, 0,
				       0, 0,
				       EDGE2);
    
    
    std::map<const Elem*, double> spectrum_x;
    std::map<const Elem*, double> spectrum_y;
    std::map<const Elem*, double> spectrum_z;

    Tensor1 polariz_x(0); polariz_x(1) = 1.0;
    Tensor1 polariz_y(0); polariz_y(2) = 1.0;
    Tensor1 polariz_z(0); polariz_z(3) = 1.0;
  
    
    // calculate_spectrum( *_energy_mesh, Gamma,polariz, spectrum ) ;
    
    calculate_spectrum( *_energy_mesh, Gamma, polariz_x, spectrum_x ) ;
    calculate_spectrum( *_energy_mesh, Gamma, polariz_y, spectrum_y ) ;
    calculate_spectrum( *_energy_mesh, Gamma, polariz_z, spectrum_z ) ;
    
    
    //  Emin max , Gamma, polariz,  =  given in  optics model ???
    
    std::string dimension;
    double area_dim_factor = 1;
    
    
    std::vector<std::string> names(3);
    
    names[0] = "power_density_k0_Px[W/eV]";
    names[1] = "power_density_k0_Py[W/eV]";
    names[2] = "power_density_k0_Pz[W/eV]";
    
    
    std::vector<double> results;
    
    {
      MeshBase::const_element_iterator       elem_it  = 
	_energy_mesh->active_elements_begin();
      const MeshBase::const_element_iterator elem_end = 
	_energy_mesh->active_elements_end();
      int n = 0;
      for(;elem_it != elem_end; ++elem_it)
	n++;
      
      results.reserve(n*3);
    }
    

    MeshBase::const_element_iterator       elem_it  = 
      _energy_mesh->active_elements_begin();
    const MeshBase::const_element_iterator elem_end = 
      _energy_mesh->active_elements_end();
    int point = 0;
  
    for(;elem_it != elem_end; ++elem_it)
    {
      const Elem* el = *elem_it;
      double value;
      
      
      //--x - polarization
      value = spin_deg * spectrum_x[el];
      value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
      value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
      value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
      results[3*point + 0] = value;
      
      //--y - polarization
      value =  spin_deg * spectrum_y[el];
      value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
      value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
      value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
      results[3*point + 1] = value;
      
      //--z - polarization
      value =  spin_deg * spectrum_z[el];
      value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
      value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
      value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
      results[3*point + 2] = value;
      
      
      point++;
    }  
    
    std::string filename(get_name() +
			 "_spectrum_k_0" + get_control().get_filename_suffix());
    
    std::string format = get_options().get_option("output_format", "grace");
    
    DataOutput data_output(*_energy_mesh, format);
    data_output.set_output_directory(get_control().get_output_dir());
    
    data_output.write_cell_data(filename, results, names);

  }

}
//========================================================================================

void OpticsTB::calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz,
				  std::map<const Elem*, double>& spectrum )
{


  spectrum.clear();


  std::vector<double> fs_eigen_values;
  std::vector<double> is_eigen_values;

  std::vector<double> fs_populations;
  std::vector<double> is_populations;


  double     trans_energy, f1, f2;


  unsigned int n1 =  _initial_eigen_state_numbers.size();
  unsigned int n2 =  _final_eigen_state_numbers.size();


  _initial_state_model->get_eigenvalues(_initial_state_particle, is_eigen_values);  

  _final_state_model->get_eigenvalues(_final_state_particle, fs_eigen_values);


  _initial_state_model->get_populations(_initial_state_particle, is_populations);

  _final_state_model->get_populations(_final_state_particle, fs_populations);

  
  // loop on  eigenstates

  for (unsigned int i = 0; i < n1; i++)  // "upper" states
  {
      
    for (unsigned int j = 0; j < n2; j++)  // "lower" states
    {

      trans_energy =  is_eigen_values[_initial_eigen_state_numbers[i]] 
	              - fs_eigen_values[ _final_eigen_state_numbers[j]];


      f1 = is_populations[_initial_eigen_state_numbers[i]];   // occupation for  electron

      f2 = fs_populations[_final_eigen_state_numbers[j]]; // occupation for  holes


      std::complex<double> Me = Px_matrix[i][j] * polariz(1) 
	+ Py_matrix[i][j] * polariz(2) 
	+ Pz_matrix[i][j] * polariz(3);
      

      MeshBase::const_element_iterator       el     = Energy.active_elements_begin();
      const MeshBase::const_element_iterator end_el = Energy.active_elements_end();

      for ( ; el!= end_el ; ++el)
      {

        const Elem* elem = *el;

        double En = elem->centroid()(0);

        double Lorenzian =  0.5*Gamma/( ( trans_energy - En) *  ( trans_energy - En)  
					+ (0.5*Gamma)*(0.5*Gamma)) * Hartree;

        double c = 1.0/fine_structure_constant;

        double omega = trans_energy/Hartree;

        //spectrum[elem] += 1 / (2 * M_PI * M_PI) * (omega * omega) /(c*c*c) 
        spectrum[elem] += 1/(2 * M_PI) * (omega * omega) /(c*c*c) 
	                  * Lorenzian * std::norm(Me) * f1 * f2;


      }


    }


  }



}

//=========================================================================================

void OpticsTB::calculate_P_matrix_elements( )
{
  
  ModelOptions options;
  unsigned int n_i =  _initial_eigen_state_numbers.size();
  unsigned int n_f =  _final_eigen_state_numbers.size();

  // assemble matrix for polarization x
  options.set_option("P_matrix", true);
  options.set_option("poldir", 1);
  
  std::cout << "Assemble Px" << std::endl;

  _initial_state_model->assemble(options);

  std::cout << "Compute matrix elements" << std::endl;
  
  // compute matrix elements
  Px_matrix.clear();
  Px_matrix.resize(n_i);
  for (unsigned int j = 0; j < n_i; j++)  Px_matrix[j].resize(n_f);

  for (unsigned int i1 = 0; i1 < n_i; i1++)
  {
    for (unsigned int i2 = 0; i2 < n_f; i2++)
    {

      unsigned int is = _initial_eigen_state_numbers[i1];
      unsigned int fs = _final_eigen_state_numbers[i2];

      Px_matrix[i1][i2] = _initial_state_model->
	calculate_matrix_element(_initial_state_particle, is, 
				 _final_state_particle, fs);

      // convert matrix element from eV*Ang to atomic units
      Px_matrix[i1][i2] *= Hartree/(bohr_radius*1e10);

    }
  }

  std::cout << "Assemble Py" << std::endl;
  // assemble matrix for polarization y
  options.set_option("P_matrix", true);
  options.set_option("poldir", 2);
  
  _initial_state_model->assemble(options);

  // compute matrix elements
  std::cout << "Compute matrix elements" << std::endl;

  Py_matrix.clear();
  Py_matrix.resize(n_i);
  for (unsigned int j = 0; j < n_i; j++)  Py_matrix[j].resize(n_f);  

  for (unsigned int i1 = 0; i1 < n_i; i1++)
  {
    for (unsigned int i2 = 0; i2 < n_f; i2++)
    {

      unsigned int is = _initial_eigen_state_numbers[i1];
      unsigned int fs = _final_eigen_state_numbers[i2];

      Py_matrix[i1][i2] = _initial_state_model->
	calculate_matrix_element(_initial_state_particle, is, 
				 _final_state_particle, fs);

      // convert matrix element from eV*Ang to atomic units
      Py_matrix[i1][i2] *= Hartree/(bohr_radius*1e10);

    }
  }

  std::cout << "Assemble Pz" << std::endl;

  // assemble matrix for polarization z
  options.set_option("P_matrix", true);
  options.set_option("poldir", 3);
  
  _initial_state_model->assemble(options);

  std::cout << "Compute matrix elements" << std::endl;

  // compute matrix elements
  Pz_matrix.clear();
  Pz_matrix.resize(n_i);
  for (unsigned int j = 0; j < n_i; j++)  Pz_matrix[j].resize(n_f);

  for (unsigned int i1 = 0; i1 < n_i; i1++)
  {
    for (unsigned int i2 = 0; i2 < n_f; i2++)
    {

      unsigned int is = _initial_eigen_state_numbers[i1];
      unsigned int fs = _final_eigen_state_numbers[i2];

      Pz_matrix[i1][i2] = _initial_state_model->
	calculate_matrix_element(_initial_state_particle, is, 
				 _final_state_particle, fs);

      // convert matrix element from eV*Ang to atomic units
      Pz_matrix[i1][i2] *= Hartree/(bohr_radius*1e10);
      
    }
  }

}

