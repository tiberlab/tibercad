#include "EigenvalueProblem.h"
#include "Constants.h"
#include "Messages.h"
#include "DataOutput.h"
#include<fstream>

void EigenvalueProblem::init_kspace(void)
{
   do_dispersion=false;

   //! Is used to parse dispersion options
   if(get_options().has_submodel("Dispersion"))
   { 
     ModelOptions::submodel_iterator it(get_options().submodels_begin("Dispersion"));

     const ModelOptions& opts = it->second;
 
     //disp_range[0] = opts.get_option("min_eigenvalue_number", 0);
     //disp_range[1] = opts.get_option("max_eigenvalue_number", 100);

     const ModelOptions kopts = parse_kspace_options(opts);
     
     _kspace = new Kspace(kopts);
     
     if(_kspace==NULL)
       throw InitFailedException("Could not initialize k-space");
     else
       Messages::info("k-space initialized");

     do_dispersion=true;
   }
}
 
ModelOptions EigenvalueProblem::parse_kspace_options(const ModelOptions& opts)
{
  ModelOptions kopts;

  kopts.set_option("mesh_units",get_mesh_units());
  
  unsigned int k_dim = 3 - get_mesh().mesh_dimension();

  kopts.set_option("k_space_dimension",k_dim);

  std::vector<unsigned int>  num_nodes;

  if (opts.find_option("k-path"))
  {

    std::string kpath = opts.get_option("k-path","");
    kopts.set_option("k-path",kpath);
    num_nodes.push_back(10);
    kopts.set_option("number_of_nodes",num_nodes);
    ModelOptions newopts;
    newopts.set_option("output_format","grace");
    set_options( newopts );
  }
 
  if (opts.find_option("number_of_nodes"))
  {
     opts.get_option("number_of_nodes",num_nodes);
     kopts.set_option("number_of_nodes", num_nodes);
  }
  else if (opts.find_option("number_of_elem"))
  {
     opts.get_option("number_of_elem",num_nodes);
     for(int i=0; i< num_nodes.size(); i++)
           if(num_nodes[i]>0) ++num_nodes[i];

     kopts.set_option("number_of_nodes", num_nodes);
  }

  if (opts.find_option("wedge"))
    kopts.set_option("wedge", opts.get_option("wedge",""));

  
  kopts.set_option("k_space_basis", opts.get_option("k_space_basis",true));


  double k_max = opts.get_option("k_max",0.1);

  kopts.set_option("k_max",k_max);


  std::vector<double> k_vector(3,0.0);   
  k_vector[0]=0.0;     k_vector[1]=0.0;     k_vector[2]=k_max; 

  opts.get_option("k1", k_vector);
  kopts.set_option("k1",k_vector);

  k_vector[0]=0.0;     k_vector[1]=k_max;     k_vector[2]=0.0; 
  opts.get_option("k2", k_vector);
  kopts.set_option("k2",k_vector);

  k_vector[0]=k_max;     k_vector[1]=0.0;     k_vector[2]=0.0; 
  opts.get_option("k3", k_vector);
  kopts.set_option("k3",k_vector);  
    

  kopts.set_option("mesh_order",opts.get_option("mesh_order","first") );

  return kopts;
}
  

void EigenvalueProblem::compute_dispersion(void)
{
  if (!do_dispersion) return;
 
  std::cout<<"Compute Dispersion ..." << std::endl;  
  const Mesh* kmesh = _kspace->get_k_mesh(); 	
  unsigned int number_of_k_points = kmesh->n_nodes();

  //if (disp_range[0] < 0) disp_range[0]=0;
  unsigned int number_of_eigs;

  std::cout<<"(EP)number of kp  " << number_of_k_points << std::endl;  

  {
    unsigned int i = 0;
    const Point&  k_point = kmesh->point(i);  

    solve_for_kpoint(k_point);
    number_of_eigs = get_num_states();
        
    std::vector<double> temp(number_of_eigs);
    _dispersion.resize(number_of_k_points, temp);  
    
    for (unsigned int j = 0 ; j <  _dispersion[0].size(); j++)
      _dispersion[0][j] = _solution[j].eigen_energy;
    
  }

  for (unsigned int i = 1; i < number_of_k_points; i++)
  {

    const Point&  k_point = kmesh->point(i);

    solve_for_kpoint(k_point);
    number_of_eigs = get_num_states();
    
    for (unsigned int j = 0 ; j < _dispersion[i].size() ; j++)
      _dispersion[i][j] = _solution[j].eigen_energy;
    

  }
	
}

void
EigenvalueProblem::plot_dispersion(void)
{

  std::vector<std::string> formats;
  get_output_format(formats);


  const Mesh* kmesh = _kspace->get_k_mesh();
  short kdim =  _kspace->mesh_dimension();

    std::cout<<"(EP) kdim: "<< kdim  << std::endl;
   

  for(short k=0; k<formats.size();k++)
  {

    std::string format = formats[k];
  
    if ((format == "grace") && (kdim > 1)) format="vtk";
 
    std::cout<<"(EP) format: "<< format  << std::endl;
   
    std::vector<double> results;
    std::vector<std::string> names;    
    
    unsigned int number_of_eigs = _dispersion[0].size();
    names.resize(number_of_eigs);

    unsigned int number_of_k_points = kmesh->n_nodes();
    results.resize( number_of_eigs * number_of_k_points );


    for (unsigned int i = 0; i < number_of_eigs ; i++)
    {
      std::ostringstream i_str;
      //The states are numbered starting from 0
      i_str << "state_number_" << i;
      names[i] = i_str.str();

      for (unsigned int j = 0; j < number_of_k_points ; j++)
        results[number_of_eigs * j + i] = _dispersion[j][i];
    }


    std::string filename(get_name() + "_dispersion"); 
                         //+ TiberCad::get_filename_suffix());

    DataOutput data_output(*kmesh, format);
    data_output.set_output_directory(get_output_directory());
    //data_output.set_filename(filename);
    
    data_output.write_nodal_data(filename, results, names);
    
  }

}

void EigenvalueProblem::do_plot(void)
{
  SimulationInterface::do_plot();

  if (do_dispersion)
  {  
      compute_dispersion();	  
      plot_dispersion();
  }
}
 
void EigenvalueProblem::solve_for_kpoint(const Point& kpoint)
{
  do_solve_for_kpoint(kpoint);
}
 
 
void EigenvalueProblem::get_eigenvalues(const std::string& particle, 
					std::vector<double>& values) const
{

  unsigned int n = _solution.size();
  unsigned int num_st = 0;
  values.reserve(n);

  for (unsigned int i = 0; i < n; i++)
  {
    if(_solution[i].particle == particle)
    {  
      num_st++;

      values.push_back( _solution[i].eigen_energy ); 
    }
  }  
 
  values.resize(num_st);

}

unsigned int EigenvalueProblem::get_num_states(void) const
{
  return _solution.size();
}


unsigned int EigenvalueProblem::get_num_states(const std::string& particle) const
{
  unsigned int num_i_states = 0;
  for(unsigned int i=0; i<_solution.size(); i++)
  {
    if(_solution[i].particle == particle) num_i_states++;  
  }
  
  return num_i_states;
}

std::vector<unsigned int>
EigenvalueProblem::get_state_indices(const std::string& particle) const
{
  unsigned int num = get_num_states(particle);	
  std::vector<unsigned int> result(num, 0);

  unsigned int num_st=0;
  for(unsigned int i=0; i<_solution.size(); i++)
  {
    if(_solution[i].particle == particle)
      { result[num_st]=i; num_st++; }
  }
  
  return result;
}


void EigenvalueProblem::get_populations(const std::string& particle, 
					std::vector<double>& values) const
{
 
  unsigned int n = _solution.size();
  unsigned int num_st = 0;
  values.reserve(n);
 
  for (unsigned int i = 0; i < n; i++)
  {
    if(_solution[i].particle == particle)
    {  
      num_st++;

      if(_solution[i].statistics == "Fermi")
      {      
	double val = Fermi(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
		       _solution[i].temperature);

	if(particle == "el" || particle == "electron")
	{
	  values.push_back(val);	  
	}	
	
	if(particle == "hl" || particle == "hole")
	{
	  values.push_back(1-val);	  
	}

      }
      else
      {
	double val = Bose(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
		       _solution[i].temperature);

	values.push_back(val);	
	
      }

    }
     
  }

  values.resize(num_st);
 
} 
 
double  EigenvalueProblem::get_population(int i) const
{
  
  if(_solution[i].statistics == "Fermi")
  {        
    double val = Fermi(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
			 _solution[i].temperature);
    
      if(_solution[i].particle == "el" || _solution[i].particle == "electron")
      {
	  return val;	  
      }	
      
      if(_solution[i].particle == "hl" || _solution[i].particle == "hole")
      {	
	  return 1-val;	  
      }
      
  }
  else
  {
    double val = Bose(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
		      _solution[i].temperature);
    
    return val;	
      
  }

}

double  EigenvalueProblem::Fermi(double Energy, double Fermi_energy, double Temperature) const
{
  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - Fermi_energy)/T_EV;
  
  double el_fermi;
  
  if (exp_arg > 35) 
    el_fermi = 0.0;
  else
    el_fermi = 1.0/(std::exp(exp_arg) + 1.0);
  
  return el_fermi;

}

double  EigenvalueProblem::Bose(double Energy, double electro_chem_pot, double Temperature) const
{
  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - electro_chem_pot)/T_EV;
  
  double bose;
  
  if (exp_arg > 35) 
    bose = 0.0;
  else
    bose = 1.0/(std::exp(exp_arg) - 1.0);
  
  return bose;

}


void EigenvalueProblem::write_states(void) const
{

  int num_st=_solution.size();

  std::cout<<std::endl;
  std::cout<<  "# T   level    stat.     pot.       pop."<<std::endl;


  for(int i=0; i< num_st; i++)
  {
    std::cout<<i<<" "<<_solution[i].particle<<" "<< std::setprecision(6)
	     <<_solution[i].eigen_energy<<" "<<_solution[i].statistics
	     <<" "<<std::setw(10)<<_solution[i].electro_chem_pot
	     <<" "<<std::setw(10)<<get_population(i)<<std::endl;
  }
  std::cout<<std::endl;
}

void EigenvalueProblem::write_states(const std::string& filename) const
{

  int num_st=_solution.size();
  std::ofstream file;
  file.open(filename.c_str());
  
  file << "# T   level    stat.     pot.       pop."<<std::endl;

  for(int i=0; i< num_st; i++)
  {
    file <<i<<" "<<_solution[i].particle<<" "<< std::setprecision(6)
	 <<_solution[i].eigen_energy<<" "<<_solution[i].statistics
	 <<" "<<std::setw(10)<<_solution[i].electro_chem_pot
	 <<" "<<std::setw(10)<<get_population(i)<<std::endl;
  }

  file.close();
}

  
void EigenvalueProblem::copy_H_to_solver( )
{
  do_copy_H_to_solver();
}

void EigenvalueProblem::copy_S_to_solver( )
{
  do_copy_S_to_solver();
}

 
