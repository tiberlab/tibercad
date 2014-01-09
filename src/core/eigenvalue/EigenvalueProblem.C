// $Id$

#include "EigenvalueProblem.h"
#include "AtomisticStructure.h"
#include "Constants.h"
#include "Messages.h"
#include "DataOutput.h"

#include "elem.h"

#include <boost/shared_ptr.hpp>
#include <fstream>

using namespace std;

namespace 
{
  double scalar_prod(const vector<Complex>& a, const vector<Complex>& b)
  {
    Complex result = 0;
    for (size_t i = 0; i < a.size(); ++i)
      result += a[i] * conj(b[i]);

    return abs(result);
  }
}

void EigenvalueProblem::init_kspace(void)
{
   ModelOptions kopts = parse_kspace_options(ModelOptions());

   _kspace = new Kspace(kopts);

   if(_kspace==NULL)
     throw InitFailedException("Could not initialize k-space");
   else
     Messages::info("k-space initialized");

   /*
   do_dispersion = false;

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

     do_dispersion = true;
   }
   */

   /*
   if(get_options().has_submodel("DOS"))
   {
     ModelOptions::submodel_iterator it(get_options().submodels_begin("DOS"));

     const ModelOptions& opts = it->second;

     //disp_range[0] = opts.get_option("min_eigenvalue_number", 0);
     //disp_range[1] = opts.get_option("max_eigenvalue_number", 100);

     const ModelOptions kopts = parse_kspace_options(opts);

     _kspace = new Kspace(kopts);

     if(_kspace==NULL)
       throw InitFailedException("Could not initialize k-space");
     else
       Messages::info("k-space initialized");
   }
   */
}
 
ModelOptions EigenvalueProblem::parse_kspace_options(const ModelOptions& opts)
{
  ModelOptions kopts(opts);

  kopts.set_option("mesh_units", get_mesh_units());
  
  unsigned int k_dim = 3 - get_mesh().mesh_dimension();

  if (opts.find_option("k_space_dimension"))
    k_dim = opts.get_option("k_space_dimension", k_dim);

  kopts.set_option("k_space_dimension", k_dim);

  //std::vector<unsigned int>  num_nodes;

  if (opts.find_option("k_path") || opts.find_option("k-path"))
  {	  
    std::string kpath = opts.get_option("k-path","");
    kpath = opts.get_option("k_path",kpath);    
    kopts.set_option("k-path",kpath);
    //num_nodes.push_back(10);
    //kopts.set_option("number_of_nodes",num_nodes);
    //ModelOptions newopts;
    //newopts.set_option("output_format","grace");
    //set_options( newopts );
  }
 
  //if (opts.find_option("number_of_nodes"))
  //{
  //   opts.get_option("number_of_nodes",num_nodes);
  //   kopts.set_option("number_of_nodes", num_nodes);
  //}
  //else if (opts.find_option("number_of_elem"))
  //{
  //   opts.get_option("number_of_elem",num_nodes);
  //   for(int i=0; i< num_nodes.size(); i++)
  //         if(num_nodes[i]>0) ++num_nodes[i];

  //   kopts.set_option("number_of_nodes", num_nodes);
  //}

  //if (opts.find_option("wedge"))
  //  kopts.set_option("wedge", opts.get_option("wedge",""));

  // these are the real space lattice vectors, in nm
  RealVectorValue a(2, 0, 0), b(0, 2, 0), c(0, 0, 2);

  // if there is an atomistic structure, we can take the lattice vectors from it
  // (they come in Angstrom!)
  if (get_atomistic_structure() != NULL)
  {
    get_atomistic_structure()->get_lattice_vectors(a, b, c);
    a *= 0.1;
    b *= 0.1;
    c *= 0.1;
  }

  switch (k_dim)
  {
    case 1:
      kopts.set_option("r1", c);
      break;

    case 2:
      kopts.set_option("r1", b);
      kopts.set_option("r2", c);
      break;

    case 3:
      kopts.set_option("r1", a);
      kopts.set_option("r2", b);
      kopts.set_option("r3", c);
      break;

    default:
      break;
  }

  //kopts.set_option("k_space_basis", opts.get_option("k_space_basis", true));

/*
  double k_max = opts.get_option("k_max", 0.1);
  kopts.set_option("k_max", k_max);


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
*/

  kopts.set_option("mesh_order", opts.get_option("mesh_order", "first"));

  return kopts;
}





Point EigenvalueProblem::get_k_point(bool relative_coord) const
{
  Point kp(_k_vector[0], _k_vector[1], _k_vector[2]);
  if (relative_coord)
  {
    _kspace->inverse_transform(kp);
  }

  return kp;
}





void EigenvalueProblem::compute_dispersion(void)
{

  if(get_options().has_submodel("Dispersion"))
  {
    Messages m;
    m.info("Compute Dispersion ...");
    m.indent();

    ModelOptions::submodel_iterator it(get_options().submodels_begin("Dispersion"));

    const ModelOptions& opts = it->second;

    //disp_range[0] = opts.get_option("min_eigenvalue_number", 0);
    //disp_range[1] = opts.get_option("max_eigenvalue_number", 100);

    const ModelOptions kopts = parse_kspace_options(opts);

    // NOTE: this is not very elegant
    Kspace* original_kspace = _kspace;
    _kspace = new Kspace(kopts);

    if(_kspace==NULL)
      throw InitFailedException("Could not initialize k-space");
    else
      Messages::info("k-space for dispersion initialized");



    const Mesh* kmesh = _kspace->get_k_mesh();
    unsigned int number_of_k_points = kmesh->n_nodes();

    //if (disp_range[0] < 0) disp_range[0]=0;
    unsigned int number_of_eigs;

    std::cout<<"(EP) number of kpoints: " << number_of_k_points << std::endl;

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

    plot_dispersion();

    _kspace = original_kspace;
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

      if ((format == "grace") && (kdim > 1)) format = "vtk";

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


void
EigenvalueProblem::process_element(const Elem* elem, unsigned int entryside,
    vector<vector<eigen_problem_solution>>& ordered_solutions)
{
  ofstream of("test.dat", ofstream::app);

  bool already_done = true;

  // choose reference node
  unsigned int ref_node;
  for (unsigned int n = 0; n < elem->n_nodes(); ++n)
  {
    if (elem->is_node_on_side(n, entryside))
    { 
      ref_node = elem->node(n);
      if (ordered_solutions[ref_node].empty())
      {
        const Point& k_point = elem->point(n);

        solve_for_kpoint(k_point);
        int number_of_eigs = get_num_states();
        ordered_solutions[ref_node].resize(number_of_eigs);

        for (unsigned int j = 0 ; j < number_of_eigs; j++)
          ordered_solutions[ref_node][j] = _solution[j];
      }

      break;
    }
  }

  // the number of the reference solutions
  // (should be usually the same as the solution size)
  int ref_size = ordered_solutions[ref_node].size();

  for (unsigned int n = 0; n < elem->n_nodes(); ++n)
  {
    unsigned int node_id = elem->node(n);
    if (ordered_solutions[node_id].empty())
    {
      const Point& k_point = elem->point(n);

      solve_for_kpoint(k_point);
      int number_of_eigs = get_num_states();
      // we take the maximum only to make everything crash if the two numbers
      // do not correspond. However, this is usually sign of a badly posed
      // simulation setup
      ordered_solutions[node_id].resize(max(number_of_eigs, ref_size));


      //of << node_id << " ";
      //k_point.write_unformatted(of, true);

      set<unsigned int> ids;

      for (unsigned int k = 0 ; k < ref_size; k++)
      {
        unsigned int idx = k;
        double max_sp = 0;
        cerr << k << " : ";
        for (unsigned int j = 0 ; j < number_of_eigs; j++)
        {
          if (!ids.count(j))
          {
            double proj = scalar_prod(ordered_solutions[ref_node][k].eigen_vector,
                _solution[j].eigen_vector);
            cerr << j << " - " << proj << " ";
            if (proj > max_sp)
            {
              max_sp = proj;
              idx = j;
            }
          }
        }
        cerr << endl;
        ids.insert(idx);
        cerr << idx << " " << max_sp << endl;
        ordered_solutions[node_id][k] = _solution[idx];
      }

      already_done = false;
    }
  }

  // go into neighbours
  if (!already_done)
  {
    for (unsigned int n = 0; n < elem->n_sides(); ++n)
    {
      if (n != entryside)
      {
        const Elem* neigh = elem->neighbor(n);
        if (neigh != NULL)
        {
          int neigh_entry;
          for (unsigned int ns = 0; ns < neigh->n_sides(); ++ns)
          {
            if (neigh->neighbor(ns) == elem)
            {
              neigh_entry = ns;
              break;
            }
          }

          process_element(neigh, neigh_entry, ordered_solutions);
        }

      }
    }
  }
}


void EigenvalueProblem::_dos_for_kpoint(const Point& k_point,
    const Point& refpoint,
    DofField& density,
    double& integrated_quantity)
{
  ModelOptions& opts = get_options().submodels_begin("DOS")->second;
  double s = opts.get_option("gaussian_width", 0.01);

  unsigned int n_energy = _energy_mesh->n_nodes();

  // calculate for k-point, but only if it is not already there
  KSolutions::iterator it(_ksolutions.find(k_point));
  if (it == _ksolutions.end())
  {
    this->solve_for_kpoint(k_point);
    it = (_ksolutions.insert(make_pair(k_point, _solution))).first;
  }

  const vector<eigen_problem_solution>& solution = it->second;

  // now order, but only if k_point != refpoint
  if (k_point != refpoint)
  {
  }

  unsigned int n_eigs = solution.size();

  double a = 1.0 / (s * sqrt(2*M_PI));

  density.resize(n_energy);

  for (unsigned int n = 0; n < n_energy; n++)
  {
    double erg =  _energy_mesh->point(n)(0);

    double sum = 0.0;

    for (unsigned int i = 0; i < n_eigs; ++i)
    {
      double ediff = (erg - solution[i].eigen_energy) / s;
      double arg = 0.5 * ediff * ediff;
      sum += exp(-arg);
    }

    density[n] = a * sum;
  }
}




void EigenvalueProblem::calculate_dos(void)
{
  if (!get_options().has_submodel("DOS")) return;

  Messages::info("Compute DOS ...");

  ModelOptions& opts = get_options().submodels_begin("DOS")->second;

  delete _energy_mesh;
  _energy_mesh = new Mesh(1);

  double emin = opts.get_option("Emin", 0);
  double emax = opts.get_option("Emax", 5);
  unsigned int num_elem = static_cast<unsigned int>((emax - emin) / opts.get_option("dE", 0.001));

  MeshTools::Generation::build_cube (*_energy_mesh,
                                     num_elem, 0, 0,
                                     emin, emax,
                                     0, 0,
                                     0, 0,
                                     EDGE2);


  //
  // The simple approach integrates in k-space with a gaussian weight
  //
  ModelOptions kopts;
  if (opts.has_submodel("k-space"))
    kopts = opts.submodels_begin("k-space")->second;
  kopts += parse_kspace_options(kopts);

  _kspace = new Kspace(kopts);

  /*
//  KspaceIntegration* kint = KspaceIntegrationTemplate<EigenvalueProblem>::create(this,
  KspaceIntegration* kint = KspaceIntegration::create(this,
      &EigenvalueProblem::_dos_for_kpoint, kopts);
  kint->init();

  kint->solve();

  //DofField doff;
  //Point kp(0);
  //double dummy;
  //_dos_for_kpoint(kp, doff, dummy);

  vector<double> results(_energy_mesh->n_nodes(),0.0);

  vector<string> names(1, "DOS");

  string filename;
  string format = opts.get_option("output_format", "grace");

  DataOutput data_output(*_energy_mesh, format);
  data_output.set_output_directory(get_output_directory());

  filename = get_name() + "_dos" + TiberCad::get_filename_suffix();
  unsigned int n_erg = _energy_mesh->n_nodes();
  for(unsigned int n = 0; n < n_erg; n++)
  {
    results[n] = kint->get_solution()[n];
  }

  data_output.write_nodal_data(filename, results, names);
  */

  /*
  for (unsigned int i = 0; i < number_of_k_points; i++)
  {
    const Point&  k_point = kmesh->point(i);

    solve_for_kpoint(k_point);
    int number_of_eigs = get_num_states();
    solutions[i].resize(number_of_eigs);

    for (unsigned int j = 0 ; j < number_of_eigs; j++)
      solutions[i][j] = _solution[j];
//      solutions[i][j] = ptr_type(new eigen_problem_solution(_solution[j]));

  }
  */




  ///*

  const Mesh* kmesh = _kspace->get_k_mesh();
  unsigned int number_of_k_points = kmesh->n_nodes();

  //typedef boost::shared_ptr<eigen_problem_solution> ptr_type;

  //vector<vector<ptr_type>> solutions(number_of_k_points);
  vector<vector<eigen_problem_solution>> solutions(number_of_k_points);

  // order the solutions according to the bands
  process_element(kmesh->elem(0), 0, solutions);


  // this plots the dispersion
  {
    std::vector<double> results;
    std::vector<std::string> names;

    unsigned int number_of_eigs = solutions[0].size();
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
        results[number_of_eigs * j + i] = solutions[j][i].eigen_energy;
    }


    std::string filename(get_name() + "_dos");

    DataOutput data_output(*kmesh, "vtk");
    data_output.set_output_directory(get_output_directory());
    //data_output.set_filename(filename);

    data_output.write_nodal_data(filename, results, names);
  }
  //*/
}

ID
EigenvalueProblem::do_remember_current_solution(ID id)
{
  map<ID, vector<eigen_problem_solution>>::iterator end(_remembered_sol.end());
  map<ID, vector<eigen_problem_solution>>::iterator it(_remembered_sol.find(id));

  if (it != end)
    it->second = _solution;
  else
  {
    if (_remembered_sol.begin() == end)
      id = 1;
    else
      id = (--end)->first + 1;

    _remembered_sol[id] = _solution;
  }


  return id;
}


void
EigenvalueProblem::do_set_to_remembered_solution(ID id)
{
  map<ID, vector<eigen_problem_solution>>::iterator end(_remembered_sol.end());
  map<ID, vector<eigen_problem_solution>>::iterator it(_remembered_sol.find(id));

  if (it != end)
    _solution = it->second;
}


void
EigenvalueProblem::do_delete_remembered_solution(ID id)
{
  map<ID, vector<eigen_problem_solution>>::iterator end(_remembered_sol.end());
  map<ID, vector<eigen_problem_solution>>::iterator it(_remembered_sol.find(id));

  if (it != end)
    _remembered_sol.erase(it);
}


void EigenvalueProblem::do_plot(void)
{

  SimulationInterface::do_plot();

  compute_dispersion();

  calculate_dos();
}
 
void EigenvalueProblem::solve_for_kpoint(const Point& kpoint)
{
  const Point oldk(_k_vector[0], _k_vector[1], _k_vector[2]);
  set_k_point(kpoint);
  do_solve_for_kpoint(kpoint);
  set_k_point(oldk);
  k_is_old();
}


void
EigenvalueProblem::do_solve_for_kpoint(const Point& k_point)
{
  reinit();
  solve();
}


void EigenvalueProblem::get_eigenvalues(const std::string& particle, 
					std::vector<double>& values) const
{

  unsigned int n = _solution.size();
  unsigned int num_st = 0;
  values.reserve(n);

  for (unsigned int i = 0; i < n; i++)
  {
    if(particle.empty() || (_solution[i].particle == particle))
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

  return 0;
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

 
