#include "QuantumDensity.h"
#include "SimulationEnvironment.h"
#include "Control.h"
#include "gnuplot_io.h"
using namespace std;


//===================================================================================//

QuantumDensity::QuantumDensity()
{
  quantum_model = NULL;

  

 
}

//============================================//
QuantumDensity:: ~QuantumDensity()
{
 
}

//============================================//
ID  QuantumDensity::convert_variable_name_to_id(const std::string& variable_name) const
{
  ID id = INVALID_ID;

  // for an empty string we return immediately
  //if (variable_name == "") return id;

  if (variable_name == "density")
    id = DENSITY;

  return id;

  
}

//============================================//
void QuantumDensity::get_solution_secure(const Elem* elem,
         const std::set<ID>& ids, std::vector<std::map<ID, double> >& values)
{
 
  

}
//============================================//
void QuantumDensity::get_solution_secure(const Elem* elem,
			 const std::vector<Point>& p, const std::set<ID>& ids,
			 std::vector<std::map<ID, double> >& values)
{

  if (ids.find(DENSITY) != ids.end())
  {
    const double coeff =  1.0/ ( (Constants::bohr_radius) * (Constants::bohr_radius) * (Constants::bohr_radius) * 1.0e6 );

    std::vector<double> density;

    get_particle_density( elem,  p,  density);

    unsigned int n = density.size();

    //values.resize(n);

    for (unsigned int i = 0; i < n; i++)
    {
      std::map<ID, double> point_map;

      point_map.insert(pair<ID, double>(DENSITY, density[i]*coeff));

      values[i] = point_map;
    }
   
  }
 
}

//============================================//




//============================================//
void QuantumDensity::get_particle_density(const Elem* element, const std::vector<Point>& quad_points, std::vector<double>& density)
{

  unsigned int n = quad_points.size();
  density.resize(n,0.0);

  map< const Elem*, double >::iterator it;
  it = real_space_density.find(element);
  if (it != real_space_density.end())
  {//if the element elem is active we put the same density for all the point
    unsigned int n = quad_points.size();
   
    for (unsigned int i = 0; i < n; i++)
      density[i] = it->second;
  }
  else
  {//if the element is not included, we have to check his children or parents
    //-------------------------------------------------------------------------
    //1) may be it has a parent the belongs to the  real_space_density map
    const Elem* el1 = element->parent();
    bool out = false;
    bool found = false;
    
    while ( !out )
    {
      if ( el1 != NULL )
      {
	it = real_space_density.find(el1);
	if ( it != real_space_density.end() )
	{
	  unsigned int n = quad_points.size();


	  for (unsigned int i = 0; i < n; i++)
	    density[i] = it->second;
	  
	  out = true;
	  found = true;
	}
	else
	{
	  el1 = el1->parent();
	}
      }
      else
      {
	out = true;
      }
    }
   
    if (!found) //2) may be it has a child that belongs to the real_space_density map
    {
      std::vector< const Elem * > active_children;
      element -> active_family_tree ( active_children, true);
      unsigned int n = active_children.size();

      for (unsigned int i1 = 0; ( i1 < n  ); i1++)
      {
	it  =  real_space_density.find(active_children[i1]);
	if (it !=  real_space_density.end() )
	{
	  //we have to check if this child contains a quadrature points
	  const unsigned int n = quad_points.size();
	  
	  for (unsigned int i = 0; i < n; i++)
	    if (active_children[i1]->contains_point(quad_points[i]) )
	    {
	      density[i] = it->second;
     
	    }
	  
	}
	
      

      }
    }
  }
}

//==================================================================================//
void QuantumDensity::build_elemental_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend)
{
 

  if (variables.find("quantum_density") != variables.end() )
  {
    legend.resize(1, "particle_density[cm^-3]");

    //------------------------------------------------------------------------------------//
    //!from map to vector
    Mesh & mesh = get_environment().get_mesh();
    
    map<const Elem*, double>::iterator it = real_space_density.begin();

    MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
    
    unsigned int    real_space_density_size = real_space_density.size();
    

    results.resize(real_space_density_size);

    unsigned int el_number = 0;

    for (; el !=end_el ; ++el)
    {
      const Elem* elem = *el; 

      it = real_space_density.find(elem);

      if (it != real_space_density.end())
	results[el_number] = it->second;
      else
	cerr << "WARNING: not all elements have a calculated density\n";

      el_number++;
    }
    //-------------------------------------------------------------------------------------//

    //now I have to transform from atomic units to [cm^-3]    

    const double coeff =  1.0/ ( (Constants::bohr_radius) * (Constants::bohr_radius) * (Constants::bohr_radius) * 1.0e6 );

    for (unsigned int i = 0; i < real_space_density_size; i++)
    {
           results[i] *= coeff;
    }
  


    //-------------------------------------------------------------------------------------//

  }



}



//=======================================================================================//
void QuantumDensity::do_init( )
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
      throw  InitFailedException("QuantumDensity: quantum_simulation " + quantum_simul_name + " does not exist");
  }
  else
  {
    throw  InitFailedException("QuantumDensity: quantum_simulation  has to be specified");
  }

 
  //--------------------------------------------------------------------------------------------//


}
 


//============================================//
void QuantumDensity::parse_options( )
{

  KspaceIntegration::parse_options();

  const ModelOptions& mod_opt = get_options();
 
  opt.Temperature             = mod_opt.get_option("Temperature", opt.Temperature);

  opt.log_output              = mod_opt.get_option("log_output", false);

  opt.intial_eigenstates_number = mod_opt.get_option("initial_eigenstates_number", 6);

 
  opt.analitic = mod_opt.get_option("analitic", true);

  opt.degeneracy  = mod_opt.get_option("degeneracy",1);

  if (opt.analitic)
  {
    if (k_dim == 1)
    {
      if (! mod_opt.find_option("k1") ) throw  InitFailedException("Kspace: k1 vectror must be defined"); 

     mod_opt.get_option("k1", k_vector1);

     if (k_vector1.size() != 3) throw  InitFailedException("Kspace: k1 vectror size must be equal to 3");

    }

    if (k_dim == 2)
    {
      if (! mod_opt.find_option("k1") ) throw  InitFailedException("Kspace: k1 vectror must be defined"); 

      mod_opt.get_option("k1", k_vector1);

      if (k_vector1.size() != 3) throw  InitFailedException("Kspace: k1 vector size must be equal to 3");


      
      
      if (! mod_opt.find_option("k2") ) throw  InitFailedException("Kspace: k2 vector must be defined"); 

      mod_opt.get_option("k2", k_vector2);
      
      if (k_vector2.size() != 3) throw  InitFailedException("Kspace: k2 vector size must be equal to 3");



    }

  }



 
}

//============================================//
void QuantumDensity:: do_solve()
{

  parse_options();

  real_space_density.clear();

  if (!opt.analitic)
  {//numerical integration
    KspaceIntegration::do_solve();
  }
  else
  {
    estimate_analitic_density();
  }

}





//================================================================//
void QuantumDensity::calculate_for_k_point(const Point& k_point, 
				     std::map<const Elem*, double>& density, 
				     double& integrated_quantity)
{

   
  vector<double> k_vector(3, 0.0);


  k_vector[0] = k_point(0);
  k_vector[1] = k_point(1);
  k_vector[2] = k_point(2);
	  

  ModelOptions quantum_model_opts;
	  

  quantum_model_opts.set_option("k_vector",  k_vector);



  quantum_model_opts.set_option("initial_eigenstates_number",opt.intial_eigenstates_number ); 
  quantum_model_opts["job"] = "density";

 
	  
  quantum_model->set_options(quantum_model_opts);


  quantum_model->solve();

  density = quantum_model->get_density();



  integrated_quantity = quantum_model->get_integrated_probability();

	 
 
}


//===============================================================================================// 

void QuantumDensity::estimate_analitic_density(void)
{




  ModelOptions quantum_model_opts;
  vector<double> k_vector(3, 0.0);

  quantum_model_opts.set_option("k_vector",  k_vector);



  quantum_model_opts.set_option("initial_eigenstates_number",opt.intial_eigenstates_number ); 
  quantum_model_opts["job"] = "density";
  quantum_model->set_options(quantum_model_opts);
 

  quantum_model->solve();

  vector<double> energy_k_0;

  quantum_model->get_eigenenergies (energy_k_0);

  unsigned int number_of_eigenstates = energy_k_0.size();

  

  vector<double> effective_mass(number_of_eigenstates);

  const Mesh& mesh = get_equation_systems().get_mesh();

  if (	k_dim == 2 )
  {
     

    
    //quantum_model_opts.set_option("initial_eigenstates_number",opt.intial_eigenstates_number ); 

    quantum_model_opts["job"] = "eigenstates";

    quantum_model_opts.set_option("number_of_eigenstates", number_of_eigenstates);
    
    double k_max1 = sqrt( k_vector1[0]*k_vector1[0] + k_vector1[1]*k_vector1[1] + k_vector1[2]*k_vector1[2]  );
    double k_max2 = sqrt( k_vector2[0]*k_vector2[0] + k_vector2[1]*k_vector2[1] + k_vector2[2]*k_vector2[2]  );
    
    vector<double> energy_k_1;
    vector<double> energy_k_2;
      
     
    quantum_model_opts.set_option("k_vector",  k_vector1);
    quantum_model->set_options(quantum_model_opts);
    
    quantum_model->solve();

    

    quantum_model->get_eigenenergies (energy_k_1);
    
    quantum_model_opts.set_option("k_vector2",  k_vector2);
    quantum_model->set_options(quantum_model_opts);
    
    quantum_model->solve();

    quantum_model->get_eigenenergies (energy_k_2);
     
    

    for (short i = 0; i < number_of_eigenstates; i++)
    {
       
      double imass1  = (2.0 * abs(energy_k_0[i] - energy_k_1[i] ) ) / Constants::Hartree /(k_max1 * k_max1);
      double imass2  = (2.0 * abs(energy_k_0[i] - energy_k_2[i] ) ) / Constants::Hartree /(k_max2 * k_max2);
      
      
 
      effective_mass[i] = 1.0/sqrt(imass1 * imass2);


       
      
    }

    quantum_model_opts.set_option("k_vector",  vector<double> (3, 0.0) );
    quantum_model->set_options(quantum_model_opts);
     
    quantum_model->solve();

    
   
    for (unsigned int i = 0; i < number_of_eigenstates; i++) 
    {
      map<const Elem*, double> state_density = quantum_model->estimate_density1D(i, effective_mass[i]);
       
       
      MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
      
      const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
      
      
      for (; el != end_el ; ++el)
      {
	const Elem* elem = *el;
	 
	real_space_density[elem] += state_density[elem] * opt.degeneracy;

      }
	 
    }
 
  }
  else if (k_dim == 1)
  {
     
    ModelOptions quantum_model_opts;
    
    vector<double> energy_k_1;
    
     quantum_model_opts.set_option("k_vector",  k_vector1);

     quantum_model_opts.set_option("initial_eigenstates_number",opt.intial_eigenstates_number ); 

     quantum_model_opts["job"] = "eigenstates";

     quantum_model_opts.set_option("number_of_eigenstates",number_of_eigenstates);

     quantum_model->set_options(quantum_model_opts);
    
     quantum_model->solve();

     quantum_model->get_eigenenergies (energy_k_1);
     
     double k_max = sqrt( k_vector1[0]*k_vector1[0] + k_vector1[1]*k_vector1[1] + k_vector1[2]*k_vector1[2]  );

     for (unsigned int i = 0; i < number_of_eigenstates; i++)
     {
       
       double imass  = (2.0 * abs(energy_k_0[i] - energy_k_1[i] ) ) / Constants::Hartree /(k_max * k_max);
       
       effective_mass[i] = 1.0/imass;

       quantum_model_opts.set_option("k_vector",  vector<double> (3, 0.0) );
       quantum_model->set_options(quantum_model_opts);
      
       quantum_model->solve();

    

       for (unsigned int i = 0; i < number_of_eigenstates; i++)
       {
	 map<const Elem*, double> state_density = quantum_model->estimate_density2D(i, effective_mass[i]);

       
	 MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
	 const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

       

	 for (; el != end_el ; ++el)
	 {
	   const Elem* elem = *el;
	    
	   real_space_density[elem] += state_density[elem] * opt.degeneracy;
	  
	    
	 }

       }

       

     }

     

   }


}


//==================================================================================================================//
