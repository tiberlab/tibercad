#include "QuantumDensity.h"
#include "SimulationEnvironment.h"
#include "Control.h"
#include "gnuplot_io.h"
using namespace std;



void QuantumDensity::get_particle_density(const Elem* element, const std::vector<double>& quad_points, std::vector<double> density)
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

//===================================================================================//

QuantumDensity::QuantumDensity()
{
  quantum_model = NULL;

  

 
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

 
}



//============================================//
QuantumDensity:: ~QuantumDensity()
{
 
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


//=================================================================// 
