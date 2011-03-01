// $Id$

#include "SimulationEnvironment.h"
#include "EnvelopFunctionApprox.h"
#include "ModelOptions.h"
#include "EFAbulkModel.h"
#include "Material.h"
#include "Boundary.h"
#include "TiberMath.h"
#include <gnuplot_io.h>
#include "SimulationOptions.h"

#include "EigenSolver.h"


#include <edge_edge2.h>


#include "Messages.h"



extern "C"
{
   void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& LDA,
	       double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info);
};

using namespace std;
using namespace Constants;

Device*   EnvelopFunctionApprox:: _device;


//---------------------------------------------------------------------------------//


/*
inline void EnvelopFunctionApprox::get_electric_potential(const Elem* elem, const std::vector<Point>& q_point,
						   std::vector<double> electric_potential) const
{
  poisson_equation->get_solution (elem, q_point, potential_ID, electric_potential);
}
*/
//---------------------------------------------------------------------------------//

inline double EnvelopFunctionApprox::get_band_edge( const Elem* elem) const
{
  vector<double> values;
  vector<Point> qp(1, elem->centroid());

  if (!poisson_equation->get_solution(elem, band_edge_ID, values, qp))
    throw ModelErrorException("efaschroedinger needs solved Poisson equation");

  return values[0];
}

//---------------------------------------------------------------------------------//
inline double EnvelopFunctionApprox::get_electric_potential(const Elem* elem, const Point&  qpoint) const
{
  vector<double> values;
  vector<Point> qp(1, qpoint);

  poisson_equation->get_solution(elem, potential_ID, values, qp);

  return values[0];
}
//---------------------------------------------------------------------------------//

inline double EnvelopFunctionApprox::get_electro_chem_potential(const Elem* elem) const
{
  vector<double> values;
  vector<Point> qp(1, elem->centroid());

  poisson_equation->get_solution (elem, electro_chem_pot_ID, values, qp);

  return values[0];
}




//---------------------------------------------------------------------------------//
void EnvelopFunctionApprox::build_integrated_quantities (std::vector< double > &values)
{



  values.resize(0);




  if (plot_solution("EigenEnergy"))
  {
    unsigned int n = solution.size();
    values.resize(n);
    for (unsigned int i = 0; i < n; i++)
    {
      values[i] = solution[i].eigen_energy;
    }
  }


  if (plot_solution("Occupation"))
  {
    unsigned int n = solution.size();
    unsigned int m = values.size();
    values.resize(n + m);
    for (unsigned int i = 0; i < n; i++)
    {
      values[i + m] = Fermi_statistics_probability(solution[i].eigen_energy, solution[i].Fermi_energy,solution[i].Temperature);
    }
  }

}

//---------------------------------------------------------------------------------//

void EnvelopFunctionApprox::build_integrated_quantities_description (
    std::vector< std::string > &legend,
    std::vector< std::string > &description)
{

  legend.resize(0);

  if (plot_solution("EigenEnergy"))
  {
    unsigned int n = solution.size();
    legend.resize(n);
    for (unsigned int i = 0; i < n; i++)
    {


      ostringstream temp;
      temp << i ;
      legend[i] = temp.str();
    }
  }

  description.resize(1);
  description[0] = "Eigen energy [eV]";

  if (plot_solution("Occupation"))
  {
    unsigned int n = solution.size();
    unsigned int m = legend.size();
    legend.resize(m + n);
    for (unsigned int i = 0; i < n; i++)
    {
      ostringstream temp;
      temp << i ;
      legend[i+m] =  temp.str() ;
    }

    description.push_back("Occupation Probability");

  }


}




//---------------------------------------------------------------------------------//

void EnvelopFunctionApprox::build_elemental_results (const std::set< std::string > &variables,
					 std::vector< double > &results, std::vector< std::string > &legend)
{
 /*
  double parallel_mass = 0.067* 0.95 +  0.026 * 0.05;
  const set<string>::const_iterator varend(variables.end());
  if (variables.find("test_density") != varend)
  {

    results = estimate_density( parallel_mass);

    legend.resize(1);
    legend[0] = string("density");
  }
*/

}

//---------------------------------------------------------------------------------//
void EnvelopFunctionApprox::build_nodal_results(const std::set<std::string>& variables,
						std::vector<double>& results, std::vector<std::string>& legend)
{

  if (! (opt.job == BULKEIGENSTATES || opt.job == BULKDENSITY ) )
  {

    const set<string>::const_iterator varend(variables.end());

    const MeshBase& mesh1 = system->get_mesh();


    const bool output_psi = variables.find("EigenFunctions") != varend ;
    const bool output_level = variables.find("EnergyLevels") != varend ;

    MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
    const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();

    unsigned int number_of_points = 0;
    unsigned int num_var = 0;
    const unsigned int num_functions = solution.size();
    unsigned int temp = 0;

    for ( ; nd != nd_el ; ++nd)  number_of_points++;

    if (output_psi)  num_var += num_functions;

    if (output_level)  num_var += num_functions;

    if (output_psi || output_level)
    {

      legend.resize(num_var);
      results.resize(number_of_points * num_var);


      if (output_psi)
      {
	temp = num_functions;
	for (unsigned int i = 0; i < num_functions; i++)
	{
	  std::ostringstream i_str;
	  i_str << "state_number_" << i ; //The states are numbered starting from 0
	  legend[i] = i_str.str();

	  std::vector<double> prob_data(number_of_points, 0.0);



	  prepare_probability_function(i,  prob_data);

	  for (unsigned int i1 = 0; i1<number_of_points; i1++)
	    results[num_var * i1 + i] = prob_data[i1];


	}
      }

      if (output_level)
      {
	for (unsigned int i = 0; i < num_functions; i++)
	{
	  std::ostringstream i_str;
	  i_str << "energy_number_" << i ; //The states are numbered starting from 0
	  legend[i + temp] = i_str.str();


	  for (unsigned int i1 = 0; i1<number_of_points; i1++)
	    results[num_var * i1 + i + temp] = solution[i].eigen_energy;

	}
      }



    }

  }


}

//---------------------------------------------------------------------------------//



void EnvelopFunctionApprox::get_integrated_quantities(const std::set<std::string>& names,
						      std::vector<double>& values)
{


}

//---------------------------------------------------------------------------------//

void EnvelopFunctionApprox::get_eigenenergies(std::vector<double>& values) const
{


  unsigned int n = solution.size();
  values.resize(n);
  for (unsigned int i = 0; i < n; i++)
  {
    values[i] = solution[i].eigen_energy;
  }


}

//---------------------------------------------------------------------------------//

void EnvelopFunctionApprox::get_occupations(std::vector<double>& values) const
{


  unsigned int n = solution.size();
  values.resize(n);

  for (unsigned int i = 0; i < n; i++)
  {
    values[i ] = Fermi_statistics_probability(solution[i].eigen_energy, solution[i].Fermi_energy, solution[i].Temperature);
  }

}




//====================================================================================//

const std::vector< EnvelopFunctionApprox::eigen_propblem_solution >& EnvelopFunctionApprox::get_solution() const
{
  return(solution);
}


//====================================================//
PhysicalModel* EnvelopFunctionApprox::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{


  EFAbulkModel* model = dynamic_cast<EFAbulkModel*> ( PhysicalModelInterface::create("EFAmodel", options) );

  if (model == NULL)
    throw ModelErrorException("EnvelopFunctionApprox: cannot create EFAbulkModel");

  return(model);

}
//=====================================================//
BoundaryProperties* EnvelopFunctionApprox::create_boundary_model(const ModelOptions& options) const  throw (ModelErrorException)
{

  return NULL;
}


//====================================================//
double EnvelopFunctionApprox::get_band_edge() const
{



  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();




  vector <double> electric_potential(1, 0.0);


  double band_edge;




  for (; el != end_el ; ++el )
  {
    const Elem* elem = *el;





    double temp = get_band_edge( elem);



    if ( elem == *(mesh->active_elements_begin()) )
    {
      band_edge = temp;
    }





    if (opt.particle == "el")
    {

      if (band_edge > temp)
	band_edge = temp;
    }
    else
    {


      if (band_edge < temp)
	band_edge = temp;

    }





  }





  return(band_edge);




}

//===================================================//
EnvelopFunctionApprox::EnvelopFunctionApprox(const ModelOptions& options)
 : FEMEigenvalueProblem(options)
{
  poisson_equation = NULL;

  temperature_simulation = NULL;

  strain = NULL;

  _bulk_mat_element = NULL;
}
//===================================================//
void EnvelopFunctionApprox::parse_options()
{

  FEMEigenvalueProblem::parse_options();

  const ModelOptions& mod_opt = get_options();



  opt.particle                = mod_opt.get_option("particle", "el");


  opt.log_output              = mod_opt.get_option("log_output", false);


  opt.output_type             = mod_opt.get_option("output_type","GMV");

  opt.eigen_number_increase_factor = mod_opt.get_option("eigen_number_increase_factor",1.2);



  opt.relative_density_tolerance =  mod_opt.get_option("relative_density_tolerance", 1e-2);

  //-------------------------------------------------------------------------------------------//
  //strain model
  std::string  strain_model_name = mod_opt.get_option("strain_model_name","no_strain");

  if (strain_model_name != "no_strain" )
  {
    opt.consider_strain = true;

    strain = dynamic_cast< Macrostrain* > ( find_simulation( strain_model_name )   );

    if (strain == NULL)
      throw InitFailedException("Unknown strain model" + strain_model_name);


  }
  else
  {
    opt.consider_strain = false;
  }


  //-------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------//
  //poisson
  std::string  poisson_model_name = mod_opt.get_option("poisson_model_name","no_poisson");
  if ( poisson_model_name != "no_poisson" )
  {
    opt.consider_potential = true;

    poisson_equation  = find_simulation ( poisson_model_name );

    if (poisson_equation == NULL)
      throw InitFailedException( "Unknown poisson model " + poisson_model_name);

    potential_ID = poisson_equation->get_solution_id("ElPotential");



    if (potential_ID ==  INVALID_ID)
      throw InitFailedException( "Unknown variable ");


    if (opt.particle == "el")
    {
      electro_chem_pot_ID = poisson_equation->get_solution_id("eQFermi");
      band_edge_ID = poisson_equation->get_solution_id("Ec");

    }
    else if (opt.particle=="hl")
    {
      band_edge_ID = poisson_equation->get_solution_id("Ev");
      electro_chem_pot_ID = poisson_equation->get_solution_id("hQFermi");

    }




  }
  else
  {
    opt.consider_potential = false;
  }





  //--------------------------------------------------------------------------------------------//
  //as default, we  estimate spectrum shift only in electric potential is defined
  opt.estimate_spectrum_shift =  opt.consider_potential;

  opt.estimate_spectrum_shift = mod_opt.get_option("estimate_spectrum_shift",  opt.estimate_spectrum_shift);

  if ( !opt.consider_potential && opt.estimate_spectrum_shift) throw InitFailedException( "EnvelopeFunctionApprox: cannot estimate spectrum shift without electric potential");

  //--------------------------------------------------------------------------------------------//

  if (opt.estimate_spectrum_shift)
  {
    opt.spectrum_shift = get_band_edge();

  }
    //--------------------------------------------------------------------------------------------//


  //k-vector
  std::vector<double> k_vec(3, 0.0);
  mod_opt.get_option("k_vector",k_vec);
  if (k_vec.size() == 3)
  {

    for (short i = 0; i < 3; i++) k_vector[i] = k_vec[i];


  }
  else
    throw InitFailedException( "EnvelopeFunctionApprox: k_vector size must be equal to 3 instead of " + k_vec.size());

  //--------------------------------------------------------------------------------//
  std::string  job_name = mod_opt.get_option("job","eigenstates");
  if (job_name == "eigenstates")
    opt.job = EIGENSTATES;
  else if (job_name == "density")
    opt.job = DENSITY;
  else if (job_name == "bulk")
    opt.job = BULKEIGENSTATES;
  else if (job_name == "bulkdensity")
    opt.job = BULKDENSITY;
  else
    throw InitFailedException( "EnvelopeFunctionApprox: Incorrect job " + job_name );
  //---------------------------------------------------------------------------------//
  std::string  heat_model_name = mod_opt.get_option("heat_model","");

  if ( heat_model_name != "" )
  {
    temperature_simulation  = find_simulation ( heat_model_name );

    if (temperature_simulation == NULL)
      throw InitFailedException( "Unknown heat model " + heat_model_name);

    temperature_ID = temperature_simulation->get_solution_id("temperature");

  }


  temp_interface.set_simulation(heat_model_name);

  //default value for temperature
  opt.Temperature = mod_opt.get_option("Temperature", SimulationOptions::temperature);




  //---------------------------------------------------------------------------------//

  opt.first_state = mod_opt.get_option("first_state", 0);
  opt.initial_eigenstates_number = mod_opt.get_option("initial_eigenstates_number", 6);

  //---------------------------------------------------------------------------------//

  std::set<const Node*> used_nodes;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();
  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const unsigned int n = elem->n_nodes();
    for (unsigned int i = 0; i < n; i++)
    {
      const Node* nd = elem->get_node(i);
      used_nodes.insert(nd);
    }
  }

  number_of_nodes = used_nodes.size();

  opt.convergent_density = mod_opt.get_option("convergent_density", true);


  opt.local_occupation = mod_opt.get_option("local_occupation", true);



  //--------------------------------------------------
  //point for bulk calculation
  if (opt.job == BULKEIGENSTATES || opt.job == BULKDENSITY)
  {
    if (mod_opt.find_option("bulk_point"))
    {
      vector<double> point;

      mod_opt.get_option("bulk_point",point);
      for (short i = 0; i < 3; i++)  _bulk_point(i) = point[i];
    }
    else
    {
      throw InitFailedException( "You have to specify a bulk_point");
    }
  }


}




//===================================================//


//===================================================//
void EnvelopFunctionApprox::do_init( )
{

  FEMEigenvalueProblem::do_init();


  SimulationEnvironment& si = get_environment();

  _device = &( si.get_device() );




  es = &(get_equation_systems());

  mesh = &(es->get_mesh());

  system_name = get_equation_system_name ( );

  es->add_system<LinearImplicitSystem> (system_name);

  system = &( es->get_system<LinearImplicitSystem>( system_name ) );

  dim = mesh->mesh_dimension();

  //--------------------------------------------------------------------------------------------------------//
  opt.number_of_bands = calculate_number_of_bands( );


  //--------------------------------------------------------------------------------------------------------//
  //add variables
  psi_name.clear();
  for (short i = 0; i < opt.number_of_bands; i++)
    {
      std::ostringstream var_str;
      var_str << "psi" << i ;
      string name = var_str.str();
      psi_name.push_back(name);

      system->add_variable(name,FIRST);
    }

  //add matrixes
  //---------------------------------------------------------------------------------------------------------//

  DofMap& dof_map = system->get_dof_map();

  system->add_matrix("Ham_real"); //add matrix for a real part of the Hamiltonian

  Ham_real = & (system->get_matrix("Ham_real"));

  system->add_matrix("Ham_imag");//add matrix for an imaginary part of the Hamiltonian

  Ham_imag = &(  system->get_matrix("Ham_imag") );

  system->add_matrix("S_real"); //add matrix for S matrix

  S_real = &( system->get_matrix("S_real") );



  //---------------------------------------------------------------------------------------------------------//

  Scaling& scaling = get_scaling();

  scaling.set_length_scaling(Constants::bohr_radius);

  scaling.set_calc_mesh_units(_device->get_mesh_units());



   //--------------------------------------------------------------------------------------------------------//
   MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
   const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

   bool temp = true;

   for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;
      short n1 = elem->n_nodes();
      for (short i1 = 0; i1 < n1 ; i1++)
	{


	  const Point& p = elem->point(i1);
	  for (unsigned i = 0; i < 3; i++)
	  {

	    if (temp)
	    {
	      min_coord[i] = p(i);
	      max_coord[i] = p(i);
	      temp = false;
	    }
	    else
	    {

	      if (min_coord[i] < p(i)) min_coord[i] = p(i);
	      if (max_coord[i] > p(i)) max_coord[i] = p(i);
	    }
	  }

	}

    }

   //-------------------------------------------------------------------------------------------------------//

   system->init();



   //------------------------------------------------------------------------------------------------------//
   //peiodicity can not be changed between runs because that will require cleaning of the DOF constraint table

   /*
     const ModelOptions& mod_opt = get_options();

     opt.periodicity[0]          = mod_opt.get_option("x-periodicity", false);
     opt.periodicity[1]          = mod_opt.get_option("y-periodicity", false);
     opt.periodicity[2]          = mod_opt.get_option("z-periodicity", false);
   */

   //------------------------------------------------------------------------------------------------------//
   //kp bands map
   {
     MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
     const Elem* elem = *el;
      const ID subdomain = elem->subdomain_id();
      const Material* mat = _device->get_material(subdomain);

      EFAbulkHamiltonian* element_hamiltonian =
	(  dynamic_cast<EFAbulkModel*> (  mat ->get_model(get_id()) )  )->get_Hamiltonian_model();

      opt.kp_bands = element_hamiltonian->get_kp_bands_map();

   }
   //------------------------------------------------------------------------------------------------------//


}
//===========================================================//
void EnvelopFunctionApprox::do_solve()
{

 parse_options();

 if (opt.job == BULKEIGENSTATES )
 {


   solve_bulk();



 }
 else
 {

   if  ( opt.job ==   EIGENSTATES || opt.job == DENSITY )
   {
     if (solver_opt.Dirichlet_bc_everywhere)
       apply_diriclet_bc_at_all_boundaries();
     else
       create_dirichlet_dofs();


     make_constraints(); //creates a copy of them


     make_nodes_periodic();

     apply_periodic_bc();

     make_new_dofs();
   }


   if ( opt.job ==   EIGENSTATES )
     solve_eigen_value_problem( solver_opt.number_of_eigenstates);
   else if ( opt.job == DENSITY || opt.job == BULKDENSITY)
     calculate_convergent_density();

 }

}
//===========================================================//
void EnvelopFunctionApprox::calculate_Hamiltonian_and_S(void)
{


  Ham_real->zero();
  Ham_imag->zero();
  S_real->zero();

  //material list
  //assemble_material_list();

 vector<unsigned int> psivar(opt.number_of_bands);
 //get numbers of variables
 for (unsigned int i = 0; i < opt.number_of_bands; i++)
 {
   psivar[i] = system->variable_number(psi_name[i]);
 }

 DofMap& dof_map = system->get_dof_map();

 FEType fe_type = dof_map.variable_type(psivar[0]); //all the variable have the same FE representation



 AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.

  fe -> attach_quadrature_rule (&qrule);

 // Here we define some references to cell-specific data that
 // will be used to assemble the linear system.
 //
 // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  //------------------------------------------------------------
 std::vector<unsigned int> dof_indices_component;

 std::vector<unsigned int> dof_indices;

  //-------------------------------------------------------------
  //matrixes to built the system

  DenseMatrix<Number> ham_real;
  DenseMatrix<Number> ham_imag;
  DenseMatrix<Number> s_real;


  DenseSubMatrix<Number> ham_real_sub(ham_real);
  DenseSubMatrix<Number> ham_imag_sub(ham_imag);
  DenseSubMatrix<Number> s_real_sub(s_real);


  vector<double> box_volume(Ham_real->n(), 0.0);

  if (solver_opt.discretization_method == BIM)
  {
    MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh->active_elements_end();
    for ( ; el != end_el ; ++el)
    {//el
      const Elem* elem = *el;
      dof_map.dof_indices (elem, dof_indices);
      for (unsigned int band1 = 0; band1 < opt.number_of_bands; band1++)
      {//band1
	dof_map.dof_indices (elem, dof_indices_component, psivar[band1]);
	const unsigned int n_psi_dofs = dof_indices_component.size();
	for (unsigned int p1=0; p1<n_psi_dofs; p1++)
	{
	  double box_part_volume;

	  Elem* box_part_elem;

	  if (dim == 1)
	  {

	    Edge2 box_part_1D;


	    Node point1(elem->point(0),999);
	    Node p_center(elem->centroid(),1000);

	    if (p1 == 0)
	    {
	      box_part_1D.set_node(0) = &point1;
	      box_part_1D.set_node(1) = &p_center;

	    }
	    else
	    {
	      box_part_1D.set_node(1) = &point1;
	      box_part_1D.set_node(0) = &p_center;
	    }

	    box_part_elem = &box_part_1D;


	    {//volume calculation

	      FEType fe_type (FIRST , LAGRANGE);

	      //AutoPtr<FEBase> fe (FEBase::build(dim,  fe_type));

	      AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

	      QGauss qrule (dim, FIRST);

	      fe -> attach_quadrature_rule (&qrule);



	      const std::vector<Real>& JxW = fe->get_JxW();

	      fe->reinit(box_part_elem);


	      box_part_volume = 0.0;
	      for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
		box_part_volume += JxW[qp];

	    }

	  }






	  // double v = box_part_1D.volume();

	  //  box_part_volume = 0.5 * std::abs(p1(0) - p2(0)) ;



	  box_volume[dof_indices_component[p1]] += box_part_volume;



	}
      }

    }

  }



  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  Tensor2Sym strain_crystal_system(0);
  double electric_potential = 0;

  EFAbulkHamiltonian* element_hamiltonian;

  unsigned int el_number = 0;







  for ( ; el != end_el ; ++el)
    {//el
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      const ID subdomain = elem->subdomain_id();
      const Material* mat = _device->get_material(subdomain);



      element_hamiltonian = (  dynamic_cast<EFAbulkModel*> (  mat ->get_model(get_id()) )  )->get_Hamiltonian_model();

      element_hamiltonian->set_temperature(temp_interface.get_temperature( elem, elem->centroid()));

      element_hamiltonian->set_k_vector(k_vector);

      element_hamiltonian->calculate_Hamiltonian_k_par();


      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs   = dof_indices.size();

      ham_real.resize(n_dofs, n_dofs);
      ham_imag.resize(n_dofs, n_dofs);
      s_real.resize(n_dofs, n_dofs);

      // complex<double> operator_sign = Complex(0.0, -1.0);
      if (solver_opt.discretization_method == FEM)
      {//FEM
	fe->reinit (elem);

	for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{//qp
	  //--------------------------------------------------------------------------------
	  /*
	    We assume that strain and electric potential may be different for different quadrature points
	    It is done for a sake of a multiscale generalization
	  */
	  if (opt.consider_strain)
	  {

	    strain_crystal_system = strain->get_strain_crystal(elem , q_point[qp]);
	  }


	  if (opt.consider_potential)
	  {
	    electric_potential = get_electric_potential( elem, q_point[qp] );
	  }




	  element_hamiltonian->apply_strain_and_potential(strain_crystal_system, electric_potential);

	  //------------------------------------------------------------------------------------------


	  std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
	    model_Ham = ( element_hamiltonian->get_Hamiltonian() );




	  for (unsigned int band1 = 0; band1 < opt.number_of_bands; band1++)
	  {//band1
	    dof_map.dof_indices (elem, dof_indices_component, psivar[band1]);
	    const unsigned int n_psi_dofs = dof_indices_component.size();

	    for (unsigned int band2 = 0; band2 < opt.number_of_bands; band2++)
	    {//band2

	      //Hamiltonian


	      ham_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
	      ham_imag_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
	      for (unsigned int p1=0; p1<n_psi_dofs; p1++)
	      {
		for (unsigned int p2=0; p2<n_psi_dofs; p2++)
		{
		  complex<double> value = (0.0, 0.0);
		  //constant
		  value += JxW[qp] * phi[p1][qp] * phi[p2][qp] * model_Ham[band1][band2].constant ;


		  //linear left

		  for (short i = 0; i < dim; i++)
		  {
		    value -= JxW[qp]* dphi[p1][qp](i) * phi[p2][qp] * model_Ham[band1][band2].linear_left[i]
		      * Complex(0.0, -1.0);
		  }
		  //linear right

		  for (short i = 0; i < dim; i++)
		  {
		    value += JxW[qp]* dphi[p2][qp](i) * phi[p1][qp] * model_Ham[band1][band2].linear_right[i]
		      * Complex(0.0, -1.0);

		  }

		  //quadratic

		  for (short i = 0; i < dim; i++)
		    for (short j = 0; j < dim; j++)
		    {
		      value -= JxW[qp] * dphi[p1][qp](i) * dphi[p2][qp](j)*model_Ham[band1][band2].quad[i][j]
			* Complex(0.0,-1.0) * Complex(0.0, -1.0);

		    }






		  ham_real_sub(p1,p2) += value.real();
		  ham_imag_sub(p1,p2) += value.imag();

		}
	      }




	      //S-matrix
	      if (band1 == band2)
	      {
		s_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
		for (unsigned int p1=0; p1<n_psi_dofs; p1++)
		{
		  for (unsigned int p2=0; p2<n_psi_dofs; p2++)
		  {
		    s_real_sub(p1,p2) += JxW[qp] * phi[p1][qp] * phi[p2][qp];
		  }
		}

	      }
	      //--------------------------------------------------------------------------//

	    }
	  }



	}

      }

      if (solver_opt.discretization_method == BIM)
      {//BIM

	Point center = elem->centroid();

	if (opt.consider_strain)
	{

	  strain_crystal_system = strain->get_strain_crystal(elem, center);
	}


	if (opt.consider_potential)
	{

	  electric_potential = get_electric_potential( elem, center );
	}


	element_hamiltonian->apply_strain_and_potential(strain_crystal_system, electric_potential);




	std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
	  model_Ham = ( element_hamiltonian->get_Hamiltonian() );


 	double box_part_volume;
	{
	  Point p1 = elem->point(0);
	  Point p2 = elem->point(1);

	  box_part_volume = 0.5 * std::abs(p1(0) - p2(0)) ;
	}

	for (unsigned int band1 = 0; band1 < opt.number_of_bands; band1++)
	{//band1
	  dof_map.dof_indices (elem, dof_indices_component, psivar[band1]);
	  const unsigned int n_psi_dofs = dof_indices_component.size();

	  for (unsigned int band2 = 0; band2 < opt.number_of_bands; band2++)
	  {//band2

	    ham_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
	    ham_imag_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);

	    //hamiltonian

	    for (unsigned int p1=0; p1<n_psi_dofs; p1++)
	    {//p1

	      Elem* box_part_elem;
	      vector< Point > qpoints;

	      Point q1(0.0, 0.0, 0.0);
	      qpoints.push_back(q1);
	      qpoints.push_back( FEInterface::inverse_map ( dim, fe_type, elem, elem->point(p1)) );

	      fe->reinit (elem, &qpoints);

	      for (unsigned int p2=0; p2<n_psi_dofs; p2++)
	      {	//p2
		complex<double> value = (0.0, 0.0);
		//constant (zero order)

		if (p1 == p2)
		{
		  value +=  model_Ham[band1][band2].constant * box_part_volume;
		  if (band1 == band2) value += - opt.spectrum_shift/Hartree * box_part_volume;
		}


	       //linear left

		for (short i = 0; i < dim; i++)
		{
		  double surface ;
		  if (p1 == 0)
		    surface = 1.0;
		  else
		    surface = -1.0;

		  for (unsigned int qp=0; qp < dim; qp++)
		    value +=   phi[p2][qp] * model_Ham[band1][band2].linear_left[i] * surface
		      * Complex(0.0, -1.0);

		}

		//linear right

		for (short i = 0; i < dim; i++)
		  value +=  dphi[p2][dim](i) * model_Ham[band1][band2].linear_right[i]
		    * Complex(0.0, -1.0)  * box_part_volume ;


		//second order
		vector<double> n1[3];

		for (unsigned int qp=0; qp < dim; qp++)
		{
		  for (short i = 0; i < dim; i++)  //x,y,z
		    for (short j = 0; j < dim; j++) //x,y,z
		    {


		      double surface ;

		      if (p1 == 0)
			surface = 1.0;
		      else
		     	surface = -1.0;


		      value +=   dphi[p2][qp](j)  * surface *
			model_Ham[band1][band2].quad[i][j] * Complex(0.0,-1.0) * Complex(0.0, -1.0);

		    }
		}





		ham_real_sub(p1,p2) += value.real()/box_volume[dof_indices_component[p1]];
		ham_imag_sub(p1,p2) += value.imag()/box_volume[dof_indices_component[p1]];


	      }
	    }




	  }
	}

      }




      if (solver_opt.discretization_method == FEM) ham_real.add( - opt.spectrum_shift/Hartree, s_real);//apply spectrum shift.

      vector<unsigned int> dof_indices_tmp;

      if (solver_opt.discretization_method == FEM)
      {
	dof_indices_tmp = dof_indices;

	dof_map.constrain_element_matrix(s_real, dof_indices_tmp);
	S_real->add_matrix(s_real,dof_indices_tmp);
      }

      dof_indices_tmp = dof_indices;

      dof_map.constrain_element_matrix(ham_real, dof_indices_tmp);
      Ham_real->add_matrix(ham_real,dof_indices_tmp);

      dof_indices_tmp = dof_indices;

      dof_map.constrain_element_matrix(ham_imag, dof_indices_tmp);
      Ham_imag->add_matrix(ham_imag,dof_indices_tmp);



      el_number++;

    }


//this is only to test
/*
  Ham_real->print_matlab("ham_r_matlab.m");
  Ham_imag->print_matlab("ham_i_matlab.m");
  S_real->print_matlab("s.m");
*/




  copy_H_matrix_to_solver( );



  if (solver_opt.discretization_method == FEM)  copy_S_matrix_to_solver( );





  //  dof_map.print_dof_constraints();

}
//============================================================//
void EnvelopFunctionApprox::copy_S_matrix_to_solver()
{

  // 2010-11-08 It seems this is just wasting memory
  //EigenSolver::init_S_matrix(number_of_new_dofs);

  int size_matrix = S_real->m();

  int non_zeros_number[number_of_new_dofs];

  PetscMatrix<Number>* p_matrix = static_cast<PetscMatrix<Number>* >(S_real);

  p_matrix->close();


  //----------preallocate memory------------------------------------------------------

  for (int row = 0 ; row < size_matrix; row++)
  {
    if (new_dofs[row].independent)
    {
      int ierr = 0;
      const  PetscScalar *petsc_row_vals;
      const PetscInt *petsc_cols;
      int n_cols = 0;

      ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      vector<unsigned int> column_vector;
      vector<Complex> row_values;

      non_zeros_number[new_dofs[row].new_number] = 0;

      for (int col = 0; col < n_cols; col++)
      {
	if ((new_dofs[petsc_cols[col]].independent) &&  (petsc_row_vals[col] != 0.0))
	{

	  non_zeros_number[new_dofs[row].new_number]++;


	}
      }

      ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

    }
  }


  EigenSolver::preallocate_S_matrix(number_of_new_dofs,  non_zeros_number);

  //--------------assebmle data--------------------------------------------------------


  for (int row = 0 ; row < size_matrix; row++)
  {
    if (new_dofs[row].independent)
    {
      int ierr = 0;
      const  PetscScalar *petsc_row_vals;
      const PetscInt *petsc_cols;
      int n_cols = 0;

      ierr = MatGetRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      vector<unsigned int> column_vector;
      vector<Complex> row_values;

      for (int col = 0; col < n_cols; col++)
      {
	if ((new_dofs[petsc_cols[col]].independent) &&  (petsc_row_vals[col] != 0.0))
	{



	  double value = petsc_row_vals[col];
	  double zero = 0.0;

	  column_vector.push_back(new_dofs[petsc_cols[col]].new_number);
	  row_values.push_back(Complex(value, zero));




	}
      }

      EigenSolver::insert_S_row( new_dofs[row].new_number, column_vector, row_values);

      ierr = MatRestoreRow(p_matrix->mat(), row ,&n_cols, &petsc_cols,&petsc_row_vals);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

    }
  }

  EigenSolver::finalize_S_assembly();

}



//============================================================//



//=============================================================//
double EnvelopFunctionApprox::get_new_spectrum_shift(void)
{
  double st_shift_value ;


  read_SLEPC_solution(1);

  assert(solution.size() == 1);


  st_shift_value = (solution[0].eigen_energy - opt.spectrum_shift)/Hartree;

  if (opt.particle == "el")
    st_shift_value -= 0.01/Hartree;
  else
    st_shift_value += 0.01/Hartree;

  return st_shift_value;

}

//=========================================================================//

void EnvelopFunctionApprox::read_SLEPC_solution(unsigned int number_of_ev )
{//
  /*
  1) Read all eigenvalues
  2) Sort the eigenvalues and select those we want
  3) Read eigenvectors that correspond to the eigenvalues we want
  4) normalize eigenfunctions
  5) calculate fermi energy for each state
  */


  // const short int_size = sizeof(int);

  //const short double_size = sizeof(double);

  // string fname_eigvals = "eigvals_SLEPC.out";



  // char buffer[int_size];
  // char buffer_double[double_size];
  //unsigned int dummy;

  //unsigned long long fict;


  //--------------------------------------------------------------------
  //how many solutions do we have from SLEPC?
  unsigned int number_of_converged_solutions;

  number_of_converged_solutions = EigenSolver::number_of_converged_eigenvalues();

  if (opt.log_output)
  {
    ostringstream os;
    os << "Number of converged solutions: "
      << number_of_converged_solutions;
    Messages::info(os.str());
  }



  //--------------------------------------------------------------------
  //read eigenvalues

  vector<EnvelopFunctionApprox::eigen_energy>   ev(number_of_converged_solutions);

  for (unsigned ind = 0; ind < number_of_converged_solutions; ind++)
  {

    ev[ind].energy =  EigenSolver::get_eigenvalue(ind) * Hartree + opt.spectrum_shift;


    ev[ind].global_number = ind;

    if (opt.log_output)
    {
      ostringstream os;
      os << ev[ind].global_number  << "    " << ev[ind].energy;
      Messages::info(os.str());
    }



  }

  //---------------------------------------------------------------------
  //sorting of the solutions


  if (opt.particle == "el") sort( ev.begin(), ev.end(), compare_eigen_energy_electrons1);

  if (opt.particle == "hl") sort( ev.begin(), ev.end(), compare_eigen_energy_holes1 );


  //----------------------------------------------------------------------
  //let us find the ground electron or hole level
  unsigned int ground_state_index = 0;
  bool finish = false;



  for (unsigned int i = 0; (i < number_of_converged_solutions && (!finish) ); i++)
  {
    if (opt.particle == "el")
    {
      if (ev[i].energy > opt.spectrum_shift)
      {
	ground_state_index = i;
	finish = true;
      }
    }
    else
    {
      if (ev[i].energy < opt.spectrum_shift)
      {
	ground_state_index = i;
	finish = true;
      }
    }



  }


  if (!finish)
  {
    throw ModelErrorException("EnvelopFunctionApprox: ground state is not found. Correct spectrum_shift or increase number of states");
  }

  unsigned int solution_size;
  if (number_of_converged_solutions - ground_state_index < number_of_ev)
    solution_size = number_of_converged_solutions - ground_state_index;
  else
    solution_size = number_of_ev;


  {
    EnvelopFunctionApprox::eigen_propblem_solution temp1;
    temp1.eigen_energy = 0;
    temp1.Fermi_energy = 0;
    temp1.Temperature = opt.Temperature;
    temp1.eigen_vector.resize(number_of_all_dofs, Complex(0.0, 0.0));

    solution.clear();
    solution.resize(solution_size, temp1);
  }













  map<unsigned int, unsigned int>  global_to_sol_index;
  map<unsigned int, unsigned int>  :: iterator it;


  for (unsigned int i = ground_state_index; i < ground_state_index + solution_size ; i++)
  {
    global_to_sol_index.insert( make_pair(ev[i].global_number, i - ground_state_index )  );
    solution[i - ground_state_index].eigen_energy = ev[i].energy;
  }



  //--------------------------------------------------------------------
  //read eigenvectors



  //read solutions - only independent dofs


  //----------------------------------------------------------------------
  for (unsigned int ind = 0; ind < number_of_converged_solutions; ind++)
    {

      vector<Complex> temp;

      EigenSolver::get_eigen_vector( ind, temp);

      it = global_to_sol_index.find(ind);

      if (  it  !=  global_to_sol_index.end() )
      {
	unsigned int solution_number = it->second;


	//-----------------------------------------------------------------------------
	//put independent dofs in the eigenvectors that may contain also non independent dofs
	for (unsigned j = 0; j < number_of_all_dofs; j++)
	{
	  if (new_dofs[j].independent)
	  {

	    solution[solution_number].eigen_vector[j] = temp[new_dofs[j].new_number];


	  }
	}


	//put constrained dofs




	for (unsigned int j = 0; j < number_of_all_dofs; j++)
	{



	  DofConstraints :: iterator it;

	  it = my_dof_constraints.find(j);


	  if (it != my_dof_constraints.end() )
	  {

	    DofConstraintRow constr_row = it->second;

	    DofConstraintRow::iterator  c =  constr_row.begin();


	    for ( ; c != constr_row.end() ; ++c )
	    {


	      solution[solution_number].eigen_vector[j] += ( c->second ) * solution[solution_number].eigen_vector[(c->first)];
	    }



	  }


	}



      }

      //------------------------------------------------------------------------

    }



  //normalization
  for (unsigned int i = 0; i < solution_size; i++)
    {
      const double norm = eigenstate_norm(i);


      const unsigned int n1 =  solution[i].eigen_vector.size();


      for (unsigned int j = 0; j < n1; j++)
	solution[i].eigen_vector[j] /= Complex(norm, 0.0);
    }





  //Fermi energy calculation

  if (poisson_equation != NULL)
    for (unsigned int i = 0; i < solution_size; i++)
    {
      solution[i].Fermi_energy = calculate_fermi_averaged( i);

    }


  //Temperature calculation
  if (temperature_simulation != NULL)
    for (unsigned int i = 0; i < solution_size; i++)
    {
      solution[i].Temperature = calculate_temperature_averaged(i);

    }




}



//=======================================================================//

void EnvelopFunctionApprox::output_eigen_function(unsigned int state_number,  const std::string& filename)
{
  //===========================================

  DofMap& dof_map = system->get_dof_map();

  const MeshBase& mesh1 = system->get_mesh();


  MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
  const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();


  unsigned int number_of_points = 0;

  for ( ; nd != nd_el ; ++nd)  number_of_points++;




  //vector of names
  vector<string> output_names(psi_name.size() * 2);


  unsigned int number_output_data = output_names.size() ;

  for (short i = 0; i < number_output_data/2; i++)
    {
      output_names[i*2    ] = psi_name[i] + "_r";
      output_names[i*2 + 1] = psi_name[i] + "_i";
    }




  //vector
  unsigned int  point_index = 0;

  unsigned int output_size = mesh->n_nodes()  * opt.number_of_bands * 2;

  vector<double>  psi_data(output_size);


  MeshBase::const_element_iterator it = mesh1.active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh1.active_local_elements_end();

  std::vector<unsigned int> dof_indices;



  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;


      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  for (unsigned int n = 0; n < elem->n_nodes(); n++)
	    {
	      unsigned int  node_id =  elem->node(n);
	      Complex value =  (solution[state_number].eigen_vector[ dof_indices[n] ]);

	      psi_data[psi_index*2 + node_id*number_output_data] = value.real();

	      psi_data[psi_index*2 + 1  + node_id*number_output_data] = value.imag();



	    }
	}
    }





  //std :: cout << filename << "\n";


  if (opt.output_type == "GMV")
    GMVIO(mesh1).write_nodal_data(filename, psi_data, output_names);

  if (opt.output_type == "tecplot")
    TecplotIO(mesh1,false).write_nodal_data(filename,psi_data,output_names);

  if (opt.output_type == "gnuplot" || opt.output_type == "GNUplot")
    GnuPlotIO(mesh1).write_nodal_data(filename, psi_data, output_names);


}

//=======================================================================//

void  EnvelopFunctionApprox::prepare_probability_function(const unsigned int state_number, std::vector<double>& prob_data)
{

  DofMap& dof_map = system->get_dof_map();

  const MeshBase& mesh1 = system->get_mesh();


  MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
  const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();


  unsigned int number_of_points = 0;

  for ( ; nd != nd_el ; ++nd)  number_of_points++;





  prob_data.clear();

  prob_data.resize(number_of_points, 0.0);


  vector< vector<Complex> >  psi_data;
  psi_data.resize(number_of_points);
  for (unsigned int i = 0; i < number_of_points; i++) psi_data[i].resize( opt.number_of_bands, Complex(0.0, 0.0) );


  MeshBase::const_element_iterator it = mesh1.active_elements_begin();
  const MeshBase::const_element_iterator end =  mesh1.active_elements_end();

  std::vector<unsigned int> dof_indices;


  //!calculation of probability function
  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;


      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  for (unsigned int n = 0; n < elem->n_nodes(); n++)
	    {
	      unsigned int  node_id =  elem->node(n);

	      Complex value =  (solution[state_number].eigen_vector[ dof_indices[n] ]);

	      psi_data[node_id][psi_index] = value;



	    }
	}
    }

  //done
  //calculation of |psi|^2
  double t1;
  for (unsigned int i = 0; i < number_of_points; i++)
  {
    //prob_data[i] = 0;
    for (unsigned int j = 0; j < opt.number_of_bands; j++)
    {

      t1 = std::abs(psi_data[i][j]);
      prob_data[i] += t1 * t1;
    }
    // if (state_number == 0) cerr << prob_data[i] << "\n";
  }
}




void EnvelopFunctionApprox::output_probability_function(unsigned int state_number,  const std::string& filename)
{
  //===========================================

  DofMap& dof_map = system->get_dof_map();

  const MeshBase& mesh1 = system->get_mesh();

  MeshBase::const_node_iterator       nd     = mesh1.active_nodes_begin();
  const MeshBase::const_node_iterator nd_el  = mesh1.active_nodes_end();


  unsigned int number_of_points = 0;

  for ( ; nd != nd_el ; ++nd)  number_of_points++;




  //vector of names
  vector<string> output_names(1);
  output_names[0] = "|psi|^2";


  unsigned int  point_index = 0;



  vector< vector<Complex> >  psi_data;
  psi_data.resize(number_of_points);
  for (unsigned int i = 0; i < number_of_points; i++) psi_data[i].resize( opt.number_of_bands, Complex(0.0, 0.0) );

  vector<double>  probability_data(number_of_points, 0.0);

  MeshBase::const_element_iterator it = mesh1.active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh1.active_local_elements_end();

  std::vector<unsigned int> dof_indices;


  //!calculation of psi

  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;


      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  for (unsigned int n = 0; n < elem->n_nodes(); n++)
	    {
	      unsigned int  node_id =  elem->node(n);
	      Complex value =  (solution[state_number].eigen_vector[ dof_indices[n] ]);

	      psi_data[node_id][psi_index] = value;





	    }
	}
    }


  //calculation of |psi|^2
  double t1;
  for (unsigned int i = 0; i < number_of_points; i++)
    for (unsigned int j = 0; j < opt.number_of_bands; j++)
      {
	t1 = std::abs(psi_data[i][j]);
	probability_data[i] += t1 * t1;
      }


  //std :: cout << filename << "\n";



  if (opt.output_type == "GMV")
    GMVIO(mesh1).write_nodal_data(filename, probability_data, output_names);

  if (opt.output_type == "tecplot")
    TecplotIO(mesh1,false).write_nodal_data(filename, probability_data,output_names);

  if (opt.output_type == "gnuplot" || opt.output_type == "GNUplot")
    GnuPlotIO(mesh1).write_nodal_data(filename, probability_data,output_names);

}



//=======================================================================//


//=======================================================================//





//=======================================================================//
/*
void EnvelopFunctionApprox::make_new_dofs( )
{
  new_dofs.clear();


  DofMap& dof_map = system->get_dof_map();


  number_of_all_dofs  =   dof_map.n_dofs();

  new_dofs.resize(number_of_all_dofs);
  const std::set<unsigned int> :: const_iterator  n_begin = dirichlet_dofs.begin();
  const std::set<unsigned int> :: const_iterator  n_end   = dirichlet_dofs.end();
  std::set<unsigned int> :: const_iterator n_it;

  unsigned int number_it = 0;


  for (unsigned int i = 0; i < number_of_all_dofs ; i++)
    {
      if ( !( dof_map.is_constrained_dof(i) ) && (find(n_begin, n_end, i) == n_end))
	{
	    new_dofs[i].independent = true;
	    new_dofs[i].new_number = number_it;
	    number_it++;

	  }
	else
	  {

	    new_dofs[i].independent = false;
	  }
    }


  number_of_new_dofs = number_it;



}
*/
//=======================================================================//
bool EnvelopFunctionApprox::compare_eigen_energy_electrons(eigen_propblem_solution state1, eigen_propblem_solution state2)
{
  return(state1.eigen_energy < state2.eigen_energy );
}


//=======================================================================//
bool EnvelopFunctionApprox::compare_eigen_energy_holes(eigen_propblem_solution state1, eigen_propblem_solution state2)
{
  return(state1.eigen_energy > state2.eigen_energy );
}

//=======================================================================//

bool EnvelopFunctionApprox::compare_eigen_energy_holes1(eigen_energy state1, eigen_energy state2)
{
  return(state1.energy> state2.energy );
}

//=======================================================================//
bool EnvelopFunctionApprox::compare_eigen_energy_electrons1(eigen_energy state1, eigen_energy state2)
{
  return(state1.energy<  state2.energy );
}


//=======================================================================//
EnvelopFunctionApprox:: ~EnvelopFunctionApprox(void)
{
  // es->delete_system(system_name);
}

//=======================================================================//


//-----------------------------------------------------------------------------//
double  EnvelopFunctionApprox::eigenstate_norm(unsigned int state_number)
{
  double  result;

  const vector< Complex > &  eigen_vector =  solution[state_number].eigen_vector;


  DofMap& dof_map = system->get_dof_map();

  FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

  // AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  QGauss qrule (dim, SECOND);

  fe -> attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  std::vector<unsigned int> dof_indices;


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();







  Complex temp(0.0, 0.0);
  Complex eigen_f_value1, eigen_f_value2;

  for ( ; el != end_el ; ++el)
    {//el

      const Elem* elem = *el;
      fe->reinit (elem);

      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  const unsigned int n_psi_dofs = dof_indices.size();

	  for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	    {//qp
	      for (unsigned int p1=0; p1<n_psi_dofs; p1++)
		{
		  eigen_f_value1 = eigen_vector[dof_indices[p1]];
		  for (unsigned int p2=0; p2<n_psi_dofs; p2++)
		    {
		      eigen_f_value2 = eigen_vector[dof_indices[p2]];
		      temp += ( JxW[qp] * phi[p1][qp] * eigen_f_value1 *  phi[p2][qp] * conj(eigen_f_value2) );
		    }
		}

	    }




	}


    }



  result = sqrt( abs(temp)  );




  return(result);

}



//==========================================================//


double EnvelopFunctionApprox::calculate_fermi_averaged(unsigned int i)
{

  Complex  result(0.0,0.0);


  //-----------------------------------------------------//





  const vector< Complex >&   eigen_vector =  solution[i].eigen_vector;




  //----------------------------------------------------//


  const MeshBase* mesh = &(es->get_mesh());


  unsigned int dim = mesh->mesh_dimension();



  system = &( es->get_system<LinearImplicitSystem>(system_name));

  DofMap& dof_map = system->get_dof_map();





  //My Jacobian







   FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

   //  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
   AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

   // A 5th order Gauss quadrature rule for numerical integration.
   QGauss qrule (dim, SECOND);

   // Tell the finite element object to use our quadrature rule.
   fe -> attach_quadrature_rule (&qrule);

   // The element Jacobian * quadrature weight at each integration point.
   const std::vector<Real>& JxW = fe->get_JxW();

   // properties at the quadrature points.
   const std::vector<Point>& q_point = fe->get_xyz();

   // The element shape functions evaluated at the quadrature points.
   const std::vector<std::vector<Real> >& phi = fe->get_phi();


   //------------------------------------------------------------
   std::vector<unsigned int> dof_indices_component;

   std::vector<unsigned int> dof_indices;

   //-------------------------------------------------------------

  //----------------------------------------------------//


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


  Complex eigen_f_value1;
  Complex eigen_f_value2;


  double chem_pot_value_eV;

  for ( ; el != end_el ; ++el)
    {//el

      const Elem* elem = *el;
      fe->reinit (elem);


      Point center = elem->centroid();

      chem_pot_value_eV = get_electro_chem_potential(elem);



      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  const unsigned int n_psi_dofs = dof_indices.size();

	  for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	    {//qp

	      for (unsigned int p1=0; p1<n_psi_dofs; p1++)
		{
		  eigen_f_value1 = eigen_vector[dof_indices[p1]];
		  for (unsigned int p2=0; p2<n_psi_dofs; p2++)

		    {
		      eigen_f_value2 = eigen_vector[dof_indices[p2]];
		      result += ( JxW[qp] * phi[p1][qp] * eigen_f_value1 *  phi[p2][qp] * conj(eigen_f_value2) ) * chem_pot_value_eV;

		    }
		}

	    }




	}


    }








  return(result.real());



}

//----//==========================================================//


double EnvelopFunctionApprox::calculate_temperature_averaged(unsigned int i)
{

  Complex  result(0.0,0.0);


  //-----------------------------------------------------//





  const vector< Complex >&   eigen_vector =  solution[i].eigen_vector;




  //----------------------------------------------------//


  const MeshBase* mesh = &(es->get_mesh());


  unsigned int dim = mesh->mesh_dimension();



  system = &( es->get_system<LinearImplicitSystem>(system_name));

  DofMap& dof_map = system->get_dof_map();





  //My Jacobian







   FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

   //  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
   AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

   // A 5th order Gauss quadrature rule for numerical integration.
   QGauss qrule (dim, SECOND);

   // Tell the finite element object to use our quadrature rule.
   fe -> attach_quadrature_rule (&qrule);

   // The element Jacobian * quadrature weight at each integration point.
   const std::vector<Real>& JxW = fe->get_JxW();

   // properties at the quadrature points.
   const std::vector<Point>& q_point = fe->get_xyz();

   // The element shape functions evaluated at the quadrature points.
   const std::vector<std::vector<Real> >& phi = fe->get_phi();


   //------------------------------------------------------------
   std::vector<unsigned int> dof_indices_component;

   std::vector<unsigned int> dof_indices;

   //-------------------------------------------------------------

  //----------------------------------------------------//


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


  Complex eigen_f_value1;
  Complex eigen_f_value2;



  double Temperature;
  vector<double> values;
  vector<Point> qp(1);

  for ( ; el != end_el ; ++el)
    {//el

      const Elem* elem = *el;
      fe->reinit (elem);


      Point center = elem->centroid();
      qp[0] = center;

      temperature_simulation->get_solution(elem, qp, temperature_ID, values);

      Temperature = values[0];




      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  const unsigned int n_psi_dofs = dof_indices.size();

	  for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	    {//qp

	      for (unsigned int p1=0; p1<n_psi_dofs; p1++)
		{
		  eigen_f_value1 = eigen_vector[dof_indices[p1]];
		  for (unsigned int p2=0; p2<n_psi_dofs; p2++)

		    {
		      eigen_f_value2 = eigen_vector[dof_indices[p2]];
		      result += ( JxW[qp] * phi[p1][qp] * eigen_f_value1 *  phi[p2][qp] * conj(eigen_f_value2) ) * Temperature;

		    }
		}

	    }




	}


    }








  return(result.real());



}

//--------------------------------------------------------------------------//


//------------------------------------------------------------------------------//

void EnvelopFunctionApprox::calculate_density( )
{
  DofMap& dof_map = system->get_dof_map();

  const MeshBase& mesh1 = system->get_mesh();


  _density.clear();


  vector<double> density_of_state;

  unsigned int number_of_eigenfunctions = opt.initial_eigenstates_number;
  //unsigned int number_of_eigenfunctions = solution.size();



  // we might specify the lowest state
  for (unsigned int i = opt.first_state; i < number_of_eigenfunctions; i++)
  {



    const double Fermi_energy = solution[i].Fermi_energy;

    const double Energy = solution[i].eigen_energy;

    double prob_factor = Fermi_statistics_probability(Energy, Fermi_energy, solution[i].Temperature); //Thermal probability


    density_of_state = calculate_cell_prob_function(i);


    const MeshBase::const_element_iterator it_end  = mesh1.active_elements_end();
    const MeshBase::const_element_iterator it_begin  = mesh1.active_elements_begin();

    MeshBase::const_element_iterator       it     =  it_begin;

    unsigned int el_number = 0;

    for( ;it !=  it_end ;++it)
    {
      const Elem* el = *it;



      if (opt.local_occupation)
      {

	Point center = el->centroid();

	double chem_pot_value_eV = get_electro_chem_potential(el);


	prob_factor = Fermi_statistics_probability(Energy,  chem_pot_value_eV, solution[i].Temperature); //Thermal probability
      }


      double temp = density_of_state[el_number] * prob_factor;

      if (it == it_begin && i == 0)
      {

	_density.insert(pair<const Elem*, double> (el, temp));
      }
      else
	_density[el] += temp;



      el_number++;
      }




    }




}

//===========================================================//
vector<double>  EnvelopFunctionApprox::calculate_cell_prob_function(unsigned int state_number)
{

  const vector< Complex > &  eigen_vector =  solution[state_number].eigen_vector;

  DofMap& dof_map = system->get_dof_map();

  const MeshBase& mesh = system->get_mesh();

  FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

  // AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  QGauss qrule (dim, FIFTH);

  fe -> attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  std::vector<unsigned int> dof_indices;


  MeshBase::const_element_iterator       el1     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el1 = mesh.active_elements_end();

  //--------------------------------------------------------------//
  //define results
  unsigned int active_el_number = 0;

  for ( ; el1 != end_el1; ++el1) active_el_number++;

  vector<double> result(active_el_number, 0);
  vector<Complex> result_complex(active_el_number, Complex(0.0, 0.0));
  //--------------------------------------------------------------//


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();






  Complex eigen_f_value1, eigen_f_value2;
  unsigned int el_number = 0;

  for ( ; el != end_el ; ++el)
    {//el

      double el_volume = 0;

      const Elem* elem = *el;
      fe->reinit (elem);



      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
      {
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  const unsigned int n_psi_dofs = dof_indices.size();

	  for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	    {//qp


	      for (unsigned int p1=0; p1<n_psi_dofs; p1++)
		{
		  eigen_f_value1 = eigen_vector[dof_indices[p1]];

		  for (unsigned int p2=0; p2<n_psi_dofs; p2++)
		    {
		      eigen_f_value2 = eigen_vector[dof_indices[p2]];

		      Complex temp =  eigen_f_value1 * conj(eigen_f_value2);

		      // result[el_number] +=  JxW[qp] * phi[p1][qp] *  phi[p2][qp] * temp.real();
		      result_complex[el_number] +=  JxW[qp] * phi[p1][qp] *  phi[p2][qp] * temp;
		    }

		}



	    }




	}

      //-------------------------------------------------//
      //has to be removed for the new libmesh

      for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
	el_volume += JxW[qp];

      result[el_number] = std::abs(result_complex[el_number])/el_volume;

      el_number++;

    }







  return(result);

}




//===========================================================//



vector<double> EnvelopFunctionApprox::calculate_prob_function(unsigned int state_number)
{
//===========================================



  DofMap& dof_map = system->get_dof_map();

  const MeshBase& mesh1 = system->get_mesh();


  unsigned int number_of_points = number_of_nodes;



  unsigned int  point_index = 0;


  vector< vector<Complex> >  psi_data;
  psi_data.resize(number_of_points);
  for (unsigned int i = 0; i < number_of_points; i++) psi_data[i].resize( opt.number_of_bands, Complex(0.0, 0.0) );

  vector<double>  probability_data(number_of_points, 0.0);

  MeshBase::const_element_iterator it = mesh1.active_local_elements_begin();
  const MeshBase::const_element_iterator end =  mesh1.active_local_elements_end();

  std::vector<unsigned int> dof_indices;


  //!calculation of psi

  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;


      for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
	{
	  dof_map.dof_indices (elem, dof_indices, psi_index);
	  for (unsigned int n = 0; n < elem->n_nodes(); n++)
	    {
	      unsigned int  node_id =  elem->node(n);
	      Complex value =  (solution[state_number].eigen_vector[ dof_indices[n] ]);

	      psi_data[node_id][psi_index] = value;





	    }
	}
    }


  //calculation of |psi|^2
  double t1;
  for (unsigned int i = 0; i < number_of_points; i++)
    for (unsigned int j = 0; j < opt.number_of_bands; j++)
      {
	t1 = std::abs(psi_data[i][j]);
	probability_data[i] += t1 * t1;
      }



  //------

  return(probability_data);

  //-----

}


//=================================================================//

double EnvelopFunctionApprox::get_integrated_probability()
{
  double result = 0;
  unsigned int number_of_eigs = solution.size();
  for (unsigned int i = opt.first_state ; i <  number_of_eigs; i++)
    {

      result += Fermi_statistics_probability(solution[i].eigen_energy, solution[i].Fermi_energy, solution[i].Temperature);
    }

  return(result);
}


//=================================================================//
void EnvelopFunctionApprox::calculate_convergent_density( )
{

  if (opt.job == BULKDENSITY)
  {
    //temporary solution

    solve_bulk();

    double total_density =  get_integrated_probability();

    _density.clear();

    _density[NULL] = total_density;




  }
  else
  {
    unsigned int number_of_states = opt.initial_eigenstates_number;

    solve_eigen_value_problem(number_of_states);



    double total_density =  get_integrated_probability();

    bool converged = false;


    for  ( ; opt.convergent_density && (!converged); )
    {


      number_of_states = (unsigned int) (number_of_states * opt.eigen_number_increase_factor) + 1;



      //TODO:

      // this is to think about why it does not work!

      //opt.solve_ev_problem_twice = false;

      //double st_shift_value = (solution[0].eigen_energy - opt.spectrum_shift)/Hartree;

      //solve_eigen_value_problem(number_of_states, st_shift_value);



      solve_eigen_value_problem(number_of_states);

      double total_density1 = get_integrated_probability();


      if (total_density1 < 1e-14)
      {
        converged = true;
      }
      else if ( abs(total_density1 - total_density)/total_density < opt.relative_density_tolerance )
      {
        converged = true;
      }

      total_density = total_density1;


    }


    calculate_density( );

  }
}

//==============================================================================//

unsigned int EnvelopFunctionApprox::get_number_of_active_cells()
{
  unsigned int result = 0;
  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  for ( ; el != end_el ; ++el)
    result++;


  return(result);



}

//===============================================================================//
short EnvelopFunctionApprox::calculate_number_of_bands(void) const
{

  short result = 0;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  EFAbulkHamiltonian* element_hamiltonian;

  const Elem* elem = *el;
  const ID subdomain = elem->subdomain_id();
  const Material* mat = _device->get_material(subdomain);


  element_hamiltonian = (  dynamic_cast<EFAbulkModel*> ( mat ->get_model(get_id()) )    ) -> get_Hamiltonian_model() ;

  std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
	    model_Ham = ( element_hamiltonian->get_Hamiltonian() );

  result = model_Ham.size();

  if (result == 0)  throw InitFailedException("EnvelopFunctionApprox: Hamiltonian is empty");

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    const ID subdomain = elem->subdomain_id();
    const Material* mat = _device->get_material(subdomain);


    element_hamiltonian = ( dynamic_cast<EFAbulkModel*> (  mat ->get_model(get_id()) ) )->get_Hamiltonian_model() ;

    const std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
      model_Ham = ( element_hamiltonian->get_Hamiltonian() );


    short result1 = model_Ham.size();

    if (result1 != result) throw InitFailedException("EnvelopFunctionApprox: Hamiltonians of different materials have different number of bands");


  }

  return(result);

}


//========================================================================================//

std::map<const Elem*, double> EnvelopFunctionApprox::estimate_density1D(unsigned int state_number, double parallel_mass)
{


  vector<double> result;
  map<const Elem*, double> result_map;





  const double T_EV = opt.Temperature * Constants::k_Boltzmann;

  const double mass_factor = 0.5 * parallel_mass/M_PI * (T_EV /Constants::Hartree);

  const double Fermi_energy = solution[state_number].Fermi_energy;

  const double Energy = solution[state_number].eigen_energy;



  double exp_arg = (Fermi_energy - Energy) / T_EV;

  if (opt.particle == "hl")
    exp_arg = -exp_arg;


  double prob_factor = (exp_arg < -20) ? exp(exp_arg) : std::log(1.0 + exp(exp_arg));



  result  = calculate_cell_prob_function(state_number);


  unsigned int n = 0;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


  for ( ; el !=end_el; ++el)
  {

    result_map[*el] = result[n] *  prob_factor * mass_factor  ;

    n++;
  }


  return(result_map);

}

//========================================================================================//

std::map<const Elem*, double> EnvelopFunctionApprox::estimate_density2D(unsigned int state_number, double parallel_mass)
{


  vector<double> result;

  map<const Elem*, double> result_map;

  const double T_EV = opt.Temperature * Constants::k_Boltzmann;

  const double Fermi_energy = solution[state_number].Fermi_energy;

  const double Energy = solution[state_number].eigen_energy;

  double prob_factor;

  if (opt.particle == "el")
  {
    prob_factor = TiberMath::fermidirac_mhalf( (Fermi_energy - Energy)  / T_EV);
  }
  else
  {
    prob_factor = TiberMath::fermidirac_mhalf( -(Fermi_energy - Energy)  / T_EV);
  }



  const double mass_factor = 0.5*sqrt( T_EV *  parallel_mass / (2.0*M_PI) );

  result  = calculate_cell_prob_function(state_number);


  unsigned int n = 0;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


  for (  ; el !=end_el; ++el )
  {

    result_map[*el] = result[n] *  prob_factor * mass_factor ;

    n++;
  }

  return(result_map);
}

//=======================================================================================//

void EnvelopFunctionApprox::solve_bulk(void)
{
  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


  bool found = false;
  const Elem* mat_elem;

  for ( ; (el != end_el) && (!found) ; ++el)
  {
    const Elem* elem = *el;
    if (elem->contains_point(_bulk_point))
    {
      found = true;
      mat_elem = elem;
      _bulk_mat_element = mat_elem;
    }

  }

  if (!found) throw SolveFailedException("Bad material point\n");

  Point qp = mat_elem->centroid();

  EFAbulkHamiltonian* element_hamiltonian;

  const ID subdomain = mat_elem->subdomain_id();

  const Material* mat = _device->get_material(subdomain);

  element_hamiltonian = (  dynamic_cast<EFAbulkModel*> (  mat ->get_model(get_id()) )  )->get_Hamiltonian_model();

  element_hamiltonian->set_k_vector(k_vector);

  element_hamiltonian->calculate_Hamiltonian_k_par();

  Tensor2Sym strain_crystal_system(0);

  double electric_potential = 0;

  if (opt.consider_strain)
  {

    strain_crystal_system = strain->get_strain_crystal(mat_elem , qp);
  }

  if (opt.consider_potential)
  {


    electric_potential = get_electric_potential( mat_elem, qp );
  }




  element_hamiltonian->apply_strain_and_potential(strain_crystal_system, electric_potential);


  std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
	    model_Ham = ( element_hamiltonian->get_Hamiltonian() );


  std::complex<double> ham_matrix[opt.number_of_bands * opt.number_of_bands ];

  for (unsigned int band1 = 0; band1 < opt.number_of_bands; band1++)
    for (unsigned int band2 = 0; band2 < opt.number_of_bands; band2++)
    {
      ham_matrix[band1 + band2 * opt.number_of_bands] = model_Ham[band1][band2].constant;


    }


  solution.clear();
  solution.resize(opt.number_of_bands);

  {
     char jobs = 'V';
     char UPLO = 'U';
     int  N = opt.number_of_bands;
     double eigvals[opt.number_of_bands];
     std::complex<double> WORK[2*opt.number_of_bands-1];
     int LWORK = 2*opt.number_of_bands-1;
     double RWORK[3*opt.number_of_bands-2];
     int info;

     zheev_(jobs, UPLO, N, ham_matrix, N, eigvals, WORK, LWORK, RWORK, info);

     if (info !=0 ) throw SolveFailedException("LAPACK problem\n");;

     for (short i = 0; i < opt.number_of_bands ; i++)
     {

       solution[i].eigen_energy = eigvals[i]*Hartree;
       solution[i].eigen_vector.resize(opt.number_of_bands);
       for (short j = 0; j < opt.number_of_bands ; j++)
       {
	 solution[i].eigen_vector[j] = ham_matrix[i * opt.number_of_bands + j];
       }
     }
  }






  unsigned int n = solution.size();

  for (unsigned int i = 0; i < n ; i++)
  {
    solution[i].Fermi_energy = 0;
    solution[i].Temperature = opt.Temperature;


    if (poisson_equation != NULL)
    	solution[i].Fermi_energy = get_electro_chem_potential(mat_elem);



    //Temperature calculation

    if (temperature_simulation != NULL)
    {

      Point center = mat_elem->centroid();

      vector<double> values;
      vector<Point> qp(1);

      qp[0] = center;

      temperature_simulation->get_solution(mat_elem, qp, temperature_ID, values);

      solution[i].Temperature = values[0];

    }

  }


}
