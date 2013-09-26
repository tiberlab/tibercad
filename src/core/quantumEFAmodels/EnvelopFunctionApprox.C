// $Id$

#include "SimulationEnvironment.h"
#include "EnvelopFunctionApprox.h"
#include "ModelOptions.h"
#include "EFAbulkModel.h"
#include "Material.h"
#include "Boundary.h"
#include "TiberMath.h"
#include "TiberLinearSystem.h"
#include <gnuplot_io.h>
#include "SimulationOptions.h"
#include "tensor.h"

#include "EigenSolver.h"


#include <edge_edge2.h>
#include <equation_systems.h>
#include <dense_submatrix.h>


#include "Messages.h"



extern "C"
{
   void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& LDA,
	       double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info);
};

using namespace std;
using namespace Constants;


//---------------------------------------------------------------------------------//


/*
inline void EnvelopFunctionApprox::get_electric_potential(const Elem* elem, const std::vector<Point>& q_point,
						   std::vector<double> electric_potential) const
{
  poisson_equation->get_solution (elem, q_point, potential_ID, electric_potential);
}
*/
//---------------------------------------------------------------------------------//

inline double EnvelopFunctionApprox::get_band_edge(const Elem* elem, const std::string& particle) const
{
  vector<double> values(elem->n_nodes());
  vector<Point> p(elem->n_nodes());
  
  for (size_t i = 0; i < elem->n_nodes(); ++i)
        p[i] = elem->point(i);

  poisson_equation->get_solution(elem, band_edge_ID, values, p);

  double bedge = values[0];

  for (size_t i = 1; i < elem->n_nodes(); ++i)
  {
     double temp = values[i];
     if(particle == "el")
        bedge = (temp < bedge) ? temp : bedge;
     else
        bedge = (temp > bedge) ? temp : bedge;
   }

  return bedge;
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

void EnvelopFunctionApprox::get_eigenenergies(std::vector<double>& values) const
{


  unsigned int n = _solution.size();
  values.resize(n);
  for (unsigned int i = 0; i < n; i++)
  {
    values[i] = _solution[i].eigen_energy;
  }


}

//---------------------------------------------------------------------------------//

void EnvelopFunctionApprox::get_occupations(std::vector<double>& values) const
{


  unsigned int n = _solution.size();
  values.resize(n);

  for (unsigned int i = 0; i < n; i++)
  {
    values[i ] = Fermi_statistics_probability(_solution[i].eigen_energy, 
                                              _solution[i].electro_chem_pot, 
                                              _solution[i].temperature);
  }

}




//====================================================================================//

//const std::vector<eigen_problem_solution >& EnvelopFunctionApprox::get_solution() const
//{
//  return(solution);
//}


//====================================================//
PhysicalModel* EnvelopFunctionApprox::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{

  ModelOptions opts(options);
  opts.set_option("particle", get_option("particle", ""));

  EFAbulkModel* model = PhysicalModelInterface::create<EFAbulkModel>("EFAmodel", mat, opts);

  if (model == NULL)
    throw ModelErrorException("efaschroedinger: cannot create bulk model");

  return model;

}



void
EnvelopFunctionApprox::do_setup_solution_variables(void)
{
  // declare solution variables
  unsigned int dim = get_mesh().mesh_dimension();
  string units("1/cm");
  if (dim == 2)
    units = "1/cm^2";
  else if (dim == 3)
    units = "1/cm^3";
  declare_solution(ProbabilityDensity, NTUPLE, NODES, units);
  add_alias("EigenFunctions", ProbabilityDensity);
  declare_solution(EigenEnergy, NTUPLE, GLOBAL, "eV");
  declare_solution(Occupation, NTUPLE, GLOBAL, "");
  declare_solution(EigenEnergyOnMesh, NTUPLE, NODES, "eV");

  if (plot_solution(EigenEnergy) && plot_solution(ProbabilityDensity))
    add_plot_variable(EigenEnergyOnMesh);

  if (plot_solution("QuantumDensity") || get_options().has_submodel("QuantumDensity"))
    _calculate_density = true;
  if (_calculate_density) declare_solution(QuantumDensity, REAL, NODES, "1/cm^3");
}


void
EnvelopFunctionApprox::get_solution_secure(map<ID, vector<double> >& values)
{
  if (values.count(EigenEnergy))
  {
    // number of states
    const unsigned int num_states = _solution.size();
    for (unsigned int sn = 0; sn < num_states; sn++)
    {
      values[EigenEnergy][sn] = _solution[sn].eigen_energy;
    }
  }

  if (values.count(Occupation))
  {
    // number of states
    const unsigned int num_states = _solution.size();
    for (unsigned int sn = 0; sn < num_states; sn++)
      values[Occupation][sn] =
          Fermi_statistics_probability(_solution[sn].eigen_energy,
                                       _solution[sn].electro_chem_pot,
                                       _solution[sn].temperature);
  }
}


void
EnvelopFunctionApprox::get_solution_secure(const Elem* elem,
    map<ID, vector<double> >& values, const vector<Point>& points)
{

  unsigned int np = points.size();

  if (values.count(ProbabilityDensity))
  {
    unsigned int dim = get_mesh().mesh_dimension();
    double scale = Constants::bohr_radius * 1e2;
   if (dim > 1)
     scale *= Constants::bohr_radius * 1e2;
   if (dim > 2)
     scale *= Constants::bohr_radius * 1e2;

    FEType fe_type = system->variable_type(0);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
    const vector<vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &points);

    DofMap& dof_map = system->get_dof_map();
    std::vector<unsigned int> dof_indices;

    // number of states
    const unsigned int num_states = _solution.size();
    values[ProbabilityDensity] = vector<double>(np * num_states, 0.0);

    // they should all be the same size
    unsigned int n_dofs = phi.size();

    vector<double> prob_dens(np, 0.0);

    for (int psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
    {
      dof_map.dof_indices(elem, dof_indices, psi_index);

      for (unsigned int n = 0; n < np; n++)
      {
        for (unsigned int sn = 0; sn < num_states; sn++)
        {
          Complex value(0.0);

          // interpolate
          for (unsigned int i = 0; i < n_dofs; i++)
          {
            value += phi[i][n] * _solution[sn].eigen_vector[dof_indices[i]];
          }

          // calculate probability density
          double tmp = abs(value);

          // the number of components should be set already to the right number
          values[ProbabilityDensity][num_states * n + sn] += tmp * tmp / scale;
        }
      }
    }
  }

  if (values.count(EigenEnergyOnMesh))
  {
    // number of states
    const unsigned int num_states = _solution.size();

    for (unsigned int n = 0; n < np; n++)
      for (unsigned int sn = 0; sn < num_states; sn++)
        values[EigenEnergyOnMesh][num_states * n + sn] = _solution[sn].eigen_energy;
  }

  if (values.count(QuantumDensity))
  {
    TiberLinearSystem& qdens_sys = get_equation_system<TiberLinearSystem>(0);
    NumericVector<Number>& qdens = *qdens_sys.solution;

    FEType fe_type = qdens_sys.variable_type(0);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
    const vector<vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &points);

    DofMap& dof_map = qdens_sys.get_dof_map();
    std::vector<unsigned int> dof_indices;
    dof_map.dof_indices(elem, dof_indices, 0);
    unsigned int n_dofs = phi.size();

    for (unsigned int n = 0; n < np; n++)
    {
      double value = 0;

      // the phi^2 factor comes from the fact that the more correct interpolation is the square
      // of the basis function, because the probability densities are the square of the states
      // NOTE: maybe one should check if this gives really a better result
      // NOTE: 2011-12-01 the above turned out to be wrong: phi is used only as linear interp. !
      for (unsigned int i = 0; i < n_dofs; i++)
         value += phi[i][n] * qdens(dof_indices[i]);

      values[QuantumDensity][n] = value;
    }

  }
}





//====================================================//
double EnvelopFunctionApprox::get_band_edge(const std::string& particle) 
{

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  assert(el != end_el);

  double band_edge = get_band_edge(*el, particle);
  ++el;


  for (; el != end_el ; ++el )
  {
    const Elem* elem = *el;

    double temp = get_band_edge(elem, particle);

    if (particle == "el")
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

  return (band_edge);

}

//===================================================//
EnvelopFunctionApprox::EnvelopFunctionApprox(const ModelOptions& options)
 : FEMEigenvalueProblem(options),
   _calculate_density(false)
{
  poisson_equation = NULL;

  _bulk_mat_element = NULL;

  has_solution_vector(false);


}





//===================================================//
void EnvelopFunctionApprox::parse_options()
{

  FEMEigenvalueProblem::parse_options();

  // This is done to allow the user to set 'number_of_eigenstates' in the principal block
  // useful when Solver{} has only default values and is not defined
  if (has_option("number_of_eigenstates"))
  {
    solver_opt.number_of_eigenstates = get_option("number_of_eigenstates", 6);
  }


  // check the quadrature rule
  {
    string qrule = get_option("quadrature_rule", "gauss");
    if (qrule == "gauss")
      _quadrature_type = QGAUSS;
    else if (qrule == "trapez")
    {
      _quadrature_type = QTRAP;
      // this is not BIM, but the flag will make it solve a non-
      // generalized problem
      solver_opt.discretization_method = BIM;
    }
    else
      throw InitFailedException("Unknown quadrature rule");
  }
  //-------------------------------------------------------------------------------------------//

  std::string  job_name = get_option("job","eigenstates");
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

  if (opt.job == BULKEIGENSTATES || opt.job == BULKDENSITY)
  {
    if (has_option("bulk_point"))
    {
      vector<double> point;

      get_option("bulk_point",point);
      for (short i = 0; i < 3; i++)  _bulk_point(i) = point[i];
    }
    else
    {
      throw InitFailedException( "You have to specify a bulk_point");
    }
  }

  //--------------------------------------------------------------------------------//
  if (!has_option("particle"))
    throw InitFailedException( "EnvelopeFunctionApprox: option particle needed" );
    
  opt.particle  = get_option("particle","");
  
  if (opt.particle != "el" && opt.particle != "hl")
    throw InitFailedException( "EnvelopeFunctionApprox: particle must be 'el' or 'hl'" );    

  //Set degeneracy factor depending on model 

  // adjust the degeneracy (in a quirky way, I must admit...)
  //opt.degeneracy = 1; This has been already initialized in do_init()

  if (get_options().has_submodel("Physics"))
  {
    ModelOptions::const_submodel_iterator it(get_options().submodels_begin("Physics"));
    const ModelOptions& opts = it->second;

    string model = opts.get_option("model", "kp");
    if ((model == "conduction_band") || (model == "single_band"))  opt.degeneracy *= 2;
  }

  //possible user override
  opt.degeneracy = get_option("degeneracy", opt.degeneracy); 

  // for degeneracy = 1 we assure that number of states is even,
  // so we take both spin states
  if ((opt.degeneracy == 1) && (solver_opt.number_of_eigenstates % 2 == 1))
  {
    solver_opt.number_of_eigenstates += 1;
    Messages::warning("Number of eigenstates increased by 1 because of spin pairing");
  }


  //-------------------------------------------------------------------------------------------//
  //Strain model
  std::string strain_model_name = get_option("strain_model_name","");
  strain_model_name = get_option("strain_simulation", strain_model_name);
  _strain_interface.set_simulation(strain_model_name);


  //-------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------//
  //Poisson model
  opt.consider_potential = false;
  opt.consider_potential_bulk = false;

  std::string  poisson_model_name = get_option("poisson_model_name","");
  poisson_model_name = get_option("poisson_simulation", poisson_model_name);
  if ( poisson_model_name != "" )
  {
    opt.consider_potential = true;

    opt.consider_potential_bulk = get_option("potential_in_bulk",true);

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
    throw InitFailedException("Needed poisson_simulation");
    opt.consider_potential = false;
  }
  //---------------------------------------------------------------------------------//
  //Heat model
  std::string heat_model_name = get_option("heat_model","");
  heat_model_name = get_option("temperature_simulation", heat_model_name);

  _temp_interface.set_simulation(heat_model_name);


  //--------------------------------------------------------------------------------------------//
  //Spectrum Shift
  //as default, we  estimate spectrum shift only in electric potential is defined
  const ModelOptions& sol_opt = get_solver_options();

  opt.estimate_spectrum_shift =  opt.consider_potential;

  //opt.estimate_spectrum_shift = sol_opt.get_option("estimate_guess",  opt.estimate_spectrum_shift);

  if (sol_opt.find_option("guess")) opt.estimate_spectrum_shift = false;

  if ( !opt.consider_potential && opt.estimate_spectrum_shift) 
    throw InitFailedException( "EnvelopeFunctionApprox: cannot estimate guess without electric potential");

  if (!sol_opt.find_option("guess") && !opt.estimate_spectrum_shift)
    throw InitFailedException( "EnvelopeFunctionApprox: value for guess required");


  //--------------------------------------------------------------------------------------------//
  // k-vector
  if (has_option("k_vector"))
  {

    RealVectorValue k_vec(3,0.0);
    get_parameter("k_vector", k_vec);
    set_k_point(k_vec);

    Messages::warning("k-vector given, will skip density calculation.");
    _calculate_density = false;
  }

  //---------------------------------------------------------------------------------//
  // Options for converged density  (NOT USED NOW)
  //---------------------------------------------------------------------------------//
  opt.convergent_density = false; //get_option("convergent_density", true);
  opt.initial_eigenstates_number = sol_opt.get_option("initial_eigenstates_number", 
                                                     solver_opt.number_of_eigenstates );

  opt.eigen_number_increase_factor = get_option("eigen_number_increase_factor",1.2);

  opt.relative_density_tolerance =  get_option("relative_density_tolerance", 1e-2);


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

  opt.local_occupation = get_option("local_occupation", true);

  //--------------------------------------------------------------------------------------------//
  // Block QuantumDensity  //
  opt.first_state = 0; 
  opt.k_val = 0.01;
  opt.assume_paraboloid = false;

  if (get_options().has_submodel("QuantumDensity"))
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("QuantumDensity"));
    ModelOptions& opts = it->second;

    opt.first_state = opts.get_option("first_state", opt.first_state);
    opt.k_val = opts.get_option("k_value", opt.k_val);
    opt.assume_paraboloid = opts.get_option("assume_diagonal_mass_matrix", opt.assume_paraboloid);
  }


  

}



//===================================================//
void EnvelopFunctionApprox::do_init( )
{

  FEMEigenvalueProblem::do_init();

  es = &(get_equation_systems());
  
  mesh = &(es->get_mesh());
  
  system_name = get_equation_system_name ( );
  
  es->add_system<LinearImplicitSystem> (system_name);
  
  system = &( es->get_system<LinearImplicitSystem>( system_name ) );

  dim = mesh->mesh_dimension();

  //--------------------------------------------------------------------------------------------------------//
  //add variables
  opt.number_of_bands = calculate_number_of_bands( );

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

  system->add_matrix("H_real"); //add matrix for a real part of the Hamiltonian

  _H_real = & (system->get_matrix("H_real"));

  system->add_matrix("H_imag");//add matrix for an imaginary part of the Hamiltonian

  _H_imag = &(  system->get_matrix("H_imag") );

  system->add_matrix("S_real"); //add matrix for S matrix

  _S_real = &( system->get_matrix("S_real") );

  //system->add_matrix("S_imag"); //add matrix for S matrix

  //_S_imag = &( system->get_matrix("S_imag") );


  //peiodicity can not be changed between runs because that will require cleaning of the DOF constraint table
  solver_opt.periodicity[0]          = get_option("x-periodicity", false);
  solver_opt.periodicity[1]          = get_option("y-periodicity", false);
  solver_opt.periodicity[2]          = get_option("z-periodicity", false);
  
   
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
  //---------------------------------------------------------------------------------------------------------//


  Scaling& scaling = get_scaling();

  scaling.set_length_scaling(Constants::bohr_radius);

  scaling.set_calc_mesh_units(get_mesh_units());

  system->init();

  
  
  // We add a second system just to contain the density
  create_equation_system("linear");
  TiberLinearSystem& linsys = get_equation_system<TiberLinearSystem>(0);
  linsys.add_variable("qdens", libMeshEnums::FIRST);
  linsys.init();
  
  //------------------------------------------------------------------------------------------------------//
  //kp bands map
  {
    MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
    const Elem* elem = *el;
    
    EFAbulkHamiltonian* element_hamiltonian =
      get_bulk_model<EFAbulkModel>(elem)->get_Hamiltonian_model();
    
    opt.kp_bands = element_hamiltonian->get_kp_bands_map();
    
    opt.degeneracy = element_hamiltonian->get_degeneracy();
  }
  //------------------------------------------------------------------------------------------------------//
 

  parse_options();

  EigenvalueProblem::init_kspace();

}



//===========================================================//
void EnvelopFunctionApprox::do_solve()
{

 // check unused tags in solver_options
 const ModelOptions& sol_opt = get_solver_options();
 sol_opt.find_option("simulation"); // remove simulation name 
 sol_opt.check_unused();

 if (opt.job == BULKEIGENSTATES )
 {
   solve_bulk();
 }
 else
 { 


   if (_calculate_density && _k_vector[0] == 0.0 
                          && _k_vector[1] == 0.0 
                          && _k_vector[2] == 0.0 )

   {
     estimate_spectrum_shift();
     apply_bc();
     calculate_density_analytic();
   }
   else
   {
     Point k_vec;
     for (short i = 0; i < 3; i++) k_vec(i)=_k_vector[i];
     solve_for_kpoint(k_vec);
   }

 
   // we have to redeclare the solution variables to adjust the number
   // of eigenstates
   const unsigned int num_states = _solution.size();
   unsigned int dim = get_mesh().mesh_dimension();
   string units("1/cm");
   if (dim == 2)
     units = "1/cm^2";
   else if (dim == 3)
     units = "1/cm^3";
   declare_solution(ProbabilityDensity, NTUPLE, NODES, units, num_states);
   declare_solution(EigenEnergy, NTUPLE, GLOBAL, "eV", num_states);
   declare_solution(Occupation, NTUPLE, GLOBAL, "", num_states);
   declare_solution(EigenEnergyOnMesh, NTUPLE, NODES, "eV", num_states);
 }


}

//===========================================================//
void EnvelopFunctionApprox::do_solve_for_kpoint(const Point& k_point)
{
   

  if (opt.job == BULKEIGENSTATES )
  {
    solve_bulk();
  }
  else
  { 
    estimate_spectrum_shift();
    
    apply_bc(); 
    
    if (verbose() > 0)
    {
      ostringstream os;
      os << "(EFA) Solving for k = ( "<<
        _k_vector[0]<<" "<<_k_vector[1]<<" "<<_k_vector[2]<< " )";
      Messages::info(os.str());
    }
    
    //assemble();
    
    solve_eigen_value_problem(solver_opt.number_of_eigenstates, 
                              solver_opt.spectrum_shift/Hartree);

  }
}

//===========================================================//
void EnvelopFunctionApprox::calculate_Hamiltonian_and_S(void)
{


  _H_real->zero();
  _H_imag->zero();
  _S_real->zero();

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
  //QGauss qrule (dim, FIFTH);
  AutoPtr<QBase> qrule(QBase::build(_quadrature_type, dim, SECOND));

  // Tell the finite element object to use our quadrature rule.

  fe->attach_quadrature_rule (qrule.get());

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

  double initval = 1.0;
  if (_quadrature_type == QTRAP)
  {
    // is this correct, or should we use n_local_dofs()?
    _sqrt_S_inv.resize(dof_map.n_dofs());
    _sqrt_S_inv.assign(_sqrt_S_inv.size(), 0.0);
    //vector<double> test(_sqrt_S_inv);


    double scale = Constants::bohr_radius * 1e9;
    if (dim > 1)
      scale *= Constants::bohr_radius * 1e9;
    if (dim > 2)
      scale *= Constants::bohr_radius * 1e9;
    scale = 1.0 / scale;

    MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh->active_elements_end();
    for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;
      int n_nodes = elem->n_nodes();
      fe->reinit(elem);

      for (unsigned int b = 0; b < opt.number_of_bands; b++)
      {
        dof_map.dof_indices(elem, dof_indices_component, psivar[b]);
        const unsigned int n_dofs = dof_indices_component.size();

        for (unsigned int qp=0; qp < (*qrule).n_points(); qp++)
          for (unsigned int p = 0; p < n_dofs; p++)
            _sqrt_S_inv[dof_indices_component[p]] += JxW[qp] * phi[p][qp] * phi[p][qp];
      }
    }

    // build the square root and invert
    for (size_t i = 0; i < _sqrt_S_inv.size(); i++)
      _sqrt_S_inv[i] = 1.0 / sqrt(_sqrt_S_inv[i]);

  }


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  double electric_potential = 0;

  EFAbulkHamiltonian* element_hamiltonian;



  for ( ; el != end_el ; ++el)
    {//el
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      element_hamiltonian = get_bulk_model<EFAbulkModel>(elem)->get_Hamiltonian_model();

      element_hamiltonian->set_temperature(_temp_interface.get_temperature( elem, elem->centroid()));

      element_hamiltonian->set_k_vector(_k_vector);

      element_hamiltonian->calculate_Hamiltonian_k_par();


      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs   = dof_indices.size();


      ham_real.resize(n_dofs, n_dofs);
      ham_imag.resize(n_dofs, n_dofs);
      if (solver_opt.discretization_method == FEM)
        s_real.resize(n_dofs, n_dofs);

      fe->reinit (elem);

      for (unsigned int qp=0; qp < (*qrule).n_points(); qp++)
      {//qp
        //--------------------------------------------------------------------------------
        /*
	    We assume that strain and electric potential may be different for different quadrature points
	    It is done for a sake of a multiscale generalization
         */
        Tensor2Sym strain_crystal_system(0);
        _strain_interface.get_crystal_strain(elem, q_point[qp], strain_crystal_system);


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
                               * (-1.0); //Complex(0.0,-1.0) * Complex(0.0, -1.0);

                  }






                ham_real_sub(p1,p2) += value.real();
                ham_imag_sub(p1,p2) += value.imag();

              }
            }




            //S-matrix
            if (solver_opt.discretization_method == FEM)
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


      //if (solver_opt.discretization_method == FEM) 
      //  ham_real.add( - solver_opt.spectrum_shift/Hartree, s_real);//apply spectrum shift.


      if (_quadrature_type == QTRAP)
      {
        // apply S^-1/2 H S^-1/2
        for (unsigned int i = 0; i < n_dofs; i++)
        {
          for (unsigned int j = 0; j < n_dofs; j++)
          {
            double scale =
                _sqrt_S_inv[dof_indices[i]] * _sqrt_S_inv[dof_indices[j]];
            ham_real(i, j) *= scale;
            ham_imag(i, j) *= scale;
            // this is not needed, as we do not solve a generalized problem
            // in this case
            if (solver_opt.discretization_method == FEM)
              s_real(i, j) *= scale;
          }
        }
      }


      vector<unsigned int> dof_indices_tmp;

      if (solver_opt.discretization_method == FEM)
      {
	dof_indices_tmp = dof_indices;

	dof_map.constrain_element_matrix(s_real, dof_indices_tmp);
	_S_real->add_matrix(s_real,dof_indices_tmp);
      }

      dof_indices_tmp = dof_indices;

      dof_map.constrain_element_matrix(ham_real, dof_indices_tmp);
      _H_real->add_matrix(ham_real,dof_indices_tmp);

      dof_indices_tmp = dof_indices;

      dof_map.constrain_element_matrix(ham_imag, dof_indices_tmp);
      _H_imag->add_matrix(ham_imag,dof_indices_tmp);



    }


//this is only to test
  /*
  _H_real->print_matlab("ham_r_matlab.m");
  //_H_imag->print_matlab("ham_i_matlab.m");
  _S_real->print_matlab("s.m");
  */




  copy_H_to_solver( );



  if (solver_opt.discretization_method == FEM)  copy_S_to_solver( );





  //  dof_map.print_dof_constraints();

}


//============================================================//
void EnvelopFunctionApprox::estimate_spectrum_shift(void)
{
 if (opt.estimate_spectrum_shift)
 {
   solver_opt.spectrum_shift = get_band_edge(opt.particle);
   if (opt.particle == "el") solver_opt.spectrum_shift -= 0.05;
   if (opt.particle == "hl") solver_opt.spectrum_shift += 0.05;    

   std::cout<<"  (EFA) Estimated guess (eV): " << solver_opt.spectrum_shift << std::endl;
 }

}
//=============================================================//
double EnvelopFunctionApprox::get_new_spectrum_shift(void)
{
  double st_shift_value ;

  int v = verbose();
  verbose() = 0;
  read_SLEPC_solution(1);
  verbose() = v;

  assert(_solution.size() == 1);


  //st_shift_value = (_solution[0].eigen_energy - opt.spectrum_shift)/Hartree;
  st_shift_value = (_solution[0].eigen_energy)/Hartree;

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

  if (number_of_converged_solutions < number_of_ev)
    throw SolveFailedException(get_name() + " obtained less eigenvalues than requested.");



  //--------------------------------------------------------------------
  //read eigenvalues
  //store also eigenvalue index for sorting

  vector<EigenvalueProblem::eigen_state>  ev(number_of_converged_solutions);
  double shift = EigenSolver::get_shift() * Hartree;

  for (unsigned ind = 0; ind < number_of_converged_solutions; ind++)
  {
    ev[ind].energy =  EigenSolver::get_eigenvalue(ind) * Hartree; // + shift;
    ev[ind].index = ind;
  }

  //---------------------------------------------------------------------
  //sorting of the solutions

  if (opt.particle == "el") sort( ev.begin(), ev.end(), EigenvalueProblem::compare_eigen_energy_electrons);

  if (opt.particle == "hl") sort( ev.begin(), ev.end(), EigenvalueProblem::compare_eigen_energy_holes);

  if (verbose() > 1)
  {
    Messages m;
    ostringstream os;
    os << "converged eigenenergies (" << number_of_converged_solutions
        << "):";
    m.info(os.str());
    m.indent();

    os.str("");
    for (unsigned int i = 0; i < number_of_converged_solutions; ++i)
    {
      os << ev[i].energy << " ";
      if (i%8 == 7)
        os << "\n";
    }
    Messages::info(os.str());
    m.newline();
  }


  //----------------------------------------------------------------------
  //let us find the ground electron or hole level
  unsigned int ground_state_index = 0;
  bool finish = false;


  //cout<<"  Check solutions against shift (eV) "<<shift<< endl;

  for (unsigned int i = 0; (i < number_of_converged_solutions && (!finish) ); i++)
  {
    if (opt.particle == "el")
    {
      if (ev[i].energy > shift)
      {
	ground_state_index = i;
	finish = true;
      }
    }
    else
    {
      if (ev[i].energy < shift)
      {
	ground_state_index = i;
	finish = true;
      }
    }
  }


  if (!finish)
  {
    std::cout<<"Found "<< number_of_converged_solutions << " eigenvalues" << std::endl;
    for (unsigned int i = 0; (i < number_of_converged_solutions && (!finish) ); i++)
    {  
      std::cout<<ev[i].energy<<" ";
      if (i%8 == 7) std::cout << "\n";
    }
    throw ModelErrorException("EnvelopFunctionApprox: ground state is not found. Correct spectrum_shift or increase number of states");
  }

  unsigned int solution_size;
  if (number_of_converged_solutions - ground_state_index < number_of_ev)
    solution_size = number_of_converged_solutions - ground_state_index;
  else
    solution_size = number_of_ev;

  if (solution_size < number_of_ev)
  {
    ostringstream os;
    os << "kp found only " << solution_size << " eigenvalues instead of " << number_of_ev;
    Messages::warning(os.str());
  }


  //--------------------------------------------------------------------
  //read eigenvectors

  //read solutions - only independent dofs

  //----------------------------------------------------------------------
  {
    EnvelopFunctionApprox::eigen_problem_solution temp1;
    temp1.eigen_energy = 0;
    temp1.electro_chem_pot = 0;
    temp1.temperature = SimulationOptions::temperature;
    temp1.eigen_vector.resize(number_of_all_dofs, Complex(0.0, 0.0));

    _solution.clear();
    _solution.resize(solution_size, temp1);
  }


  // set up a map to recover the right eigen-pairs after sorting
  map<unsigned int, unsigned int>  global_to_sol_index;
  map<unsigned int, unsigned int>  :: iterator it;
  

  for (unsigned int i = ground_state_index; i < ground_state_index + solution_size ; i++)
  {
    global_to_sol_index.insert( make_pair( i - ground_state_index, ev[i].index )  );
    _solution[i - ground_state_index].eigen_energy = ev[i].energy;
  }
   //---------------------------------------------------------------------- 
  for (unsigned int ind = 0; ind < solution_size; ind++)
  {
    
    vector<Complex> temp;
    
    it = global_to_sol_index.find(ind);
    
    if (  it  !=  global_to_sol_index.end() )
    {
      unsigned int solution_number = it->second;
      
      EigenSolver::get_eigen_vector(solution_number, temp);




      //-----------------------------------------------------------------------------
      //put independent dofs in the eigenvectors that may contain also non independent dofs
      for (unsigned j = 0; j < number_of_all_dofs; j++)
      {
        if (new_dofs[j].independent)
        {
          
          _solution[ind].eigen_vector[j] = temp[new_dofs[j].new_number];
          
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
            
            _solution[ind].eigen_vector[j] += ( c->second ) * 
              _solution[ind].eigen_vector[(c->first)];
          }
          
        }
        
      }
      
    }
    
      //------------------------------------------------------------------------

  }

  // apply transformation if needed
  transform_eigenstates();
  

  //normalization
  for (unsigned int i = 0; i < solution_size; i++)
  {
    const double norm = eigenstate_norm(i);
    
    const unsigned int n1 =  _solution[i].eigen_vector.size();
    
    for (unsigned int j = 0; j < n1; j++)
      _solution[i].eigen_vector[j] /= Complex(norm, 0.0);

//    if (!check_confinement(_solution[i].eigen_vector))
//    {
//      ostringstream os;
//      os << "State " << i << " is not confined!";
//      Messages::warning(os.str());
//    }
  }
  
  


  //Fermi energy calculation
  
  if (poisson_equation != NULL)
    for (unsigned int i = 0; i < solution_size; i++)
    {
      _solution[i].particle = opt.particle;
      
      _solution[i].statistics = "Fermi";
      
      _solution[i].electro_chem_pot = calculate_fermi_averaged(i);
    }
  
  
  //Temperature calculation
  if (_temp_interface.has_simulation())
    for (unsigned int i = 0; i < solution_size; i++)
    {
      _solution[i].temperature = calculate_temperature_averaged(i);
      
    }
  

}




bool
EnvelopFunctionApprox::check_confinement(const vector<Complex>& state)
{
  bool confined = true;

  IDSet reg_ids;
  get_environment().get_device().extract_physical_regions(
      get_option("check_confinement", ""), reg_ids);

  if (reg_ids.empty()) return confined;



  DofMap& dof_map = system->get_dof_map();

  FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

  AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  QGauss qrule (dim, SECOND);

  fe -> attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  std::vector<unsigned int> dof_indices;


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  Complex sum(0.0, 0.0);

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    if (!reg_ids.count(elem->subdomain_id()))
      continue;

    fe->reinit (elem);

    for (short psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
    {
      dof_map.dof_indices (elem, dof_indices, psi_index);
      const unsigned int n_psi_dofs = dof_indices.size();

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {//qp
        for (unsigned int p1=0; p1<n_psi_dofs; p1++)
        {
          Complex eigen_f_value1(state[dof_indices[p1]]);
          for (unsigned int p2=0; p2<n_psi_dofs; p2++)
          {
            Complex eigen_f_value2(state[dof_indices[p2]]);
            Complex tmp = JxW[qp] * phi[p1][qp] * eigen_f_value1 *
                phi[p2][qp] * conj(eigen_f_value2);
            sum += tmp;
          }
        }
      }
    }
  }



  double confinement_threshold = get_option("confinement_threshold", 0.2);
  confined = (sqrt(abs(sum)) > confinement_threshold) ? true : false;

  return confined;
}


void
EnvelopFunctionApprox::plot_globaldata(void)
{

  string outdir = get_output_directory();

  string filename(outdir + "/" + get_output_filename() + ".dat");
  ofstream file;
  file.open(filename.c_str());
  if (file.good())
  {
    // header
    file << "# EFA eigenstates (" << get_name() << ")\n";
    file << "# Particle: ";
    if (opt.particle == "el")
      file << "electron";
    else
      file << "hole";
    file << "\n#\n";
    file << "# Index" << setw(12) << "EigenEnergy" << setw(15) << "Occupation"
        << setw(12) << "FermiLevel" << setw(12) << "Temperature" << "\n";

    for (unsigned int i = 0; i < _solution.size(); i++)
    {
      file << setw(7) << i << " "
          << setw(11) << _solution[i].eigen_energy << " "
          << setw(14) << Fermi_statistics_probability(_solution[i].eigen_energy,
              _solution[i].electro_chem_pot, _solution[i].temperature) << " "
          << setw(11) << _solution[i].electro_chem_pot << " "
          << setw(11) << _solution[i].temperature << "\n";
    }
  }

}


//=======================================================================//


EnvelopFunctionApprox:: ~EnvelopFunctionApprox(void)
{
  // es->delete_system(system_name);
}

//=======================================================================//



void EnvelopFunctionApprox::transform_eigenstates(void)
{
  // if we use QTRAP (diagonal overlap), we must transform the eigenstate
  // with S^-1/2 as we transformed the system to have unit overlap
  if (_quadrature_type != QTRAP) return;

  size_t num_states = _solution.size();

  for (size_t s = 0; s < num_states; s++)
  {
    vector<Complex>& eigvec =  _solution[s].eigen_vector;

    size_t n_dofs = eigvec.size();
    if (n_dofs != _sqrt_S_inv.size())
      cerr << "eigvec : " << eigvec.size() << " " << "S : " << _sqrt_S_inv.size() << endl;

    for (size_t i = 0; i < n_dofs; i++)
      eigvec[i] *= _sqrt_S_inv[i];

  }

}

//-----------------------------------------------------------------------------//
double  EnvelopFunctionApprox::eigenstate_norm(unsigned int state_number)
{
  double  result;

  const vector< Complex > &  eigen_vector =  _solution[state_number].eigen_vector;



  DofMap& dof_map = system->get_dof_map();

  FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

  // AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  //AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type)  );
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
		      temp += ( JxW[qp] * phi[p1][qp] * eigen_f_value1 * 
                                      phi[p2][qp] * conj(eigen_f_value2) );
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





  const vector< Complex >&   eigen_vector =  _solution[i].eigen_vector;




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
		      result += ( JxW[qp] * phi[p1][qp] * eigen_f_value1 *  
                                phi[p2][qp] * conj(eigen_f_value2) ) * chem_pot_value_eV;

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





  const vector< Complex >&   eigen_vector =  _solution[i].eigen_vector;




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

  for ( ; el != end_el ; ++el)
    {//el

      const Elem* elem = *el;
      fe->reinit (elem);


      Point center = elem->centroid();
      Temperature = _temp_interface.get_temperature(elem, center);




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
		      result += ( JxW[qp] * phi[p1][qp] * eigen_f_value1 *  
                                  phi[p2][qp] * conj(eigen_f_value2) ) * Temperature;

		    }
		}

	    }




	}


    }








  return(result.real());



}

//--------------------------------------------------------------------------//
void EnvelopFunctionApprox::do_assemble(const ModelOptions& options)
{
  calculate_Hamiltonian_and_S();
}




//===========================================================//





//=================================================================//

  
void EnvelopFunctionApprox::calculate_density_analytic(void)
{
  unsigned int dim = get_mesh().mesh_dimension();

  unsigned int num_states = solver_opt.number_of_eigenstates;


  Point kvector_0;
  Point kvector_1;
  Point kvector_2;
  if (dim < 3)
    kvector_1(2) = opt.k_val;
  if (dim < 2)
    kvector_2(1) = opt.k_val;


  // for now, this does only analytic integration

  Messages m;
  m.info("(EFA) Calculating quantum density");
  m.indent();

  double a_B =  Constants::bohr_radius;
  double scaling =  1.0 / (a_B * a_B * a_B * 1.0e6 );

  vector<double> energy_k_0;
  vector<double> energy_k_1;
  vector<double> energy_k_2;
  vector<double> energy_k_3;
  vector<double> effective_mass(num_states);

  bool solve_twice = solver_opt.solve_ev_problem_twice;
  double spectrum_shift = solver_opt.spectrum_shift/Constants::Hartree;

  if (solver_opt.solve_ev_problem_twice)
  {
    solver_opt.solve_ev_problem_twice = false;

    set_k_point(kvector_0);

    if (verbose() > 0)
      Messages::info("Solve to obtain spectrum shift ... ", false);

    int v = verbose();
    verbose() = 0;
    solve_eigen_value_problem(1, spectrum_shift);
    verbose() = v;

    get_eigenenergies(energy_k_0);

    spectrum_shift = get_new_spectrum_shift();
  }


  // [0 0 1]
  if (dim < 3)
  {
    ostringstream os;
    os << "Solving for k1 = ( ";
    kvector_1.write_unformatted(os, false);
    os << ")";
    if (verbose() > 0) m.info(os.str());

    set_k_point(kvector_1);
    solve_eigen_value_problem(num_states, spectrum_shift);
    get_eigenenergies(energy_k_1);
  }

  // [0 1 0]
  if (dim < 2)
  {
    ostringstream os;
    os << "Solving for k2 = ( ";
    kvector_2.write_unformatted(os, false);
    os << ")";
    if (verbose() > 0) m.info(os.str());

    set_k_point(kvector_2);
    solve_eigen_value_problem(num_states, spectrum_shift);
    get_eigenenergies(energy_k_2);

    // [0 1 1]
    if (!opt.assume_paraboloid)
    {
      kvector_2(2) = opt.k_val;

      ostringstream os;
      os << "Solving for k3 = ( ";
      kvector_2.write_unformatted(os, false);
      os << ")";
      if (verbose() > 0) m.info(os.str());

      set_k_point(kvector_2);
      solve_eigen_value_problem(num_states, spectrum_shift);
      get_eigenenergies(energy_k_3);
    }
  }

  if (dim < 3)
  {
    ostringstream os;
    os << "Solving for k0 = ( ";
    kvector_0.write_unformatted(os, false);
    os << ")";
    if (verbose() > 0) m.info(os.str());
  }

  set_k_point(kvector_0);
  // we solve with +1 state to be used in dd to define the classical boundary
  solve_eigen_value_problem(num_states + 1, spectrum_shift);
  get_eigenenergies(energy_k_0);

  solver_opt.solve_ev_problem_twice = solve_twice;


  double Eh_k2 = Constants::Hartree * opt.k_val * opt.k_val;

  if (dim == 1)
  {
    for (unsigned int i = opt.first_state; i < num_states; i++)
    {
      double imass11 = 2.0 * (energy_k_0[i] - energy_k_1[i]) / Eh_k2;
      double imass22 = 2.0 * (energy_k_0[i] - energy_k_2[i]) / Eh_k2;
      double imass12 = 0;

      if (!opt.assume_paraboloid)
      {
        imass12 = (energy_k_0[i] - energy_k_3[i]) / Eh_k2
            - 0.5 * (imass11 + imass22);
      }
      double det = abs(imass11 * imass22 - imass12 * imass12);

      effective_mass[i] = 1.0 / sqrt(det);
    }
  }
  else if (dim == 2)
  {
    for (unsigned int i = opt.first_state; i < num_states; i++)
    {
      double imass = 2.0 * abs(energy_k_0[i] - energy_k_1[i]) / Eh_k2;
      effective_mass[i] = 1.0 / imass;
    }
  }
  if ((verbose() > 1) && (dim < 3))
  {
    m.info("effective masses:");
    m.indent();

    stringstream os;
    for (unsigned int i = 0; i < num_states; ++i)
    {
      os << effective_mass[i] << " ";
      if (i%8 == 7)
        os << "\n";
    }
    Messages::info(os.str());
    m.newline();
  }


  FEType fe_type = system->variable_type(0);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

  DofMap& dof_map = system->get_dof_map();
  std::vector<unsigned int> dof_indices;

  // The qdens_sys system contains the nodal quantum density
  TiberLinearSystem& qdens_sys = get_equation_system<TiberLinearSystem>();
  DofMap& dof_map_qdens = qdens_sys.get_dof_map();
  std::vector<unsigned int> dof_indices_qdens;
  NumericVector<Number>& qdens = *qdens_sys.solution;
  qdens.zero();

  // we need the connectivity of the nodes to not double count
  vector<int> connectivity(qdens.size(), 0.0);
  {
    MeshBase::const_element_iterator el = get_mesh().active_elements_begin();
    const MeshBase::const_element_iterator end_el = get_mesh().active_elements_end();
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map_qdens.dof_indices(elem, dof_indices_qdens, 0);
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
        connectivity[dof_indices_qdens[n]]++;
    }

  }


  for (unsigned int i = opt.first_state; i < num_states; i++)
  {
    double fermi_energy = _solution[i].electro_chem_pot;
    double kT = _solution[i].temperature * Constants::k_Boltzmann;

    double energy = _solution[i].eigen_energy;
    double mass_factor = 1.0;
    double dos_factor = 1.0;

    if (dim == 1)
    {
      // 1D is correct and tested
      mass_factor = effective_mass[i] * kT / (2.0 * M_PI * Constants::Hartree);

      double exp_arg = (fermi_energy - energy) / kT;
      if (opt.particle == "hl")
        exp_arg = -exp_arg;

      dos_factor = (exp_arg < -20) ? exp(exp_arg) : log(1.0 + exp(exp_arg));
    }
    else if (dim == 2)
    {
      // 2D is correct and tested
      mass_factor = sqrt(kT *  effective_mass[i] / (2.0 * M_PI * Constants::Hartree));
      double exp_arg = (fermi_energy - energy) / kT;
      if (opt.particle == "hl")
        exp_arg = -exp_arg;

      dos_factor = TiberMath::fermidirac_mhalf(exp_arg);
    }
    else if (dim == 3)
      dos_factor = Fermi_statistics_probability(energy, fermi_energy,
          _solution[i].temperature);

    dos_factor *= mass_factor * scaling * opt.degeneracy;


    MeshBase::const_element_iterator el = get_mesh().active_elements_begin();
    const MeshBase::const_element_iterator end_el = get_mesh().active_elements_end();
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map_qdens.dof_indices(elem, dof_indices_qdens, 0);

      for (int psi_index = 0; psi_index < opt.number_of_bands; psi_index++)
      {
        dof_map.dof_indices(elem, dof_indices, psi_index);

        for (unsigned int n = 0; n < elem->n_nodes(); n++)
        {
          double psi = abs(_solution[i].eigen_vector[dof_indices[n]]);
          double val = dos_factor * psi * psi / connectivity[dof_indices_qdens[n]];
          qdens.add(dof_indices_qdens[n], val);
        }

      }
    }
  }
  qdens.close();
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
  element_hamiltonian = get_bulk_model<EFAbulkModel>(elem)->get_Hamiltonian_model();

  std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
	    model_Ham = ( element_hamiltonian->get_Hamiltonian() );

  result = model_Ham.size();

  if (result == 0)  throw InitFailedException("EnvelopFunctionApprox: Hamiltonian is empty");

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    element_hamiltonian = get_bulk_model<EFAbulkModel>(elem)->get_Hamiltonian_model();


    const std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
      model_Ham = ( element_hamiltonian->get_Hamiltonian() );


    short result1 = model_Ham.size();

    if (result1 != result)
    {
      std::string mess = "EnvelopFunctionApprox: Hamiltonians of different";
      mess += " materials have different number of bands";
      throw InitFailedException(mess);
    }

  }

  return(result);

}


//=======================================================================================//

void EnvelopFunctionApprox::solve_bulk(void)
{
  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();


  bool found = false;
  const Elem* mat_elem;
  std::cout<<"bulk calculation"<<std::endl;
  for ( ; (el != end_el) && (!found) ; ++el)
  {
    const Elem* elem = *el;
    if (elem->contains_point(_bulk_point))
    {
      found = true;
      mat_elem = elem;
      _bulk_mat_element = mat_elem;
      break;
    }

  }

  if (!found) throw SolveFailedException("Bad material point\n");

  Point qp = mat_elem->centroid();

  EFAbulkHamiltonian* element_hamiltonian;

  std::cout<<"elem H"<<std::endl;

  element_hamiltonian = get_bulk_model<EFAbulkModel>(mat_elem)->get_Hamiltonian_model();

  element_hamiltonian->set_k_vector(_k_vector);

  element_hamiltonian->calculate_Hamiltonian_k_par();

  std::cout<<"strain"<<std::endl;
  Tensor2Sym strain_crystal_system(0);
  _strain_interface.get_crystal_strain(mat_elem, qp, strain_crystal_system);

  std::cout<<"(EP) strain exx "<<strain_crystal_system(1,1)<<std::endl;
  std::cout<<"(EP) strain eyy "<<strain_crystal_system(2,2)<<std::endl;
  std::cout<<"(EP) strain ezz "<<strain_crystal_system(3,3)<<std::endl;

  std::cout<<"pot"<<std::endl;

  double electric_potential = 0;
  if (opt.consider_potential_bulk)
      electric_potential = get_electric_potential( mat_elem, qp );


  element_hamiltonian->apply_strain_and_potential(strain_crystal_system, electric_potential);


  std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
	    model_Ham = ( element_hamiltonian->get_Hamiltonian() );

  std::cout<<"number of bands: " << opt.number_of_bands<<std::endl;

  std::complex<double> ham_matrix[opt.number_of_bands * opt.number_of_bands ];

  for (unsigned int band1 = 0; band1 < opt.number_of_bands; band1++)
    for (unsigned int band2 = 0; band2 < opt.number_of_bands; band2++)
    {
      ham_matrix[band1 + band2 * opt.number_of_bands] = model_Ham[band1][band2].constant;

    }


  _solution.clear();
  _solution.resize(opt.number_of_bands);

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

       _solution[i].eigen_energy = eigvals[i]*Hartree;
       _solution[i].eigen_vector.resize(opt.number_of_bands);
       for (short j = 0; j < opt.number_of_bands ; j++)
       {
	 _solution[i].eigen_vector[j] = ham_matrix[i * opt.number_of_bands + j];
       }
     }
  }






  unsigned int n = _solution.size();

  for (unsigned int i = 0; i < n ; i++)
  {
    _solution[i].electro_chem_pot = 0;

    if (poisson_equation != NULL)
    	_solution[i].electro_chem_pot = get_electro_chem_potential(mat_elem);

    _solution[i].temperature =
        _temp_interface.get_temperature(mat_elem, mat_elem->centroid());
  }


}
