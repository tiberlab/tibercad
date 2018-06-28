// $Id$

#include "SimulationEnvironment.h"
#include "EnvelopFunctionApprox.h"
#include "ModelOptions.h"
#include "EFAbulkModel.h"
#include "Material.h"
#include "Boundary.h"
#include "TiberMath.h"
#include "TiberLinearSystem.h"
#include "SimulationOptions.h"
#include "tensor.h"

#include "EigenSolver.h"
#include "quadrature_gauss.h"


#include <edge_edge2.h>
#include <equation_systems.h>
#include <dense_submatrix.h>
#include <quadrature_gauss.h>


#include "Messages.h"


extern "C"
{
   void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& LDA,
	       double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info);
}

using namespace std;
using namespace Constants;
using namespace libMesh;


namespace
{
  bool compare_eigen_energy(const EigenvalueProblem::eigen_state& state1,
      const EigenvalueProblem::eigen_state& state2)
  {
    if ((state1.particle == "hl") && (state1.particle == state2.particle))
      return (state2.energy < state1.energy);

    return(state1.energy< state2.energy);
  }
}

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
  /*
  if (poisson_equation==NULL) 
  {
    Messages::warning("trying to find band-edge without drift-diffusion: return 0");
    return 0.0;
  }
  */
  if (((particle == "el") && (_cb_edge.second == INVALID_ID)) ||
      ((particle == "hl") && (_vb_edge.second == INVALID_ID)))
    return 0.0;

  // 2014-09-05 this gives sometimes wrong results, try with centroid
  //vector<double> values(elem->n_nodes());
  //vector<Point> p(elem->n_nodes());

  //for (size_t i = 0; i < elem->n_nodes(); ++i)
  //      p[i] = elem->point(i);
  
  vector<double> values(1);
  vector<Point> p(1, elem->centroid());

  bool electron = (particle == "el") || (particle == "Ec");

  if (electron)
    _cb_edge.first->get_solution(elem, _cb_edge.second, values, p);
  else
    _vb_edge.first->get_solution(elem, _vb_edge.second, values, p);

  double bedge = values[0];

  for (size_t i = 1; i < values.size(); ++i)
  {
     double temp = values[i];
     if (electron)
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

  _el_pot.first->get_solution(elem, _el_pot.second, values, qp);

  return values[0];
}
//---------------------------------------------------------------------------------//

inline double EnvelopFunctionApprox::get_el_electro_chem_potential(const Elem* elem) const
{
  vector<double> values;
  vector<Point> qp(1, elem->centroid());

  _el_elchem.first->get_solution(elem, _el_elchem.second, values, qp);

  return values[0];
}


inline double EnvelopFunctionApprox::get_hl_electro_chem_potential(const Elem* elem) const
{
  vector<double> values;
  vector<Point> qp(1, elem->centroid());

  _hl_elchem.first->get_solution(elem, _hl_elchem.second, values, qp);

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
    values[i] = Fermi_statistics_probability(_solution[i].eigen_energy,
                                             _solution[i].electro_chem_pot,
                                             _solution[i].temperature,
                                             _solution[i].particle);
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
  declare_solution(EnvelopeFunctions, NTUPLE, NODES, "");
  add_alias("EigenFunctions", ProbabilityDensity);
  declare_solution(EigenEnergy, NTUPLE, GLOBAL, "eV");
  declare_solution(Occupation, NTUPLE, GLOBAL, "");
  declare_solution(EigenEnergyOnMesh, NTUPLE, NODES, "eV");

  if (plot_solution(EigenEnergy) && plot_solution(ProbabilityDensity))
    add_plot_variable(EigenEnergyOnMesh);

  if (plot_solution("eDensity") ||
      plot_solution("hDensity") ||
      plot_solution("QuantumDensity") ||
      get_options().has_submodel("QuantumDensity"))
    _calculate_density = true;
  if (_calculate_density)
  {
    declare_solution(eDensity, REAL, NODES, "1/cm^3");
    declare_solution(hDensity, REAL, NODES, "1/cm^3");
  }
  if (plot_solution("QuantumDensity"))
  {
    add_plot_variable(eDensity);
    add_plot_variable(hDensity);
  }
}


void
EnvelopFunctionApprox::do_set_to_remembered_solution(ID id)
{
  FEMEigenvalueProblem::do_set_to_remembered_solution(id);
  redeclare_solutions();
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
                                       _solution[sn].temperature,
                                       _solution[sn].particle);
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
    UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));
    const vector<vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &points);

    DofMap& dof_map = system->get_dof_map();
    std::vector<unsigned int> dof_indices;

    // number of states
    const unsigned int num_states = _solution.size();
    //values[ProbabilityDensity] = vector<double>(np * num_states, 0.0);
    values[ProbabilityDensity].clear();
    values[ProbabilityDensity].resize(np * num_states, 0.0);

    // they should all be the same size
    unsigned int n_dofs = phi.size();

    for (int psi_index = 0; psi_index < number_of_bands; psi_index++)
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

  if (values.count(EnvelopeFunctions))
  {
    unsigned int dim = get_mesh().mesh_dimension();
    double scale = Constants::bohr_radius * 1e2;
    if (dim > 1)
      scale *= Constants::bohr_radius * 1e2;
    if (dim > 2)
      scale *= Constants::bohr_radius * 1e2;

    FEType fe_type = system->variable_type(0);
    UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));
    const vector<vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &points);

    DofMap& dof_map = system->get_dof_map();
    std::vector<unsigned int> dof_indices;

    // number of states
    const unsigned int num_states = _solution.size();
    //values[ProbabilityDensity] = vector<double>(np * num_states, 0.0);
    values[EnvelopeFunctions].clear();
    values[EnvelopeFunctions].resize(np * num_states * 2 * number_of_bands, 0.0);

    // they should all be the same size
    unsigned int n_dofs = phi.size();

    for (int psi_index = 0; psi_index < number_of_bands; psi_index++)
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

          // the number of components should be set already to the right number
          unsigned int index = 2 * (number_of_bands * num_states * n + number_of_bands * sn + psi_index);
          values[EnvelopeFunctions][index] = real(value) / scale;
          values[EnvelopeFunctions][index + 1] = imag(value) / scale;
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

  bool do_edens = values.count(eDensity);
  bool do_hdens = values.count(hDensity);

  if (do_edens || do_hdens)
  {
    TiberLinearSystem& qdens_sys = get_equation_system<TiberLinearSystem>(0);
    libMesh::NumericVector<Number>& qdens = *qdens_sys.current_local_solution;

    FEType fe_type = qdens_sys.variable_type(0);
    UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));
    const vector<vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &points);

    DofMap& dof_map = qdens_sys.get_dof_map();
    std::vector<unsigned int> dof_indices_el, dof_indices_hl;
    dof_map.dof_indices(elem, dof_indices_el, 0);
    dof_map.dof_indices(elem, dof_indices_hl, 1);
    unsigned int n_dofs = phi.size();

    for (unsigned int n = 0; n < np; n++)
    {
      double value_e = 0;
      double value_h = 0;

      // the phi^2 factor comes from the fact that the more correct interpolation is the square
      // of the basis function, because the probability densities are the square of the states
      // NOTE: maybe one should check if this gives really a better result
      // NOTE: 2011-12-01 the above turned out to be wrong: phi is used only as linear interp. !
      for (unsigned int i = 0; i < n_dofs; i++)
      {
         value_e += phi[i][n] * qdens(dof_indices_el[i]);
         value_h += phi[i][n] * qdens(dof_indices_hl[i]);
      }

      if (do_edens)
        values[eDensity][n] = value_e;
      if (do_hdens)
        values[hDensity][n] = value_h;
    }

  }
}





//====================================================//
double EnvelopFunctionApprox::get_band_edge(const std::string& particle) 
{
  double band_edge = 0;

  if (_job == BULKEIGENSTATES)
  {
    if ((particle == "el") && (_cb_edge.second != INVALID_ID))
      _cb_edge.first->get_solution(_bulk_mat_element,
          _cb_edge.second, band_edge, _bulk_point);
    else if (_vb_edge.second != INVALID_ID)
      _vb_edge.first->get_solution(_bulk_mat_element,
          _vb_edge.second, band_edge, _bulk_point);

    // TODO when the bulk point is not on our processor???
  }
  else
  {
    MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

    bool condband = (particle == "el") || (particle == "Ec");

    band_edge = (condband ? 1.0 : -1.0 ) * numeric_limits<double>::max();


    for (; el != end_el ; ++el )
    {
      const Elem* elem = *el;

      double temp = get_band_edge(elem, particle);


      if (condband)
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

    // get the extremum over all MPI processes
    if (condband)
      this->get_solver_communicator().min(band_edge);
    else
      this->get_solver_communicator().max(band_edge);
  }

  return (band_edge);

}

//===================================================//
EnvelopFunctionApprox::EnvelopFunctionApprox(const ModelOptions& options)
 : FEMEigenvalueProblem(options),
   _calculate_density(false),
   poisson_equation(NULL),
   _bulk_mat_element(NULL)
{

  has_solution_vector(false);

}





//===================================================//
void EnvelopFunctionApprox::parse_options()
{

  FEMEigenvalueProblem::parse_options();

  const ModelOptions& solopts = get_solver_options();

  // this is useful afterwards
  char singleband = 0x0;
  if (get_options().has_submodel("Physics"))
  {
    ModelOptions::const_submodel_iterator it(get_options().submodels_begin("Physics"));
    const ModelOptions& opts = it->second;

    string model = opts.get_option("model", "8x8");

    if (model == "conduction_band")
      singleband = 'c';
    else if (model == "single_band")
    {
      string particle = opts.get_option("particle", "el");
      if (particle == "el")
        singleband = 'c';
      else if (particle == "hl")
        singleband = 'v';
      else
      {
        throw InitFailedException("In " +  get_name() + ": \'" + particle +
          " is unknown particle for 'single_band' model.");
      }
    }
    else if (model == "valence_band")
      singleband = 'v';
    else if (model == "2x2")
      singleband = 'b';
    else if (model == "6x6")
      singleband = 'v';

    // a quirky way to adjust degeneracy for certain singleband models
    switch (singleband)
    {
      case 'v':
        if (model == "6x6")
          break;

      case 'b':
      case 'c':
      case '1':
        opt.degeneracy *= 2;

      default:
        break;
    }

  }

  opt.num_hl_states = opt.num_el_states = 0;

  if (singleband != 'c')
  {
    opt.num_hl_states = get_option("number_of_eigenstates", 0);
    opt.num_hl_states = solopts.get_option("number_of_eigenstates", opt.num_hl_states);
    opt.num_hl_states = get_option("num_valence_eigenvalues", opt.num_hl_states);
    opt.num_hl_states = solopts.get_option("num_valence_eigenvalues", opt.num_hl_states);
    opt.num_hl_states = get_option("num_hole_states", opt.num_hl_states);
    opt.num_hl_states = solopts.get_option("num_hole_states", opt.num_hl_states);
  }

  if (singleband != 'v')
  {
    opt.num_el_states = get_option("number_of_eigenstates", 0);
    opt.num_el_states = solopts.get_option("number_of_eigenstates", opt.num_el_states);
    opt.num_el_states = get_option("num_conduction_eigenvalues", opt.num_el_states);
    opt.num_el_states = solopts.get_option("num_conduction_eigenvalues", opt.num_el_states);
    opt.num_el_states = get_option("num_electron_states", opt.num_el_states);
    opt.num_el_states = solopts.get_option("num_electron_states", opt.num_el_states);
  }

  solver_opt.number_of_eigenstates = opt.num_el_states + opt.num_hl_states;

  if (solver_opt.number_of_eigenstates == 0)
  {
    throw InitFailedException("In " +  get_name() + ": you must give at least "
        " one of 'num_electron_states' or 'num_hole_states'");
  }

  if (singleband && (singleband != 'b'))
  {
    if ((opt.num_hl_states > 0) && (opt.num_el_states > 0))
      throw InitFailedException("In " +  get_name() + ": you cannot ask for both hole and "
          " electron states when using 'single_band' model.");
  }

  //possible user override
  opt.degeneracy = get_option("degeneracy", opt.degeneracy);

  // for degeneracy = 1 we assure that number of states is even,
  // so we take both spin states
  if ((opt.degeneracy == 1) && (opt.num_el_states % 2 == 1))
  {
    opt.num_el_states += 1;
    Messages::warning("Number of electron eigenstates increased by 1 because of spin pairing");
  }
  if ((opt.degeneracy == 1) && (opt.num_hl_states % 2 == 1))
  {
    opt.num_hl_states += 1;
    Messages::warning("Number of hole eigenstates increased by 1 because of spin pairing");
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

  std::string  job_name = get_option("job", "eigenstates");
  if (job_name == "eigenstates")
    _job = EIGENSTATES;
  else if (job_name == "bulk")
    _job = BULKEIGENSTATES;
  else
    throw InitFailedException( "EnvelopeFunctionApprox: Incorrect job " + job_name );

  
  if (_job == BULKEIGENSTATES)
  {
    if (has_option("bulk_point"))
    {
      vector<double> point;

      get_option("bulk_point", point);
      point.resize(3, 0);
      for (short i = 0; i < 3; i++)  _bulk_point(i) = point[i];
    }
    else
    {
      throw InitFailedException( "You have to specify a bulk_point");
    }

    // get the bulk point's element

    // TODO to be adjusted for parallel

    MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

    bool found = false;

    for ( ; (el != end_el) && (!found) ; ++el)
    {
      const Elem* elem = *el;
      if (elem->contains_point(_bulk_point))
      {
        found = true;
        _bulk_mat_element = elem;
        break;
      }

    }

    if (!found) throw SolveFailedException("Bad bulk material point\n");
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
  opt.consider_strain_bulk = true;

  std::string  poisson_model_name = get_option("poisson_model_name","");
  poisson_model_name = get_option("poisson_simulation", poisson_model_name);

  if ( poisson_model_name != "" )
  {
    opt.consider_potential = true;

    opt.consider_potential_bulk = get_option("potential_in_bulk",true);

    opt.consider_strain_bulk = get_option("strain_in_bulk",true);

    poisson_equation  = find_simulation ( poisson_model_name );

    if (poisson_equation == NULL)
      throw InitFailedException( "Unknown poisson model " + poisson_model_name);

    potential_ID = poisson_equation->get_solution_id("ElPotential");
    _el_pot = make_pair(poisson_equation, potential_ID);


    if (potential_ID ==  INVALID_ID)
      throw InitFailedException( "Unknown variable ");


    el_electro_chem_pot_ID = poisson_equation->get_solution_id("eQFermi");
    _el_elchem = make_pair(poisson_equation, el_electro_chem_pot_ID);

    cb_band_edge_ID = poisson_equation->get_solution_id("Ec");
    _cb_edge = make_pair(poisson_equation, cb_band_edge_ID);

    hl_electro_chem_pot_ID = poisson_equation->get_solution_id("hQFermi");
    _hl_elchem = make_pair(poisson_equation, hl_electro_chem_pot_ID);

    vb_band_edge_ID = poisson_equation->get_solution_id("Ev");
    _vb_edge = make_pair(poisson_equation, vb_band_edge_ID);


  }
  else
  {

  _el_pot = find_solution_provider(
      get_option("electrostatic_potential", ""), "ElPotential");

  if (_el_pot.second != INVALID_ID)
  {
    opt.consider_potential = true;
    opt.consider_potential_bulk = get_option("potential_in_bulk",true);
    opt.consider_strain_bulk = get_option("strain_in_bulk",true);
  }
  else if (_el_pot.first != NULL)
  {
    throw InitFailedException( "Could not find electrostatic potential source "
        + get_option("electrostatic_potential", ""));
  }

  _el_elchem = find_solution_provider(
      get_option("el_electrochemical_potential", ""), "eQFermi");

  if ((_el_elchem.second == INVALID_ID) && (_el_elchem.first != NULL))
  {
    throw InitFailedException( "Could not find electrochemical potential source "
        + get_option("el_electrochemical_potential", ""));
  }


  _hl_elchem = find_solution_provider(
      get_option("hl_electrochemical_potential", ""), "hQFermi");

  if ((_hl_elchem.second == INVALID_ID) && (_hl_elchem.first != NULL))
  {
    throw InitFailedException( "Could not find electrochemical potential source "
        + get_option("hl_electrochemical_potential", ""));
  }


  _cb_edge = find_solution_provider(
      get_option("cb_edge", ""), "Ec");

  if ((_cb_edge.second == INVALID_ID) && (_cb_edge.first != NULL))
  {
    throw InitFailedException( "Could not find conduction band edge source "
        + get_option("cb_edge", ""));
  }

  _vb_edge = find_solution_provider(
      get_option("vb_edge", ""), "Ev");

  if ((_vb_edge.second == INVALID_ID) && (_vb_edge.first != NULL))
  {
    throw InitFailedException( "Could not find valence band edge source "
        + get_option("vb_edge", ""));
  }
  }

  //---------------------------------------------------------------------------------//
  //Heat model
  std::string heat_model_name = get_option("heat_model","");
  heat_model_name = get_option("temperature_simulation", heat_model_name);

  _temp_interface.set_simulation(heat_model_name);


  //--------------------------------------------------------------------------------------------//
  //Spectrum Shift
  //as default, we  estimate spectrum shift only if electric potential is defined
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

  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
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
  opt.k_val = 0.188; //0.01;
  opt.assume_paraboloid = false;

  opt.analytic_k_int = true;
  if (get_options().has_submodel("k_integration"))
    opt.analytic_k_int = false;

  opt.analytic_k_int = get_option("analytic_k_integration", opt.analytic_k_int);

  if (get_options().has_submodel("QuantumDensity"))
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("QuantumDensity"));
    ModelOptions& opts = it->second;

    opt.first_state = opts.get_option("first_state", opt.first_state);
    opt.k_val = opts.get_option("k_value", opt.k_val);
    opt.assume_paraboloid = opts.get_option("assume_diagonal_mass_matrix", opt.assume_paraboloid);
    opt.analytic_k_int = opts.get_option("analytic_k_integration", opt.analytic_k_int);
  }


  

}



void EnvelopFunctionApprox::do_print_info(void)
{
   Messages m;
   ostringstream os;

   os << "Number of bands: " << get_number_of_bands();
   m.info(os.str());

   os.str("");
   os << "Degeneracy: " << opt.degeneracy;
   m.info(os.str());
   m.newline();

}


//===================================================//
void EnvelopFunctionApprox::do_init( )
{
  FEMEigenvalueProblem::do_init();

  es = &(get_equation_systems());
  
  system_name = get_equation_system_name ( );
  
  es->add_system<LinearImplicitSystem> (system_name);
  
  system = &( es->get_system<LinearImplicitSystem>( system_name ) );

  dim = get_mesh().mesh_dimension();

//  set<ID> region_ids;
//  this->get_region_ids(region_ids);
//  set<unsigned short> reg_ids_s;
//  for (auto&& it : region_ids)
//  {
//    reg_ids_s.insert(static_cast<unsigned short>(it));
//  }

  //--------------------------------------------------------------------------------------------------------//
  //add variables
  number_of_bands = calculate_number_of_bands( );

  psi_name.clear();
  for (short i = 0; i < number_of_bands; i++)
    {
      std::ostringstream var_str;
      var_str << "psi" << i ;
      string name = var_str.str();
      psi_name.push_back(name);

      system->add_variable(name, FIRST, &get_region_ids());
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
  
   
  pair<Point, Point> bbox(get_environment().get_bounding_box());
  for (unsigned i = 0; i < 3; i++)
  {
    min_coord[i] = bbox.first(i);
    max_coord[i] = bbox.second(i);
  }


  Scaling& scaling = get_scaling();

  scaling.set_length_scaling(Constants::bohr_radius);

  //scaling.set_calc_mesh_units(get_mesh_units());
  system->init();

  
  // We add a second system just to contain the density
  create_equation_system("linear");
  TiberLinearSystem& linsys = get_equation_system<TiberLinearSystem>(0);
  linsys.add_variable("edens", libMeshEnums::FIRST,
                               libMeshEnums::LAGRANGE,
                               &(this->get_region_ids()));
  linsys.add_variable("hdens", libMeshEnums::FIRST,
                               libMeshEnums::LAGRANGE,
                               &(this->get_region_ids()));
  linsys.init();
  

  //------------------------------------------------------------------------------------------------------//
  //kp bands map
  {
    MeshBase::const_element_iterator  el     = this->active_local_elements_begin();
    MeshBase::const_element_iterator  end_el = this->active_local_elements_end();

    if (el != end_el)
    {
      const Elem* elem = *el;

      EFAbulkHamiltonian* element_hamiltonian =
          get_bulk_model<EFAbulkModel>(elem)->get_Hamiltonian_model();

      band_map = element_hamiltonian->get_kp_bands_map();

      // TODO there is some mess with the degeneracy, I believe
      opt.degeneracy = element_hamiltonian->get_degeneracy();
    }

    // if we don't have active pieces, band_map is size zero, so let's distribute
    // the band map from any other process.
    vector<map<short,short>::size_type> band_map_sizes(this->get_solver_communicator().size());
    this->get_solver_communicator().allgather(band_map.size(), band_map_sizes);

    int proc_id = 0;
    while ((proc_id < band_map_sizes.size()) && (band_map_sizes[proc_id] == 0))
      proc_id++;

    if (proc_id >= band_map_sizes.size())
      throw InitFailedException("efaschroedinger: cannot identify bulk bands");

    // now we can broadcast from processor
    this->get_solver_communicator().broadcast(band_map, proc_id);
    this->get_solver_communicator().broadcast(opt.degeneracy, proc_id);
  }
  //------------------------------------------------------------------------------------------------------//
 

  parse_options();

  init_kspace(ModelOptions());

  // Initialize identity permutation: does not have constrained dofs 
  //if (solver_opt.Dirichlet_bc_everywhere)
  EigenvalueProblem::init_permutation(dof_map.n_dofs());
}



//===========================================================//
void EnvelopFunctionApprox::do_solve()
{

  // check unused tags in solver_options
  const ModelOptions& sol_opt = get_solver_options();
  sol_opt.find_option("simulation"); // remove simulation name
  sol_opt.check_unused();

  if (_solution.size() !=
      get_solution_descriptor(EigenEnergy).n_components())
  {
    redeclare_solutions();
  }
   RealVectorValue k_vec(0.0);
   _calculate_density = true;
   if (has_option("k_vector"))
   {
     get_parameter("k_vector", k_vec);
     Messages::warning("k-vector given, will skip density calculation.");
     _calculate_density = false;
   }
   set_k_point(k_vec);

  if (_calculate_density && (get_k_point().norm() == 0.0))
  {
    estimate_spectrum_shift();
    //apply_bc();

    if (opt.analytic_k_int)
    {
      calculate_density_analytic();
    }
    else
    {
      DofField dens;
      integrate_density(dens);

      TiberLinearSystem& qdens_sys = get_equation_system<TiberLinearSystem>();
      NumericVector<Number>& qdens = *qdens_sys.solution;
      for (unsigned int i = 0; i < qdens.size(); ++i)
        qdens.set(i, dens[i]);

      qdens.close();
    }
  }
  else
  {
    solve_for_kpoint(get_k_point());
  }
}



void
EnvelopFunctionApprox::redeclare_solutions(void)
{
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
  declare_solution(EnvelopeFunctions, NTUPLE, NODES, "", num_states * 2 * number_of_bands);
  declare_solution(EigenEnergy, NTUPLE, GLOBAL, "eV", num_states);
  declare_solution(Occupation, NTUPLE, GLOBAL, "", num_states);
  declare_solution(EigenEnergyOnMesh, NTUPLE, NODES, "eV", num_states);
}




//===========================================================//
void EnvelopFunctionApprox::do_solve_for_kpoint(const Point& k_point)
{
   
  if (verbose() > 0)
  {
    ostringstream os;
    os << "(EFA) Solving for k = (";
    get_k_point().write_unformatted(os, false);
    os << ") /nm";
    Messages::info(os.str());
  }

  if ( _job == BULKEIGENSTATES )
  {
    solve_bulk();
  }
  else
  { 

    estimate_spectrum_shift();
    
    //apply_bc();
    

    initialize_solution_container(opt.num_el_states + opt.num_hl_states);
    
    solve_eigen_value_problem((opt.num_el_states + opt.num_hl_states) / 2 + 1,
                              solver_opt.spectrum_shift/Hartree);


    redeclare_solutions();
  }
}

//===========================================================//
void EnvelopFunctionApprox::calculate_Hamiltonian_and_S(void)
{

  apply_bc();

  _H_real->zero();
  _H_imag->zero();
  _S_real->zero();
 
  _haveS = (solver_opt.discretization_method == FEM);
  
  apply_bc();

  //material list
  //assemble_material_list();

 vector<unsigned int> psivar(number_of_bands);
 //get numbers of variables
 for (unsigned int i = 0; i < number_of_bands; i++)
 {
   psivar[i] = system->variable_number(psi_name[i]);
 }

 DofMap& dof_map = system->get_dof_map();

 FEType fe_type = dof_map.variable_type(psivar[0]); //all the variable have the same FE representation



 UniquePtr<libMesh::FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  // A 5th order Gauss quadrature rule for numerical integration.
  //QGauss qrule (dim, FIFTH);
 UniquePtr<libMesh::QBase> qrule (libMesh::QBase::build(_quadrature_type, dim, SECOND));

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

    MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
    for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;
      int n_nodes = elem->n_nodes();
      fe->reinit(elem);

      for (unsigned int b = 0; b < number_of_bands; b++)
      {
        dof_map.dof_indices(elem, dof_indices_component, psivar[b]);
        const unsigned int n_dofs = dof_indices_component.size();

        for (unsigned int qp=0; qp < (*qrule).n_points(); qp++)
          for (unsigned int p = 0; p < n_dofs; p++)
            _sqrt_S_inv[dof_indices_component[p]] += JxW[qp] * phi[p][qp] * phi[p][qp];
      }
    }

    this->get_solver_communicator().sum(_sqrt_S_inv);

    // build the square root and invert
    for (size_t i = 0; i < _sqrt_S_inv.size(); i++)
      _sqrt_S_inv[i] = 1.0 / sqrt(_sqrt_S_inv[i]);

  }


  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  double electric_potential = 0;

  EFAbulkHamiltonian* element_hamiltonian;

  // a temporary array to pass k point to models
  double k_vector[3] = {  get_k_point()(0) * 1e9 * Constants::bohr_radius,
                          get_k_point()(1) * 1e9 * Constants::bohr_radius,
                          get_k_point()(2) * 1e9 * Constants::bohr_radius};

  for ( ; el != end_el ; ++el)
    {//el
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      element_hamiltonian = get_bulk_model<EFAbulkModel>(elem)->get_Hamiltonian_model();

      element_hamiltonian->set_temperature(_temp_interface.get_temperature( elem, elem->centroid()));

      element_hamiltonian->set_k_vector(k_vector);

      element_hamiltonian->calculate_Hamiltonian_k_par();


      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs   = dof_indices.size();


      ham_real.resize(n_dofs, n_dofs);
      ham_imag.resize(n_dofs, n_dofs);
      if (_haveS) s_real.resize(n_dofs, n_dofs);

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




        for (unsigned int band1 = 0; band1 < number_of_bands; band1++)
        {//band1
          dof_map.dof_indices (elem, dof_indices_component, psivar[band1]);
          const unsigned int n_psi_dofs = dof_indices_component.size();

          for (unsigned int band2 = 0; band2 < number_of_bands; band2++)
          {//band2

            // Reposition the submatrix relative to the indices of the block
            // corresponding to band1-band2 inside the whole dense block
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
            if (_haveS)
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


      if (_quadrature_type == QTRAP)
      {
        // apply S^-1/2 H S^-1/2
        for (unsigned int i = 0; i < n_dofs; i++)
        {
          for (unsigned int j = 0; j < n_dofs; j++)
          {
            double scale = _sqrt_S_inv[dof_indices[i]] * _sqrt_S_inv[dof_indices[j]];
            ham_real(i, j) *= scale;
            ham_imag(i, j) *= scale;
            // this is not needed, as we do not solve a generalized problem
            // in this case
            if (_haveS) s_real(i, j) *= scale;
          }
        }
      }

      double penalty = 1e6;

      // penalty method for Dirichlet nodes
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        double sign = ham_real(i, i) / abs(ham_real(i,i));
        if (dirichlet_dofs.count(dof_indices[i]))
        {
          ham_real(i, i) += sign * penalty;

          //for (unsigned int j = 0; j < n_dofs; j++)
          //{
          //  ham_real(i, j) /= penalty;
          //  ham_imag(i, j) /= penalty;
          //  s_real(i, j) /= penalty;
          //}
        }
      }


      {
        vector<unsigned int> new_dof_indices;
        vector<unsigned int> dof_indices_tmp;
        
        //set<unsigned int> constrained;


        dof_indices_tmp = dof_indices;
        dof_map.constrain_element_matrix(ham_real, dof_indices_tmp, false);
        new_dof_indices.resize(dof_indices_tmp.size());

        dof_indices_tmp = dof_indices;
        dof_map.constrain_element_matrix(ham_imag, dof_indices_tmp, false);

        for (unsigned int i=0; i< new_dof_indices.size(); i++)
        {
          DofConstraints::iterator it = my_dof_constraints.find(dof_indices_tmp[i]);
          if (it != my_dof_constraints.end())
          {
            double arg = get_k_point() * get_periodicity_vector(dof_indices_tmp[i]);
            Complex phase = exp(Complex(0.0, 1.0)*arg);
            //cerr << "arg = " << arg << " " << get_k_point() << " " << get_periodicity_vector(dof_indices_tmp[i]) << endl;
            //constrained.insert(i);
            double sign = ham_real(i,i) / abs(ham_real(i,i));
            ham_real(i,i) += sign * penalty;
            for (int j = 0; j < ham_real.n(); j++)
            {
              DofConstraintRow::iterator constr(
                  (it->second).find(dof_indices_tmp[j]));
              if (constr != (it->second).end())
              {
                ham_real(i,j) -= sign * penalty * real(phase) * constr->second;
                ham_imag(i,j) -= sign * penalty * imag(phase) * constr->second;
              }
            }
          }

          new_dof_indices[i] = _perm[dof_indices_tmp[i]];
        }
        
        if (_haveS)
        {
          dof_indices_tmp = dof_indices;
          dof_map.constrain_element_matrix(s_real, dof_indices_tmp, false);
        }

        /*
        set<unsigned int>::iterator it = constrained.begin();
        for ( ; it != constrained.end(); ++it)
        {
          unsigned int i = *it;
          for (int j = 0; j < ham_real.n(); j++)
          {
            ham_real(i,j) /= penalty;
            ham_imag(i,j) /= penalty;
            if (_haveS)
              s_real(i,j) /= penalty;
          }
        }
        */

        _H_real->add_matrix(ham_real, new_dof_indices);
        _H_imag->add_matrix(ham_imag, new_dof_indices);
        if (_haveS)
          _S_real->add_matrix(s_real, new_dof_indices);
      }


    }

  
  _H_real->close();
  _H_imag->close();
  if (_haveS)
    _S_real->close();
  //_H_real->print_matlab("Hr.m");
  //_H_imag->print_matlab("Hi.m");
  //if (_haveS)
  //  _S_real->print_matlab("S.m");

  //  dof_map.print_dof_constraints();

}


//============================================================//
void EnvelopFunctionApprox::estimate_spectrum_shift(void)
{
 if (opt.estimate_spectrum_shift)
 {
   double Ec = get_band_edge("el");
   double Ev = get_band_edge("hl");

   // TODO treat electron and hole only case

   ostringstream os;
   os <<"Maximum of Ev (eV) :" <<Ev<<std::endl;
   os <<"Minimum of Ec (eV) :" <<Ec<<std::endl;
   Messages::info(os.str());
   os.str("");

   if ((Ec - Ev) <= 0.0)
   {
     Messages::warning("Your system apparently does not have a global gap: "
         "cannot find reasonable guess.");
     Messages::warning("Will use mean value of band edges");
     solver_opt.spectrum_shift = (Ec + Ev) / 2.0;
   }
   else
   {
     int states = opt.num_el_states + opt.num_hl_states;
     double frac = opt.num_el_states / states;
     double margin = 0.1 * (Ec - Ev);
     solver_opt.spectrum_shift = Ev + margin + frac * (Ec - 2*margin - Ev);
   }
   os<<"(EFA) Estimated guess (eV): " << solver_opt.spectrum_shift << std::endl;
   Messages::info(os.str());
 }

}
//=============================================================//
double EnvelopFunctionApprox::get_new_spectrum_shift(void)
{
/*
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
  */

  return solver_opt.spectrum_shift / Constants::Hartree;

}




bool EnvelopFunctionApprox::read_SLEPC_solution(void)
{//
  /*
  1) Read all eigenvalues
  2) Sort the eigenvalues and select those we want
  3) Read eigenvectors that correspond to the eigenvalues we want
  4) normalize eigenfunctions
  5) calculate fermi energy for each state
  */


  //--------------------------------------------------------------------
  //how many solutions do we have from SLEPC?
  unsigned int number_of_converged_solutions;

  number_of_converged_solutions = EigenSolver::number_of_converged_eigenvalues();

  //--------------------------------------------------------------------
  //read eigenvalues
  //store also eigenvalue index for sorting

  vector<EigenvalueProblem::eigen_state>  ev(number_of_converged_solutions);
  double shift = EigenSolver::get_shift() * Hartree;

  // if we have already solutions, we should use their energy levels to discriminate
  // between el and hl
  if (opt.num_hl_states < _solution.size() &&
      _solution[opt.num_hl_states].eigen_vector.size() > 0)
  {
    shift = _solution[opt.num_hl_states].eigen_energy;
  }
  else if ((opt.num_hl_states > 0) &&
           (_solution[opt.num_hl_states - 1].eigen_vector.size() > 0))
  {
    shift = _solution[opt.num_hl_states - 1].eigen_energy;
  }

  for (unsigned ind = 0; ind < number_of_converged_solutions; ind++)
  {
    ev[ind].energy =  EigenSolver::get_eigenvalue(ind) * Hartree; // + shift;
    ev[ind].index = ind;

    if (ev[ind].energy > shift)
      ev[ind].particle = "el";
    else
      ev[ind].particle = "hl";
  }

  // sorting of the solutions
  // we sort both electrons and holes by distance from the ground state

  sort(ev.begin(), ev.end(), compare_eigen_energy);

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


  // find the first electron state

  unsigned int first_el_index = number_of_converged_solutions;
  bool finish = false;

  for (unsigned int i = 0; i < number_of_converged_solutions; i++)
  {
    if (ev[i].particle == "el")
    {
      first_el_index = i;
      break;
    }
  }


  //--------------------------------------------------------------------
  //read eigenvectors

  // The idea is that arriving here the _solution structure is set up,
  // but all eigenvectors are empty or already calculated, valid eigenstates


  // the first num_hl_states in _solution are holes, the upper
  // num_el_states ones are electrons


  const int n_states = _solution.size();

  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < number_of_converged_solutions; i++)
  {
    // calculate the index in the final solution structure
    int index = static_cast<int>(opt.num_hl_states + i) - first_el_index;
    if (ev[i].particle == "hl")
      index = i;

    if (((ev[i].particle == "hl") && (index >= first_el_index)) ||
        ((ev[i].particle == "el") && (index < first_el_index)) ||
        ((index < 0) || (index >= n_states)))
        continue;


    // we need a small delta to decide if two states may be degenerate
    // TODO adjust it automatically
    const double delta = 1e-5;
    // if this is true, then we have to check linear dependency with neighbouring
    // states in the interval +/- delta
    int check_linear_dependency = false;

    // look for the first available slot
    if (ev[i].particle == "el")
    {
      while ((index < n_states) &&
             (_solution[index].eigen_vector.size() > 0))
      {
        index++;
      }

      if (index > opt.num_hl_states)
      {
        if ((index >= n_states) ||
            (ev[i].energy < (_solution[index - 1].eigen_energy - delta))) // go to the next state
          continue;


        if (ev[i].energy < (_solution[index - 1].eigen_energy + delta))
        {
          // we may have found a degenerate eigenvalue
          check_linear_dependency = true;
        }
      }
    }
    else // if (ev[i].particle == "hl")
    {
      while ((index < static_cast<int>(opt.num_hl_states)) &&
             (_solution[index].eigen_vector.size() > 0))
      {
        index++;
      }

      if (index > 0)
      {
        if ((index >= static_cast<int>(opt.num_hl_states)) ||
            (ev[i].energy > _solution[index - 1].eigen_energy + delta)) // go to the next state
          continue;


        if (ev[i].energy > _solution[index - 1].eigen_energy - delta)
        {
          // we may have found a degenerate eigenvalue
          check_linear_dependency = true;
        }
      }
    }
    //check_orthogonality = true;

    // we found a (potentially) valid slot and fill it
    _solution[index].eigen_energy = ev[i].energy;
    _solution[index].particle = ev[i].particle;
    _solution[index].statistics = "Fermi";
    //_solution[index].eigen_vector.resize(number_of_all_dofs, Complex(0.0, 0.0));
    
    //vector<Complex> temp;
    
    unsigned int solution_number = ev[i].index;

    //EigenSolver::get_eigen_vector(solution_number, temp);
    EigenSolver::get_eigen_vector(solution_number, _solution[index].eigen_vector);

    this->get_solver_communicator().allgather(_solution[index].eigen_vector);

    //-----------------------------------------------------------------------------
    //put independent dofs in the eigenvectors that may contain also non independent dofs
    for (int j = 0; j < number_of_all_dofs; j++)
    {
      //_solution[index].eigen_vector[j] = temp[j];
      //if (new_dofs[j].independent)
      //{
      //  _solution[index].eigen_vector[j] = temp[new_dofs[j].new_number];
      //}
    }

/*
    //put constrained dofs

    for (unsigned int j = 0; j < number_of_all_dofs; j++)
    {
      DofConstraints::iterator it(my_dof_constraints.find(j));

      if (it != my_dof_constraints.end() )
      {

        DofConstraintRow constr_row((it->second).first);

        DofConstraintRow::iterator  c =  constr_row.begin();

        for ( ; c != constr_row.end() ; ++c )
        {
          _solution[index].eigen_vector[j] += ( c->second ) *
              _solution[index].eigen_vector[(c->first)];
        }
      }
    }
*/

    //
    // apply transformation if needed
    //
    transform_eigenstate(_solution[index].eigen_vector);



    //
    // Now, we check that we did not take by chance a linearly dependent eigenstate
    // for a second time. If so, we delete the eigenvector and go to the
    // next state.
    //
    // TODO is this check correct? It should, because linearly independent eigenvectors
    //      should be so even if it is a generalized EVP. But actually, it should not be
    //      needed because now we add a deflation space.
    //
    if (check_linear_dependency)
    {
      vector<Complex> tempvec(_solution[index].eigen_vector);

      //cerr << "orthogonality (" << index << ") :";
      int ind = index - 1;

      while ((ind >= 0) && (ind < n_states))
      {
        //if (_solution[ind].eigen_vector.size() == 0) break;
        if ((_solution[index].eigen_energy < (_solution[ind].eigen_energy - delta)) ||
            (_solution[index].eigen_energy > (_solution[ind].eigen_energy + delta)))
          break;

        Complex norm = sqrt(scalar_product(_solution[ind], _solution[ind]));
        Complex alpha = scalar_product(_solution[ind].eigen_vector, tempvec) / norm;
        //cerr << " (" << ind << ") " << alpha << ",";
        for (size_t j = 0; j < number_of_all_dofs; j++)
          tempvec[j] -= alpha * _solution[ind].eigen_vector[j];

        ind--;
      }
      //cerr << "\n";

      double dotprod = abs(scalar_product(tempvec, tempvec));
      if (dotprod < 1e-6)
      {
        // we delete it and skip to the next state;
        _solution[index].eigen_vector.clear();
        continue;
      }
    }

    //
    //normalization
    //
    double norm = eigenstate_norm(index);

    for (unsigned int j = 0; j < number_of_all_dofs; j++)
      _solution[index].eigen_vector[j] /= Complex(norm, 0.0);



  }

  // the 1e-5 below is to not make the Hamiltonian singular,
  // and to be sure to take all states

  double Ec = shift;
  double Ev = shift;
  if (opt.estimate_spectrum_shift)
  {
    Ec = get_band_edge("el");
    Ev = get_band_edge("hl");
  }

  bool foundall = true;

  // did we find all electron eigenstates?
  int n_eig = opt.num_hl_states;
  for ( ; n_eig < n_states; n_eig++)
  {
    if (_solution[n_eig].eigen_vector.size() == 0)
    {
      foundall = false;
      break;
    }
    else
      solver_opt.spectrum_shift = _solution[n_eig].eigen_energy - 1e-5;
  }

  if (n_eig == opt.num_hl_states)
  {
    // in this case we found no electron state at all, so set shift near Ec
    if (opt.estimate_spectrum_shift)
    {
      solver_opt.spectrum_shift = Ec - 0.05;
      // If there is no gap, we leave the guess at the mean band edge energy.
      // Sooner or later we will find all states, and el/hl does not make
      // really sense here anyway.
      if ((Ec - Ev) <= 0.0)
        solver_opt.spectrum_shift = (Ec + Ev) / 2.0;
    }
    else
      solver_opt.spectrum_shift += 0.3;  //!? Rise a little the guess and restart  

  }
  // if not all are found, look for so many electron states:
  solver_opt.number_of_eigenstates = opt.num_el_states - (n_eig - opt.num_hl_states) + 1;

  if (foundall)
  {
    solver_opt.spectrum_shift = shift;

    // did we find all hole eigenstates?
    for (n_eig = static_cast<int>(opt.num_hl_states) - 1; n_eig >= 0; n_eig--)
    {
      if (_solution[n_eig].eigen_vector.size() == 0)
      {
        foundall = false;
        break;
      }
      else
        solver_opt.spectrum_shift = _solution[n_eig].eigen_energy + 1e-5;
    }

    if (n_eig == static_cast<int>(opt.num_hl_states) - 1)
    {
      // in this case we found no state at all, so set shift near Ev
      if (opt.estimate_spectrum_shift)
      {
        solver_opt.spectrum_shift = Ev + 0.05;
        if ((Ec - Ev) <= 0.0)
          solver_opt.spectrum_shift = (Ec + Ev) / 2.0;
      }
      else
        solver_opt.spectrum_shift -= 0.3;  //!? Lower a little the guess and restart

    }
    solver_opt.number_of_eigenstates = n_eig + 1;

    //redeclare_solutions();
  }

  if (foundall)
  {
    redeclare_solutions();

    for (unsigned int i = 0; i < n_states; i++)
    {
      // Fermi energy calculation
      _solution[i].electro_chem_pot = calculate_fermi_averaged(i,
          _solution[i].particle);


      //Temperature calculation
      _solution[i].temperature = calculate_temperature_averaged(i);

//    if (!check_confinement(_solution[i].eigen_vector))
//    {
//      ostringstream os;
//      os << "State " << i << " is not confined!";
//      Messages::warning(os.str());
//    }
    }
  }

  
  return foundall;

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

  UniquePtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  QGauss qrule (dim, SECOND);

  fe -> attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  std::vector<unsigned int> dof_indices;


  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  Complex sum(0.0, 0.0);

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    if (!reg_ids.count(elem->subdomain_id()))
      continue;

    fe->reinit (elem);

    for (short psi_index = 0; psi_index < number_of_bands; psi_index++)
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

  // sum up over MPI processes
  this->get_solver_communicator().sum(sum);



  double confinement_threshold = get_option("confinement_threshold", 0.2);
  confined = (sqrt(abs(sum)) > confinement_threshold) ? true : false;

  return confined;
}


/*
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
    file << "\n#\n";
    file << "# Index" << setw(12) << "EigenEnergy" << setw(15) << "Occupation"
        << setw(12) << "FermiLevel" << setw(12) << "Temperature" << "\n";

    for (unsigned int i = 0; i < _solution.size(); i++)
    {
      file << setw(7) << i << " "
          << setw(11) << _solution[i].eigen_energy << " "
          << setw(14) << Fermi_statistics_probability(_solution[i].eigen_energy,
              _solution[i].electro_chem_pot, _solution[i].temperature, _solution[i].particle) << " "
          << setw(11) << _solution[i].electro_chem_pot << " "
          << setw(11) << _solution[i].temperature << "\n";
    }
  }


}
*/


//=======================================================================//


EnvelopFunctionApprox:: ~EnvelopFunctionApprox(void)
{

  // es->delete_system(system_name);
}

//=======================================================================//



void EnvelopFunctionApprox::transform_eigenstate(vector<Complex>& eigvec)
{
  // if we use QTRAP (diagonal overlap), we must transform the eigenstate
  // with S^-1/2 as we transformed the system to have unit overlap
  if (_quadrature_type != QTRAP) return;
  if (eigvec.size() == 0) return;

  size_t n_dofs = eigvec.size();
  if (n_dofs != _sqrt_S_inv.size())
    cerr << "eigvec : " << eigvec.size() << " " << "S : " << _sqrt_S_inv.size() << endl;

  for (size_t i = 0; i < n_dofs; i++)
    eigvec[i] *= _sqrt_S_inv[i];

}

//-----------------------------------------------------------------------------//
double  EnvelopFunctionApprox::eigenstate_norm(unsigned int state_number)
{
  double  result;

  const vector< Complex > &  eigen_vector =  _solution[state_number].eigen_vector;



  DofMap& dof_map = system->get_dof_map();

  FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

  // UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  //UniquePtr<FEBase> fe (  build_finite_element(dim, fe_type)  );
  UniquePtr<libMesh::FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  QGauss qrule(dim, SECOND);

  fe -> attach_quadrature_rule (&qrule);


  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  std::vector<unsigned int> dof_indices;


  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  Complex temp(0.0, 0.0);
  Complex eigen_f_value1, eigen_f_value2;

  for ( ; el != end_el ; ++el)
    {//el

      const Elem* elem = *el;
      fe->reinit (elem);

      for (short psi_index = 0; psi_index < number_of_bands; psi_index++)
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


  this->get_solver_communicator().sum(temp);

  result = sqrt( abs(temp)  );




  return(result);

}



//==========================================================//


double EnvelopFunctionApprox::calculate_fermi_averaged(unsigned int i, const string& particle)
{

  if (((particle == "el") && (_el_elchem.second == INVALID_ID)) ||
      ((particle == "hl") && (_hl_elchem.second == INVALID_ID)))
    return 0.0;

  Complex  result(0.0, 0.0);

  const vector< Complex >&   eigen_vector =  _solution[i].eigen_vector;


  const MeshBase* mesh = &(es->get_mesh());


  unsigned int dim = mesh->mesh_dimension();


  system = &( es->get_system<LinearImplicitSystem>(system_name));

  DofMap& dof_map = system->get_dof_map();


  FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

   //  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
   UniquePtr<libMesh::FEBase> fe (build_finite_element(dim, fe_type, true));

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


  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();


  Complex eigen_f_value1;
  Complex eigen_f_value2;


  double chem_pot_value_eV;

  for ( ; el != end_el ; ++el)
    {//el

      const Elem* elem = *el;
      fe->reinit (elem);


      Point center = elem->centroid();

      if (particle == "el")
        chem_pot_value_eV = get_el_electro_chem_potential(elem);
      else
        chem_pot_value_eV = get_hl_electro_chem_potential(elem);



      for (short psi_index = 0; psi_index < number_of_bands; psi_index++)
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

  this->get_solver_communicator().sum(result);

  return(result.real());
}




double EnvelopFunctionApprox::calculate_temperature_averaged(unsigned int i)
{

  Complex result(SimulationOptions::temperature,0.0);


  //-----------------------------------------------------//


  if (_temp_interface.has_simulation())
  {

    result = 0.0;

    const vector< Complex >&   eigen_vector =  _solution[i].eigen_vector;

    //----------------------------------------------------//


    const MeshBase* mesh = &(es->get_mesh());

    unsigned int dim = mesh->mesh_dimension();

    system = &( es->get_system<libMesh::LinearImplicitSystem>(system_name));

    libMesh::DofMap& dof_map = system->get_dof_map();

    FEType fe_type = dof_map.variable_type(0); //all the variable have the same FE representation

    //  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    UniquePtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

    // A 5th order Gauss quadrature rule for numerical integration.
    QGauss qrule (dim, libMesh::SECOND);

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


    MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();


    Complex eigen_f_value1;
    Complex eigen_f_value2;


    for ( ; el != end_el ; ++el)
    {//el

      const Elem* elem = *el;
      fe->reinit (elem);


      Point center = elem->centroid();
      double Temperature = _temp_interface.get_temperature(elem, center);




      for (short psi_index = 0; psi_index < number_of_bands; psi_index++)
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
    this->get_solver_communicator().sum(result);
  }

  return(result.real());
}




//--------------------------------------------------------------------------//
void EnvelopFunctionApprox::do_assemble(const ModelOptions&)
{
  calculate_Hamiltonian_and_S();
}




//===========================================================//





//=================================================================//

  
void EnvelopFunctionApprox::calculate_density_analytic(void)
{
  unsigned int dim = get_mesh().mesh_dimension();


  // we solve with +1 state to be used in dd to define the classical boundary
  if (opt.num_el_states > 0) opt.num_el_states++;
  if (opt.num_hl_states > 0) opt.num_hl_states++;


  unsigned int num_states = opt.num_el_states + opt.num_hl_states;


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

  vector<double> energy_k_0;
  vector<double> energy_k_1;
  vector<double> energy_k_2;
  vector<double> energy_k_3;
  vector<double> effective_mass(num_states);

  double spectrum_shift = solver_opt.spectrum_shift/Constants::Hartree;

  // [0 0 1]
  if (dim < 3)
  {
    ostringstream os;
    os << "Solving for k1 = ( ";
    kvector_1.write_unformatted(os, false);
    os << ")";
    if (verbose() > 0) m.info(os.str());

    set_k_point(kvector_1);
    initialize_solution_container(num_states);
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
    initialize_solution_container(num_states);
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
      initialize_solution_container(num_states);
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
  initialize_solution_container(opt.num_el_states + opt.num_hl_states);
  solve_eigen_value_problem(opt.num_el_states + opt.num_hl_states, spectrum_shift);
  get_eigenenergies(energy_k_0);


  // 2014-01-07: before, k-points were measured in Bohr radii, but now they
  //             are in nm!
  //double Eh_k2 = Constants::Hartree * opt.k_val * opt.k_val;
  double Eh_k2 = opt.k_val * opt.k_val *
      (hbar * hbar) / (electron_mass * 1e-18 * elementary_charge);

  if (dim == 1)
  {
    for (unsigned int i = 0; i < num_states; i++)
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
    for (unsigned int i = 0; i < num_states; i++)
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
  UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));

  DofMap& dof_map = system->get_dof_map();
  std::vector<unsigned int> dof_indices;

  // The qdens_sys system contains the nodal quantum density
  TiberLinearSystem& qdens_sys = get_equation_system<TiberLinearSystem>();
  DofMap& dof_map_qdens = qdens_sys.get_dof_map();
  std::vector<unsigned int> dof_indices_qdens;
  std::vector<unsigned int> dof_indices_qdens_p;
  NumericVector<Number>& qdens = *qdens_sys.solution;
  qdens.zero();

  // we need the connectivity of the nodes to not double count
  vector<int> connectivity(qdens.size(), 0.0);
  {
    MeshBase::const_element_iterator el = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map_qdens.dof_indices(elem, dof_indices_qdens, 0);
      dof_map_qdens.dof_indices(elem, dof_indices_qdens_p, 1);
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        connectivity[dof_indices_qdens[n]]++;
        connectivity[dof_indices_qdens_p[n]]++;
      }
    }
  }

  this->get_solver_communicator().sum(connectivity);

  unsigned int start = 0;
  unsigned int stop  = num_states - 1;

  // this is for the length scaling, EFA uses Bohr radii internally
  double a_B =  Constants::bohr_radius;
  double scaling =  1.0 / (a_B * a_B * a_B * 1.0e6 );

  for (unsigned int i = start; i < stop; i++)
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
      if (_solution[i].particle == "hl")
        exp_arg = -exp_arg;

      dos_factor = (exp_arg < -20) ? exp(exp_arg) : log(1.0 + exp(exp_arg));
    }
    else if (dim == 2)
    {
      // 2D is correct and tested
      mass_factor = sqrt(kT *  effective_mass[i] / (2.0 * M_PI * Constants::Hartree));
      double exp_arg = (fermi_energy - energy) / kT;
      if (_solution[i].particle == "hl")
        exp_arg = -exp_arg;

      dos_factor = TiberMath::fermidirac_mhalf(exp_arg);
    }
    else if (dim == 3)
      dos_factor = Fermi_statistics_probability(energy, fermi_energy,
          _solution[i].temperature, _solution[i].particle);

    dos_factor *= mass_factor * scaling * opt.degeneracy;

    unsigned int particletype = 0;
    if (_solution[i].particle == "hl")
      particletype = 1;

    MeshBase::const_element_iterator el = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map_qdens.dof_indices(elem, dof_indices_qdens, particletype);

      for (int psi_index = 0; psi_index < number_of_bands; psi_index++)
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
  qdens_sys.update();

  redeclare_solutions();

  // set them back to original values
  if (opt.num_el_states > 0) opt.num_el_states--;
  if (opt.num_hl_states > 0) opt.num_hl_states--;
}




void
EnvelopFunctionApprox::do_calculate_density_at_k(DofField& density)
{
  FEType fe_type = system->variable_type(0);
  UniquePtr<FEBase> fe(build_finite_element(dim, fe_type));

  DofMap& dof_map = system->get_dof_map();
  std::vector<unsigned int> dof_indices;

  // The qdens_sys system contains the nodal quantum density
  TiberLinearSystem& qdens_sys = get_equation_system<TiberLinearSystem>();
  DofMap& dof_map_qdens = qdens_sys.get_dof_map();
  std::vector<unsigned int> dof_indices_qdens;
  //std::vector<unsigned int> dof_indices_qdens_p;
  NumericVector<Number>& qdens = *qdens_sys.solution;
  //qdens.zero();

  // we need the connectivity of the nodes to not double count
  vector<int> connectivity(qdens.size(), 0.0);
  {
    MeshBase::const_element_iterator el = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map_qdens.dof_indices(elem, dof_indices_qdens, 0);
      //dof_map_qdens.dof_indices(elem, dof_indices_qdens_p, 1);
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        connectivity[dof_indices_qdens[n]]++;
        //connectivity[dof_indices_qdens_p[n]]++;
      }
    }
  }

  this->get_solver_communicator().sum(connectivity);

  density.clear();
  density.resize(qdens.size(), 0.0);

  // this is for the length scaling, EFA uses Bohr radii internally
  double a_B =  Constants::bohr_radius;
  double scaling = opt.degeneracy / 1e6;
  switch (dim)
  {
    case 3:
      scaling /= a_B;
    case 2:
      scaling /= a_B;
    case 1:
      scaling /= a_B;
    default:
      break;
  }

  // k-space uses nm as units
  switch (get_kspace()->dimension())
  {
    case 3:
      scaling *= 1e9;
    case 2:
      scaling *= 1e9;
    case 1:
      scaling *= 1e9;
    default:
      break;
  }


  for (unsigned int i = 0; i < _solution.size(); i++)
  {
    double fermi_energy = _solution[i].electro_chem_pot;
    double kT = _solution[i].temperature * Constants::k_Boltzmann;

    double energy = _solution[i].eigen_energy;

    unsigned int particletype = 0;
    if (_solution[i].particle == "hl")
      particletype = 1;

    double dos_factor = Fermi_statistics_probability(energy, fermi_energy,
        _solution[i].temperature, _solution[i].particle);

    MeshBase::const_element_iterator el = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map_qdens.dof_indices(elem, dof_indices_qdens, particletype);

      for (int psi_index = 0; psi_index < number_of_bands; psi_index++)
      {
        dof_map.dof_indices(elem, dof_indices, psi_index);

        for (unsigned int n = 0; n < elem->n_nodes(); n++)
        {
          double psi = abs(_solution[i].eigen_vector[dof_indices[n]]);
          double val = scaling * dos_factor * psi * psi / connectivity[dof_indices_qdens[n]];
          density[dof_indices_qdens[n]] += val;
        }
      }
    }
  }
}





//==============================================================================//

unsigned int EnvelopFunctionApprox::get_number_of_active_cells()
{
  unsigned int result = 0;
  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  for ( ; el != end_el ; ++el)
    result++;

  this->get_solver_communicator().sum(result);

  return(result);



}

//===============================================================================//
short EnvelopFunctionApprox::calculate_number_of_bands(void) const
{

  // let -1 indicate absence of the model
  short result = -1;

  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    if (!this->get_environment().contains_region(elem->subdomain_id()))
      continue;

    EFAbulkHamiltonian* element_hamiltonian = get_bulk_model<EFAbulkModel>(elem)->get_Hamiltonian_model();


    const std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
      model_Ham = ( element_hamiltonian->get_Hamiltonian() );


    short result1 = model_Ham.size();

    if (result1 != result)
    {
      if (result > -1)
      {
        std::string mess = "EnvelopFunctionApprox: Hamiltonians of different"
            " materials have different number of bands";
        throw InitFailedException(mess);
      }
      else
        result = result1;
    }

  }

  vector<short> num_bands(this->get_solver_communicator().size());
  this->get_solver_communicator().allgather(result, num_bands);

  for (unsigned int i = 0; i < num_bands.size(); ++i)
  {
    if (result != num_bands[i])
    {
      if (result == -1)
        result = num_bands[i];

      break;
    }
  }

  if (!this->get_solver_communicator().verify(result))
  {
    std::string mess = "EnvelopFunctionApprox: Hamiltonians have different"
       " different number of bands on different mesh pieces";
    throw InitFailedException(mess);
  }

  return(result);

}


//=======================================================================================//

void EnvelopFunctionApprox::solve_bulk(void)
{

  Point qp = _bulk_mat_element->centroid();

  EFAbulkHamiltonian* element_hamiltonian;

  //std::cout<<"elem H"<<std::endl;

  element_hamiltonian = get_bulk_model<EFAbulkModel>(_bulk_mat_element)->get_Hamiltonian_model();

  // a temporary array to pass k point to models, ASSUMING k TO BE IN nm
  double k_vector[3] = { get_k_point()(0) * 1e9 * Constants::bohr_radius,
                         get_k_point()(1) * 1e9 * Constants::bohr_radius,
                         get_k_point()(2) * 1e9 * Constants::bohr_radius};

  element_hamiltonian->set_k_vector(k_vector);

  element_hamiltonian->calculate_Hamiltonian_k_par();

  Tensor2Sym strain_crystal_system(0);
  _strain_interface.get_crystal_strain(_bulk_mat_element, qp, strain_crystal_system);

  //std::cout<<"strain"<<std::endl;
  //std::cout<<"(EP) strain exx "<<strain_crystal_system(1,1)<<std::endl;
  //std::cout<<"(EP) strain eyy "<<strain_crystal_system(2,2)<<std::endl;
  //std::cout<<"(EP) strain ezz "<<strain_crystal_system(3,3)<<std::endl;


  double electric_potential = 0;
  if (opt.consider_potential_bulk)
  {
    electric_potential = get_electric_potential(_bulk_mat_element, qp );
    estimate_spectrum_shift();
  }


  element_hamiltonian->apply_strain_and_potential(strain_crystal_system, electric_potential);


  std::vector<std::vector<EFAbulkHamiltonian::MatrixElement> >&
	    model_Ham = ( element_hamiltonian->get_Hamiltonian() );


  std::complex<double> ham_matrix[number_of_bands * number_of_bands ];

  for (unsigned int band1 = 0; band1 < number_of_bands; band1++)
    for (unsigned int band2 = 0; band2 < number_of_bands; band2++)
    {
      ham_matrix[band1 + band2 * number_of_bands] = model_Ham[band1][band2].constant;

    }


  initialize_solution_container(number_of_bands);

  {
     char jobs = 'V';
     char UPLO = 'U';
     int  N = number_of_bands;
     double eigvals[number_of_bands];
     std::complex<double> WORK[2*number_of_bands-1];
     int LWORK = 2*number_of_bands-1;
     double RWORK[3*number_of_bands-2];
     int info;

     zheev_(jobs, UPLO, N, ham_matrix, N, eigvals, WORK, LWORK, RWORK, info);

     if (info !=0 ) throw SolveFailedException("LAPACK problem\n");;

     for (short i = 0; i < number_of_bands ; i++)
     {

       _solution[i].eigen_energy = eigvals[i]*Hartree;
       _solution[i].statistics = "Fermi";
       if (_solution[i].eigen_energy < solver_opt.spectrum_shift)
         _solution[i].particle = "hl";
       else
         _solution[i].particle = "el";
       _solution[i].eigen_vector.resize(number_of_bands);
       for (short j = 0; j < number_of_bands ; j++)
       {
         _solution[i].eigen_vector[j] = ham_matrix[i * number_of_bands + j];
       }
     }
  }

  unsigned int n = _solution.size();

  // reorder
  int max_hl = 0;

  double el_chem = 0;
  double hl_chem = 0;
  double temp = _temp_interface.get_temperature(_bulk_mat_element, _bulk_mat_element->centroid());

  if (_el_elchem.second != INVALID_ID)
    el_chem = get_el_electro_chem_potential(_bulk_mat_element);

  if (_hl_elchem.second != INVALID_ID)
    hl_chem = get_hl_electro_chem_potential(_bulk_mat_element);

  for (unsigned int i = 0; i < n ; i++)
  {
    _solution[i].temperature = temp;

    if (_solution[i].particle == "el")
      _solution[i].electro_chem_pot = el_chem;
    else
    {
      _solution[i].electro_chem_pot = hl_chem;
      max_hl = i;
    }
  }

  // reorder according to energy
  for (unsigned int i = 0; i <= max_hl ; i++)
  {
    unsigned int max_id = i;
    for (unsigned int j = i+1; j <= max_hl; j++)
      if (_solution[j].eigen_energy > _solution[max_id].eigen_energy)
        max_id = j;

    swap(_solution[i], _solution[max_id]);
  }


  // we have to redeclare the solution variables to adjust the number
  // of eigenstates
  redeclare_solutions();

}
