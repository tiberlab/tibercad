// $Id: SaveSolution.C 3735 2014-01-23 09:04:13Z maufder $

#include "SaveSolution.h"
#include "Ramp.h"
#include "Utils.h"
#include "Variable.h"
#include "Messages.h"
#include "TiberNonlinearSystem.h"
#include "SimulationEnvironment.h"

// libmesh includes
#include "node.h"
#include "mesh.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "quadrature_trap.h"
#include "equation_systems.h"
#include "mesh_refinement.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"

#include <fstream>



using namespace std;
using namespace libMesh;

SaveSolution::~SaveSolution(void)
{
}


void
SaveSolution::do_init(void)
{
  set<ID> save_id, sim_id;

  get_region_ids(save_id);
  _simulation->get_region_ids(sim_id);

  if ( save_id != sim_id )
  {
     string msg("Regions for " + get_name() + " and " + _simulation->get_name() + " need to be the same");
     throw InitFailedException(msg);
  }

  build_equation_system();
}


void
SaveSolution::do_print_info(void)
{
  string msg("Saving solutions for simulation: " + _simulation->get_name());
  Messages::info(msg);
}


void
SaveSolution::parse_options(void)
{
  ModelOptions& opts = get_options();

  // get the name of the simulation whose solutions are to be saved
  string sim;
  sim = opts.get_option("simulation", "");

  _simulation = SimulationInterface::find_simulation(sim);

  if (_simulation == NULL)
  {
    string msg("To be saved simulation " + sim + " not found");
    throw InitFailedException(msg);
  }

  // get the list of solutions to be saved
  vector<string> sols;
  opts.get_option("solutions", sols);

        vector<string>::iterator sol = sols.begin();
  const vector<string>::iterator sol_end = sols.end();

  for ( ; sol != sol_end; ++sol )
  {
    string solution = *sol;

    ID id = _simulation->get_solution_id(solution);
    if ( id == INVALID_ID )
    {
      string msg(solution + " is not a valid solution name for simulation " + sim);
      throw InitFailedException(msg);
    }

    const SolutionDescriptor& sd = _simulation->get_solution_descriptor(id);
    if ( ( sd.type() != SolutionDescriptor::REAL ) || ( sd.location() != SolutionDescriptor::NODES ) )
    {
      string msg(solution + " is not a REAL solution on NODES \nFor now SaveSolution module works only with REAL NODES");
      throw InitFailedException(msg);
    }

    _solution_ids.push_back(id);
  }
}


void
SaveSolution::do_solve(void)
{
  TiberNonlinearSystem& eq = get_equation_system<TiberNonlinearSystem>();
  const unsigned int dummy_var = eq.variable_number("dummy");
  const DofMap& dof_map = eq.get_dof_map();
  vector<unsigned int> dof_indices;

  map<ID, NumericVector<Number>&> solution_map; 

  vector<ID>::iterator sol = _solution_ids.begin();
  const vector<ID>::iterator sol_end = _solution_ids.end();

  for ( ; sol != sol_end; ++sol )
  {
    ID sol_id = *sol;
    const SolutionDescriptor& sd = _simulation->get_solution_descriptor(sol_id);
    solution_map.insert( pair<ID, NumericVector<Number>&>( sol_id, eq.get_vector(sd.name()) ) );
  }

        MeshBase::const_element_iterator el     = get_mesh().active_elements_begin();
  const MeshBase::const_element_iterator end_el = get_mesh().active_elements_end();


  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices, dummy_var);

    for (int i = 0; i < elem->n_nodes(); i++)
    {
            map<ID, NumericVector<Number>&>::iterator sol_map_it = solution_map.begin();
      const map<ID, NumericVector<Number>&>::iterator sol_map_it_end = solution_map.end();

      for ( ; sol_map_it != sol_map_it_end; ++sol_map_it )
      {
        double x = 0.0;
        _simulation->get_solution(elem, sol_map_it->first, x, elem->point(i));
        sol_map_it->second.set(dof_indices[i], x);
      }
    }
  }
}


void
SaveSolution::get_solution_secure(const Elem* elem,
                                        map<ID, std::vector<double> >& values,
                                  const vector<Point>& points)
{
  TiberNonlinearSystem& eq = get_equation_system<TiberNonlinearSystem>();
  const unsigned int dummy_var = eq.variable_number("dummy");
  const DofMap& dof_map = eq.get_dof_map();
  vector<unsigned int> dof_indices;

  const unsigned int dim = get_mesh().mesh_dimension();

  map<ID, NumericVector<Number>&> solution_map;

  vector<ID>::iterator sol = _solution_ids.begin();
  const vector<ID>::iterator sol_end = _solution_ids.end();

  for ( ; sol != sol_end; ++sol )
  {
    ID sol_id = *sol;
    const SolutionDescriptor& sd = _simulation->get_solution_descriptor(sol_id);
    solution_map.insert( pair<ID, NumericVector<Number>&>( sol_id, eq.get_vector(sd.name()) ) );
    //cout<<"Sim: "<<get_name()<<" eq.get_vector(sd.name()).size() = "<<eq.get_vector(sd.name()).size()<<endl;
    //cout<<"sol_id = "<<sol_id<<" solution_map.find(sol_id)->first = "<<solution_map.find(sol_id)->first<<" solution_map.find(sol_id)->second.size() = "<<solution_map.find(sol_id)->second.size()<<endl;
  }


  FEType fe_type = eq.variable_type(dummy_var);
  std::unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));

  const vector<vector<Real> >& phi = fe->get_phi();

  fe->reinit(elem, &points);

  dof_map.dof_indices(elem, dof_indices, dummy_var);

  unsigned int np = points.size();
  const unsigned int n_dofs = dof_indices.size();

  for (unsigned int n = 0; n < np; n++)
  {
    map<ID, std::vector<double> >::iterator values_it = values.begin();
    const map<ID, std::vector<double> >::iterator values_it_end = values.end();

    values_it = values.begin();
    for ( ; values_it != values_it_end; ++values_it )
    {
      //cout<<"values_it->first = "<<values_it->first<<"  solution_map.find(values_it->first)->first = "<<solution_map.find(values_it->first)->first<<endl;
      NumericVector<Number>& solution = solution_map.find(values_it->first)->second;

      //cout<<"values_it->first = "<<values_it->first<<" solution_map.find(values_it->first)->second.size() = "; //cout<<solution_map.find(values_it->first)->second.size()<<endl;

      double x = 0.0;

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        x += phi[i][n] * solution(dof_indices[i]);
      }
      values_it->second[n] = x;
    }


  }
}


void
SaveSolution::do_setup_solution_variables(void)
{
  parse_options();

        vector<ID>::iterator sol = _solution_ids.begin();
  const vector<ID>::iterator sol_end = _solution_ids.end();

  for ( ; sol != sol_end; ++sol )
  {
    ID sol_id = *sol;
    const SolutionDescriptor& sd = _simulation->get_solution_descriptor(sol_id);
    declare_solution_ext(sd.name(), sol_id, sd.type(), sd.location(), sd.units());
  }
}


void
SaveSolution::build_equation_system(void)
{
  //MeshBase& mesh = _simulation->get_mesh();
  //set_mesh(&mesh);

  create_equation_system("nonlinear");
  TiberNonlinearSystem& eq = get_equation_system<TiberNonlinearSystem>();

  eq.add_variable("dummy", libMeshEnums::FIRST);

  vector<ID>::iterator sol = _solution_ids.begin();
  const vector<ID>::iterator sol_end = _solution_ids.end();

  for ( ; sol != sol_end; ++sol )
  {
    ID sol_id = *sol;
    const SolutionDescriptor& sd = _simulation->get_solution_descriptor(sol_id);
    eq.add_vector(sd.name());
  }

  eq.init();
}
