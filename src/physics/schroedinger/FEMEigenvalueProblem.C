// $Id$

#include "tibercad/physics/schroedinger/FEMEigenvalueProblem.h"
#include "tibercad/geom/Boundary.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/solver/slepc/EigenSolver.h"
#include "tibercad/io/Messages.h"

#include "libmesh/equation_systems.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/fe_interface.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/petsc_matrix.h"

#include <numeric>

using std::set;
using std::vector;

using namespace libMesh;


namespace
{

  class Compare
  {
    public:
      Compare(const PetscInt *c,
              const std::vector<unsigned int>& v) : _c(c), _v(v) {};

      bool operator()(unsigned int i, unsigned int j)
      {
        return(_v[_c[i]] < _v[_c[j]]);
      }

    private:
      const PetscInt *_c;
      const std::vector<unsigned int>& _v;
  };

}

FEMEigenvalueProblem::FEMEigenvalueProblem(const ModelOptions& options)
 : EigenvalueProblem(options),
   min_coord {std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max()},
   max_coord {std::numeric_limits<double>::min(),
              std::numeric_limits<double>::min(),
              std::numeric_limits<double>::min()}
{

  es = NULL;

  system = NULL;

  _hamiltonian_size = 0;

  // dummy read
  get_solver_options().get_option("use_deflation_space", true);
}

//===============================================================//
void FEMEigenvalueProblem::make_new_dofs( )
{
  new_dofs.clear();


  DofMap& dof_map = system->get_dof_map();

  
  int number_of_all_dofs  =   dof_map.n_dofs();
  
  new_dofs.resize(number_of_all_dofs);
  const std::set<unsigned int> :: const_iterator  n_begin = dirichlet_dofs.begin();
  const std::set<unsigned int> :: const_iterator  n_end   = dirichlet_dofs.end();
  std::set<unsigned int> :: const_iterator n_it;

  unsigned int number_it = 0;


  //for (unsigned int i = 0; i < number_of_all_dofs ; i++)
  for (unsigned int i = dof_map.first_dof(); i < dof_map.end_dof() ; i++)
    {
//      if ( !( dof_map.is_constrained_dof(i) ) && (find(n_begin, n_end, i) == n_end))
	{
	    new_dofs[i].independent = true;
	    new_dofs[i].new_number = i;
	    //new_dofs[i].new_number = number_it;
	    number_it++;

	  }
//	else
//	  {
	  
//	    new_dofs[i].independent = false;
//	  }
    }

  

}
//=====================================================================//
void FEMEigenvalueProblem::make_constraints(void)
{
 
  DofMap& dof_map = system->get_dof_map();
  
 
  //----------------------------------------------------------------------------//
  //I recalculate my copy of the dof constraints because I need them explicitely!
  //  my_dof_constraints.clear(); not clear!!!

  // Look at all the variables in the system
  for (unsigned int variable_number=0; variable_number < dof_map.n_variables();
        ++variable_number)
  {
      

    MeshBase::const_element_iterator       elem_it  = this->active_local_elements_begin();
    const MeshBase::const_element_iterator elem_end = this->active_local_elements_end();
      
    for ( ; elem_it != elem_end; ++elem_it)
      FEInterface::compute_constraints (my_dof_constraints,
					dof_map,
					variable_number,
					*elem_it);
  }
 
  //-----------------------------------------------------------------------//
  //TODO: add periodic boundary conditions constraints

}

//=========================================================================//
//======================================================================//
void  FEMEigenvalueProblem::create_dirichlet_dofs( )
{
  
  

  SimulationEnvironment& se = get_environment(); 

  DofMap& dof_map = system->get_dof_map();

  MeshBase::const_element_iterator it = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end =  this->active_local_elements_end();

  dirichlet_dofs.clear();

  unsigned int number_of_variables = dof_map.n_variables();
 
  std::vector<unsigned int> dof_indices;

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    { 
      unsigned int  node_id =  elem->node_id(n);
      
      const Node* nd = elem->node_ptr(n); 

      Boundary* bd = se.get_boundary(nd); 
      
      //does a node belong to a a dirichlet boundary condition
      
      //   if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) )
       if (  (bd != NULL  ) )
	if ( bd->get_name() == "Dirichlet" || bd->get_name() == "dirichlet" )
	{
	  
	  for (short band = 0 ; band < number_of_variables ; band++)
	  {
	    dof_map.dof_indices (elem, dof_indices,band); 
	    dirichlet_dofs.insert(dof_indices[n]);
	  }
	     
	}

	  
    }
      
  }
  

}

//=======================================================================//
void FEMEigenvalueProblem::apply_dirichlet_at_all_boundaries()
{
  MeshBase::const_element_iterator it = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end =  this->active_local_elements_end();
 
 

  dirichlet_dofs.clear();

  

  DofMap& dof_map = system->get_dof_map();

  unsigned int number_of_variables = dof_map.n_variables();

  std::vector<unsigned int> dof_indices;

  

  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    unsigned int n_sides;

    // should now work always
    //if ( dim > 1 ) 
      n_sides = elem->n_sides();
    //else
    //  n_sides = elem->n_nodes();
     
	 	
    for (short i = 0; i < n_sides; i++)
    {
	 
      //check if the side is external -------------------------
      
      const Elem* el1 = elem->neighbor_ptr(i);

      bool side_is_external = false;

      if (get_environment().is_outer_boundary(ElementSide(elem, i)))
        side_is_external = true;
      
      /*
      if (el1 == NULL) 
      {
	side_is_external = true;

       
      }
      else
      {


	std::vector< const Elem * > active_family;
	
	if ( el1->has_children() )
	{
	  el1->active_family_tree (active_family);
	  
	  if (active_family.size() == 0)
	    side_is_external = true;
	  
	  //TODO
	  // 
	  //	has to be corrected because it may contain active child that does not belong to boundary
	  // 

	}
	else
	{//no children
	  if ( !(el1->active()) )
	    side_is_external = true;
	}



      }
      */
  

      //-------------------------------------------------------
	



      if (   side_is_external   )  
      {
	if (dim > 1)
	{//2D/3D
	  for (unsigned int nd = 0; nd < elem->n_nodes(); nd++)
	  {
	    if (elem->is_node_on_side(nd, i))
	    {
	      for (short band = 0 ; band <  number_of_variables; band++)
	      {
		dof_map.dof_indices (elem, dof_indices,band); 
		dirichlet_dofs.insert(dof_indices[nd]);
	      }
	    }
	  }
	  
	}
	else
	{//1D
	  for (short band = 0 ; band < number_of_variables; band++)
	  {
	    dof_map.dof_indices (elem, dof_indices,band);
	    dirichlet_dofs.insert(dof_indices[i]);
	  }
	}
      }      
    }   
  }
}



//=====================================================//
void FEMEigenvalueProblem::parse_options()
{
  ModelOptions& sol_opt = get_solver_options();

  solver_opt.number_of_eigenstates   = sol_opt.get_option("number_of_eigenstates", 6);

  solver_opt.spectrum_shift = sol_opt.get_option("guess", 0.0);

  solver_opt.Dirichlet_bc_everywhere = sol_opt.get_option("Dirichlet_bc_everywhere", true);
  solver_opt.Dirichlet_bc_everywhere = get_option("Dirichlet_bc_everywhere",
      solver_opt.Dirichlet_bc_everywhere);

  {
    std::string  method_name = get_option("discretization_method","FEM");
    if (method_name == "FEM")
      solver_opt.discretization_method = FEM;
    else if (method_name == "BIM")
      solver_opt.discretization_method = BIM;
    else
      throw InitFailedException( "FEMEigenvalueProblem: Incorrect method " + method_name );  

   
  }

  std::string preconditioner("lu");
  if (dim > 1)
    preconditioner = std::string("jacobi");

  preconditioner =  sol_opt.get_option("pc_type", preconditioner);
  sol_opt.set_option("pc_type", preconditioner);

}


void FEMEigenvalueProblem::copy_matrix_to_solver(const char matrix)
{

  // our problems are square
  int size_matrix = _H_real->n();
  _hamiltonian_size = size_matrix;
  
  PetscMatrix<Number>* H_real = nullptr;
  PetscMatrix<Number>* H_imag = nullptr;

  if (matrix == 'H')
  {
    H_real = static_cast<PetscMatrix<Number>* >(_H_real);
    H_real->close();

    H_imag = static_cast<PetscMatrix<Number>* >(_H_imag);
    H_imag->close();
  }
  else if (matrix == 'S')
  {
    H_real = static_cast<PetscMatrix<Number>* >(_S_real);
    H_real->close();
  }


  // the first local index
  const int start = H_real->row_start();
  const int end   = H_real->row_stop();

  const int local_size = end - start;

  // diagonal nonzeros
  vector<int> d_nnz(local_size);
  // off-diagonal nonzeros
  vector<int> o_nnz(local_size);

  // we strongly assume H and S have the same partitioning
  const int offset = H_real->row_start();

  for (int i = start; i < end; ++i)
  {
    PetscInt ncols_r = 0, ncols_i = 0;
    const PetscInt *row_r, *row_i;
    PetscErrorCode ierr;
    //const PetscReal *values_r, *values_i;

    ierr = MatGetRow(H_real->mat(), i, &ncols_r, &row_r, PETSC_NULLPTR);

    if (matrix == 'H')
      ierr = MatGetRow(H_imag->mat(), i, &ncols_i, &row_i, PETSC_NULLPTR);

    int diag = 0;
    int offdiag = 0;

    for (int j = 0; j < ncols_r; ++j)
    {
      if ((row_r[j] < start) || (row_r[j] >= end))
        offdiag++;
      else
        diag++;
    }

    for (int j = 0; j < ncols_i; ++j)
    {
      if ((row_i[j] < start) || (row_i[j] >= end))
        offdiag++;
      else
        diag++;
    }

    d_nnz[i - offset] = diag;
    o_nnz[i - offset] = offdiag;

    ierr = MatRestoreRow(H_real->mat(), i, &ncols_r, &row_r, PETSC_NULLPTR);

    if (matrix == 'H')
      ierr = MatRestoreRow(H_imag->mat(), i, &ncols_i, &row_i, PETSC_NULLPTR);
  }

  EigenSolver::preallocate_matrix(matrix, size_matrix, local_size, d_nnz, o_nnz);

  for (int i = start; i < end; ++i)
  {
    PetscInt ncols_r = 0, ncols_i = 0;
    const PetscInt *row_r, *row_i;
    PetscErrorCode ierr;
    const PetscReal *values_r, *values_i;

    ierr = MatGetRow(H_real->mat(), i, &ncols_r, &row_r, &values_r);

    if (matrix == 'H')
      ierr = MatGetRow(H_imag->mat(), i, &ncols_i, &row_i, &values_i);

    set<unsigned int> idset;
    for (int j = 0; j < ncols_r; ++j)
      idset.insert(row_r[j]);
    for (int j = 0; j < ncols_i; ++j)
      idset.insert(row_i[j]);

    vector<unsigned int> column_ids;
    column_ids.reserve(idset.size());
    for (set<unsigned int>::iterator it(idset.begin()); it != idset.end(); ++it)
      column_ids.push_back(*it);


    vector<Complex> complex_values;
    complex_values.reserve(column_ids.size());

    vector<unsigned int> nonzero_ids;
    nonzero_ids.reserve(column_ids.size());

    for (int j = 0, k = 0, l = 0; j < column_ids.size(); ++j)
    {
      Complex value = 0.0;
      while ((k < ncols_r) && (row_r[k] < column_ids[j]))
        ++k;

      if ((k < ncols_r) && (row_r[k] == column_ids[j]))
        value += values_r[k];

      while ((l < ncols_i) && (row_i[l] < column_ids[j]))
        ++l;

      if ((l < ncols_i) && (row_i[l] == column_ids[j]))
        value += Complex(0.0, values_i[l]);

      if (value != 0.0)
      {
        nonzero_ids.push_back(column_ids[j]);
        complex_values.push_back(value);
      }
    }


    EigenSolver::insert_matrix_row(matrix, i, nonzero_ids, complex_values);


    ierr = MatRestoreRow(H_real->mat(), i, &ncols_r, &row_r, &values_r);

    if (matrix == 'H')
      ierr = MatRestoreRow(H_imag->mat(), i, &ncols_i, &row_i, &values_i);
  }

  EigenSolver::finalize_matrix_assembly(matrix);
}


  
//========================================================================================//
void FEMEigenvalueProblem::do_copy_H_to_solver( )
{
  copy_matrix_to_solver('H');
}

//============================================================//

void FEMEigenvalueProblem::do_copy_S_to_solver()
{
  copy_matrix_to_solver('S');
}


//=======================================================================================/
void FEMEigenvalueProblem::apply_bc()
{

  if (solver_opt.Dirichlet_bc_everywhere)
  {
    apply_dirichlet_at_all_boundaries();
    make_new_dofs();
  }
  else
  {
    create_dirichlet_dofs();
    
    //make_constraints(); //creates a copy of them
    
    make_nodes_periodic();
    
    apply_periodic_bc();
    
    make_new_dofs();
  }
  
}
//=================================================================================
void FEMEigenvalueProblem::apply_periodic_bc()
{

  DofMap& dof_map = system->get_dof_map();
  //dof_map.create_dof_constraints(get_mesh());

  unsigned int number_of_variables = dof_map.n_variables();
  
  double tolerance = 1e-6;
 
  // The dimension that we are running
  //const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving


  vector<unsigned int> uvar(number_of_variables);

  for (unsigned int i = 0; i < number_of_variables; i++) 
  {
    uvar[i] = i;
  }

  unsigned int system_number=system->number();
  
 
  
  FEType fe_type = dof_map.variable_type(uvar[0]);
  
 
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
   

  

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

 
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_component;
  
  const double pos_tol = 1e-10;
  const double func_tol = 1e-10;
 
  //dof_map.print_dof_constraints();  

  for (int i = 0; i < get_mesh().mesh_dimension(); i++) //Loop over all the mesh directions
  {

    if  (solver_opt.periodicity[i]) //Check if the periodic b.c. are applied along the direction i
    {

      std::vector <const Node*>& vec =  nodes_periodic[i];

      for (unsigned int n = 0; n < vec.size(); n++) // Loop over all the nodes
      {
        const Node* node1 = vec[n];

        for (unsigned int var_index = 0 ; var_index  <  number_of_variables;  var_index ++)
        {//let us find dof for it-----------------


          // const Node& node = mesh.node(n);

          const unsigned int  n_dof = node1->dof_number(system_number,uvar[var_index],0);

          //dof is found-------------------------------



          //only if the dof is not constrained and not a dirichlet node do the job
          if (!(dof_map.is_constrained_dof(n_dof) || dirichlet_dofs.count(n_dof)))
          {
            //let us make a  point that lies at the opposite side
            Point point2(*node1);
            point2(i) = point2(i) + max_coord[i] - min_coord[i] - tolerance;


            //corresponding point is created


            //let us find an element this point belongs to and calculate the constraints

            //the most coarse element first
            unsigned int refinement_level = 0;
            MeshBase::const_element_iterator el3  = this->active_local_elements_begin();
            MeshBase::const_element_iterator end_el3 = this->active_local_elements_end();

            const Elem* elem1;
            bool found = false;

            unsigned int el_number = 0;

            for ( ; ( (el3 != end_el3) ) ; ++el3)
            {
              const Elem* elem = *el3;

              for (unsigned int ns = 0; ns < elem->n_sides(); ++ns)
              {
                ElementSide es(elem, ns);
                if (get_environment().is_boundary(es))
                {
                  if (elem->contains_point(point2))
                  {
                    elem1 = elem;
                    found = true;
                    break;

                  }
                }
                if (found) break;
              }
              el_number++;
            }

            //if (!found)  throw ModelErrorException("EnvelopFunctionApprox: Mesh periproblem");
            if (!found)
            {
              continue;
            }

            _constrained_dof_periodicity[n_dof](i) = min_coord[i] - max_coord[i];



            //active elem1 contains the opposite  point, we can constrain it now

            DofConstraintRow constraint;
            constraint.clear();


            dof_map.dof_indices (elem1, dof_indices_component, uvar[var_index]);

            std::vector<Point> point2_vec(1);

            point2_vec[0] = point2;

            std::vector<Point> point2_ref_vec(1);


            FEInterface::inverse_map (elem1->dim(), fe_type , elem1,  point2_vec,  point2_ref_vec)  ;

            fe->reinit (elem1, &point2_ref_vec);

            Point point_temp = point2_ref_vec[0];


            double sum = 0;
            for (int i1 = 0; i1 < phi.size(); i1++)
            {

              //if ( std::abs(phi[i1][0]) >  func_tol )
              {

                constraint[dof_indices_component[i1]] = phi[i1][0];
                sum += phi[i1][0];

              }
            }

            //constraint[n_dof] = sum;



            dof_map.add_constraint_row (n_dof,  constraint);
            my_dof_constraints.insert(make_pair(n_dof, constraint));

          }


        }
      }


    }
  }

  dof_map.process_constraints(get_mesh());
  system->reinit();
}
//---------------------------------------------------------------------------------------------

void FEMEigenvalueProblem::make_nodes_periodic()
{
  const double pos_tol = 1e-9;
 
  nodes_periodic.clear();

  unsigned int system_number = system->number();

  for (unsigned dir = 0; dir <= dim-1; dir++)
  {//directions
    std::vector<const Node*> temp_vec;
    
    if (solver_opt.periodicity[dir]) 
    {  
      // Loop over all the nodes
      for (unsigned int n = 0; n < get_mesh().n_nodes(); n++)
      {
	const Node* node1 = get_mesh().node_ptr(n);
	// look only at the ones used by the module
	if (node1->n_dofs(system_number) > 0)
	{		
	  if (std::abs((*node1)(dir) - min_coord[dir]) < pos_tol)
	    temp_vec.push_back(node1);
	}
      }
    }
    nodes_periodic.push_back(temp_vec);
  }
}

//--------------------------------------------------------------------------------------//
void FEMEigenvalueProblem::do_init()
{

}

//============================================================//

void FEMEigenvalueProblem::print_H(const std::string& outpath) const
{
  //std::string path = SimulationOptions::scratch_path;

  _H_real->print_matlab(outpath+"/Hr.m");
  _H_imag->print_matlab(outpath+"/Hi.m");
  if (_haveS) _S_real->print_matlab(outpath+"/Sr.m");

}


int FEMEigenvalueProblem::get_H_dim() const 
{ 
  return  _H_real->n();
}
    

int 
FEMEigenvalueProblem::get_H_nnz() const 
{ 
 

  DofMap& dof_map = system->get_dof_map();

  const std::vector<unsigned int>& n_nz = dof_map.get_n_nz();

  unsigned int row_start = _H_real->row_start();
  unsigned int row_stop = _H_real->row_stop();
  unsigned int nnz=0;

  for (unsigned int i = row_start ; i < row_stop; i++)
    nnz += n_nz[i];

  return nnz;

}



void 
FEMEigenvalueProblem::get_H_csr(std::vector<Complex>& A,
                                std::vector<int>& JA,
                                std::vector<int>& IA,
                                const std::vector<unsigned int>& perm) const

{

  const Scaling& sc = get_scaling();
  double scaling = get_H_units() * sc.get_length_scaling() / sc.get_calc_mesh_units();
  switch (get_mesh().mesh_dimension())
  {
    case 3:
      scaling *= sc.get_length_scaling() / sc.get_calc_mesh_units();
    case 2:
      scaling *= sc.get_length_scaling() / sc.get_calc_mesh_units();
    default:
      break;
  }

  get_csr(A, JA, IA, _H_real, _H_imag, scaling, perm);
}



void
FEMEigenvalueProblem::get_S_csr(std::vector<Complex>& A,
                                std::vector<int>& JA,
                                std::vector<int>& IA,
                                const std::vector<unsigned int>& perm) const
{

  // scale away length scaling in integration
  const Scaling& sc = get_scaling();
  double scale = sc.get_length_scaling() / sc.get_calc_mesh_units();
  switch (get_mesh().mesh_dimension())
  {
    case 3:
      scale *= sc.get_length_scaling() / sc.get_calc_mesh_units();
    case 2:
      scale *= sc.get_length_scaling() / sc.get_calc_mesh_units();
    default:
      break;
  }

  get_csr(A, JA, IA, _S_real, nullptr, scale, perm);
}



void
FEMEigenvalueProblem::get_csr(std::vector<libMesh::Complex>& A,
                              std::vector<int>& JA,
                              std::vector<int>& IA,
                              SparseMatrix<Real>* Mreal,
                              SparseMatrix<Real>* Mimag,
                              double scaling,
                              const std::vector<unsigned int>& perm) const
{

  if (Mreal == nullptr)
    throw RuntimeException("Real matrix part must not be "
        "nullptr in conversion to csr");

  PetscMatrix<Number>* Mr = static_cast<PetscMatrix<Number>* >(Mreal);
  Mreal->close();

  PetscMatrix<Number>* Mi = nullptr;
  if (Mimag != nullptr)
  {
    Mi = static_cast<PetscMatrix<Number>* >(Mimag);
    Mimag->close();
  }

  unsigned int row_start = Mr->row_start();
  unsigned int row_stop  = Mr->row_stop();
  unsigned int row, col, ind = 0;

  unsigned int n_rows = row_stop - row_start;

  std::vector<size_t> inv_perm(perm.size());
  if (!perm.empty())
  {
    if (n_rows != perm.size())
      throw std::runtime_error("In CSR matrix passing: "
          "permutation is incompatible with (local) matrix size.");

    for (size_t i = 0; i < perm.size(); ++i)
      inv_perm[perm[i]] = i;
  }

  IA.resize(n_rows + 1);

  IA[0] = 0;


  //DofMap& dof_map = system->get_dof_map();

  for (unsigned int i = 0; i < n_rows; ++i)
  {

    unsigned int row = perm.empty() ? i + row_start : perm[i];

    //if (dof_map.is_constrained_dof(_inv_perm[row]))
    //  continue;

    int ierr = 0;
    const PetscScalar *petsc_row_vals_real;
    const PetscScalar *petsc_row_vals_imag;
    const PetscInt *petsc_cols;
    int n_cols_real = 0;
    int n_cols_imag = 0;
    
    ierr = MatGetRow(Mr->mat(), row, &n_cols_real, &petsc_cols, &petsc_row_vals_real);
    CHKERRABORT(libMesh::GLOBAL_COMM_WORLD,ierr);

    if (Mi != nullptr)
    {
      ierr = MatGetRow(Mi->mat(), row, &n_cols_imag, &petsc_cols, &petsc_row_vals_imag);
      CHKERRABORT(libMesh::GLOBAL_COMM_WORLD,ierr);

      if (n_cols_real != n_cols_imag)
        Messages::error("n_cols_real != n_cols_imag");
    }


    std::vector<unsigned int> j_order(n_cols_real);
    std::iota(j_order.begin(), j_order.end(), 0);


    // reshuffle according to increasing sequence
    if (!perm.empty())
    {
       std::sort(j_order.begin(), j_order.end(), Compare(petsc_cols, perm));
    }

    for (unsigned int j = 0; j < n_cols_real; j++)
    {
      unsigned int jj = j_order[j];

      unsigned int col = petsc_cols[jj];
      if (!perm.empty())
        col = inv_perm[col];

      //if (dof_map.is_constrained_dof(_inv_perm[col]))
      //  continue;

      if (A.size() <= ind)
      {
        A.resize(5 * A.size() / 4);
        JA.resize(5 * A.size() / 4);
      }

      double val_r = petsc_row_vals_real[jj];
      double val_i = (Mi == nullptr) ? 0.0 : petsc_row_vals_imag[jj];

      A[ind] = scaling * Complex(val_r, val_i);
      JA[ind] = col;

      ind++;  
    }
   
    IA[i + 1]= ind;

    ierr = MatRestoreRow(Mr->mat(), row ,&n_cols_real, &petsc_cols, &petsc_row_vals_real);
    CHKERRABORT(libMesh::GLOBAL_COMM_WORLD,ierr);
      
    if (Mi != nullptr)
    {
      ierr = MatRestoreRow(Mi->mat(), row ,&n_cols_imag, &petsc_cols, &petsc_row_vals_imag);
      CHKERRABORT(libMesh::GLOBAL_COMM_WORLD,ierr);
    }
 

  }
  A.resize(ind);
  JA.resize(ind);

}


