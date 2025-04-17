// $Id$

#include "OpticsKP.h"
#include "SimulationInterface.h"
#include "KPbulkHamiltonian.h"
#include "EigenvalueProblem.h"
#include "EnvelopFunctionApprox.h"
#include "EFAbulkModel.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Messages.h"
#include "SimulationOptions.h"
#include "ModelOptions.h"

#include "libmesh/libmesh_common.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

#include <cstdlib>

#include "TiberModule.h"

using namespace std;
using namespace Constants;
using namespace libMesh;


OpticsKP::~OpticsKP()
{
}

 

//===============================================//
PhysicalModel* OpticsKP::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  ModelOptions kp8x8options = options;

  kp8x8options["model"] = "8x8";

  EFAbulkModel* model = PhysicalModel::create<EFAbulkModel>("EFAmodel", mat, kp8x8options);

  if (model == NULL)
    throw ModelErrorException("OpticsKP: cannot create EFAbulkModel");

  return(model);

}



//===============================================//
OpticsKP::OpticsKP(const ModelOptions& options)
 : Optics(options)
{
 
  has_solution_vector(false);
}



//==============================================//
void OpticsKP::do_init()
{

  Optics::do_init(); 

  _compute_intraband_terms = get_option("compute_intraband_terms", _compute_intraband_terms);

  // we need to remap some pointers in order to access some EFA stuff.
  initial_state_model = static_cast<EnvelopFunctionApprox*> (_initial_state_model); 
  final_state_model = static_cast<EnvelopFunctionApprox*> (_final_state_model); 

  //-------------------------------------------------------------------------------------------------
  // Take TiberCAD mesh 
  const MeshBase& mesh = get_mesh();

  unsigned int dim = mesh.mesh_dimension();
  double mesh_units = get_mesh_units();

  Scaling& scaling = get_scaling();

  scaling.set_length_scaling(Constants::bohr_radius);

  scaling.set_calc_mesh_units(mesh_units);

  //-------------------------------------------------------------------------------------------------
  //-Initialize a (fake) system used to build the P- matrix on FEM --------------------------
  //
  string system_name(get_equation_system_name());
  EquationSystems& es = get_equation_systems();
  es.add_system<libMesh::LinearImplicitSystem>(system_name);
  system = &(es.get_system<libMesh::LinearImplicitSystem>(system_name));

  //-------------------------------------------------------------------------------------------------
  //add variables for an 8x8 k.p problem

  psivar.resize(8);

  vector<string> psi_name;

  psi_name.clear();

  for (short i = 0; i < 8; i++)
  {
    std::ostringstream var_str;
    var_str << "psi" << i ;
    string name = var_str.str();
    psi_name.push_back(name);

    system->add_variable(name, FIRST, &(this->get_region_ids()));
  }


  for (unsigned int i = 0; i < 8; i++)
  {
    psivar[i] = system->variable_number(psi_name[i]);
  }

  //-----------------------------------------------------------------------------------------------
  //add matrixes

  system->add_matrix("Px_real");
  Px_matr_real = & (system->get_matrix("Px_real"));

  system->add_matrix("Py_real");
  Py_matr_real = & (system->get_matrix("Py_real"));

  system->add_matrix("Pz_real");
  Pz_matr_real = & (system->get_matrix("Pz_real"));

  system->add_matrix("Px_imag");
  Px_matr_imag = & (system->get_matrix("Px_imag"));

  system->add_matrix("Py_imag");
  Py_matr_imag = & (system->get_matrix("Py_imag"));

  system->add_matrix("Pz_imag");
  Pz_matr_imag = & (system->get_matrix("Pz_imag"));

  system->init();


}


//=====================================================================================================

void OpticsKP::do_compute_matrix_elements(void)
{
  assemble(); // this assembles Px Py Pz 

  unsigned int n_i =  _initial_state_numbers.size();
  unsigned int n_f =  _final_state_numbers.size();

  for (unsigned i = 0; i < 3; i++)
  {
    _P_matrix[i].resize(n_i);
    for (unsigned j = 0; j < n_i; j++)   _P_matrix[i][j].resize(n_f);
  }

  for (unsigned int i1 = 0; i1 < n_i; i1++)
    for (unsigned int i2 = 0; i2 < n_f; i2++)
    {
      unsigned int is = _initial_state_numbers[_initial_indices[i1]];
      unsigned int fs = _final_state_numbers[_final_indices[i2]];

      std::vector<libMesh::Complex> mat_el =  calculate_matrix_element(is, fs);
      for (unsigned i = 0; i < 3; i++)  _P_matrix[i][i1][i2] = mat_el[i];
    }

}
//=====================================================================================================

void OpticsKP::calculate_matrix_bulk(void)
{

  const Elem* mat_elem = initial_state_model->return_bulk_element();

  Point qp = mat_elem->vertex_average();


  EFAbulkHamiltonian* element_hamiltonian;
  element_hamiltonian = get_bulk_model<EFAbulkModel>(mat_elem)->get_Hamiltonian_model();

  element_hamiltonian->set_k_vector(_k_vector);


  KPbulkHamiltonian* element_kp_hamiltonian;
  element_kp_hamiltonian = dynamic_cast<KPbulkHamiltonian*>  (element_hamiltonian);

  const vector < vector <vector <EFAbulkHamiltonian::MatrixElement> > >&
    P = element_kp_hamiltonian->get_optical_operator() ;


  unsigned int n_i =  _initial_state_numbers.size();
  unsigned int n_f =  _final_state_numbers.size();

  for (unsigned i = 0; i < 3; i++)
  {
    _P_matrix[i].resize(n_i);
    for (unsigned j = 0; j < n_i; j++)
    {
      _P_matrix[i][j].clear();
      _P_matrix[i][j].resize(n_f, 0.0);
    }
  }



  //!number of bands in initial state
  //short    num_bands_initial = initial_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_in = initial_state_model->get_band_map();

  assert( kp_bands_map_in.size() > 0);

  //number of bands in final state
  //short    num_bands_final   = final_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_fi = final_state_model->get_band_map();

  assert( kp_bands_map_fi.size() > 0);

  map<short, short>::const_iterator  band_it;




  for (unsigned int band1 = 0; band1 < 8; band1++)
    for (unsigned int band2 = 0; band2 < 8; band2++)
    {

      band_it = kp_bands_map_in.find( band1);

      if (band_it != kp_bands_map_in.end())
      {

        short number1 = band_it->second;



        band_it = kp_bands_map_fi.find( band2);

        if (band_it != kp_bands_map_fi.end())
        {

          short number2 = band_it->second;

          for (unsigned int i1 = 0; i1 < n_i; i1++)
          {
            for (unsigned int i2 = 0; i2 < n_f; i2++)
            {
              unsigned int is = _initial_state_numbers[_initial_indices[i1]];
              unsigned int fs = _final_state_numbers[_final_indices[i2]];

              for (short pol = 0; pol < 3; pol++)
              {
                _P_matrix[pol][i1][i2] +=  P[pol][band1][band2].constant *
                    conj( _i_states[is].eigen_vector[number1] ) *
                    ( _f_states[fs].eigen_vector[number2] );

              }
            }
          }
        }
      }

    }
}

//=========================================================================//
void OpticsKP::do_assemble(const ModelOptions& opts)
{

  const MeshBase* mesh = &get_mesh();
  unsigned int dim = mesh->mesh_dimension();

  libMesh::DofMap& dof_map = system->get_dof_map();

  // this makes at least the windows version to crash:
  //system->reinit();

  libMesh::FEType fe_type = dof_map.variable_type(psivar[0]);

  std::unique_ptr<libMesh::FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  // A 5th order Gauss quadrature rule for numerical integration.
  libMesh::QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe -> attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW = fe->get_JxW();

  // properties at the quadrature points.
  const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();


  //------------------------------------------------------------
  std::vector<unsigned int> dof_indices_component;

  std::vector<unsigned int> dof_indices;

  //-------------------------------------------------------------
  Px_matr_real->zero();
  Py_matr_real->zero();
  Pz_matr_real->zero();
  Px_matr_imag->zero();
  Py_matr_imag->zero();
  Pz_matr_imag->zero();


  libMesh::DenseMatrix<Number> Px_real;
  libMesh::DenseMatrix<Number> Px_imag;

  libMesh::DenseMatrix<Number> Py_real;
  libMesh::DenseMatrix<Number> Py_imag;

  libMesh::DenseMatrix<Number> Pz_real;
  libMesh::DenseMatrix<Number> Pz_imag;


  libMesh::DenseSubMatrix<Number> Px_real_sub(Px_real);
  libMesh::DenseSubMatrix<Number> Px_imag_sub(Px_imag);

  libMesh::DenseSubMatrix<Number> Py_real_sub(Py_real);
  libMesh::DenseSubMatrix<Number> Py_imag_sub(Py_imag);

  libMesh::DenseSubMatrix<Number> Pz_real_sub(Pz_real);
  libMesh::DenseSubMatrix<Number> Pz_imag_sub(Pz_imag);


  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
  //                                                              //
  // we do not apply neither electric potential nor strain        //
  // because in our model they do not affect optical properties   //
  //                                                              //

  EFAbulkHamiltonian* element_hamiltonian;
  KPbulkHamiltonian* element_kp_hamiltonian;


  for ( ; el != end_el ; ++el)
  {//el
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem* elem = *el;

    element_hamiltonian = get_bulk_model<EFAbulkModel>(elem)->get_Hamiltonian_model();

    element_kp_hamiltonian = dynamic_cast<KPbulkHamiltonian*>  (element_hamiltonian);

    // ******************************************************* for k integration **********
 
    element_kp_hamiltonian->set_k_vector(_k_vector);

    // *********************************************************************************

 
    dof_map.dof_indices (elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();
    fe->reinit (elem);

    Px_real.resize(n_dofs, n_dofs );
    Px_imag.resize(n_dofs, n_dofs );

    Py_real.resize(n_dofs, n_dofs );
    Py_imag.resize(n_dofs, n_dofs );

    Pz_real.resize(n_dofs, n_dofs );
    Pz_imag.resize(n_dofs, n_dofs );



    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {//qp
      const vector < vector <vector <EFAbulkHamiltonian::MatrixElement> > >&
        P = element_kp_hamiltonian->get_optical_operator() ;



      for (unsigned int band1 = 0; band1 < 8; band1++)
      {//band1
        dof_map.dof_indices (elem, dof_indices_component, psivar[band1]);

        const unsigned int n_psi_dofs = dof_indices_component.size();

        for (unsigned int band2 = 0; band2 < 8; band2++)
        {//band2
          Px_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
          Px_imag_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);

          Py_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
          Py_imag_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);

          Pz_real_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);
          Pz_imag_sub.reposition(psivar[band1]*n_psi_dofs, psivar[band2]*n_psi_dofs, n_psi_dofs, n_psi_dofs);

          //Integration----------------------------------------------------------------
          for (unsigned int p1=0; p1<n_psi_dofs; p1++)
          {//p1
            for (unsigned int p2=0; p2<n_psi_dofs; p2++)
            {//p2


              vector<libMesh::Complex> value(3, libMesh::Complex(0.0,0.0));

              //----constant part-----------------------------------
              double temp = JxW[qp] * phi[p1][qp] * phi[p2][qp];

              for (short pol = 0; pol < 3; pol++)
                value[pol] += temp * P[pol][band1][band2].constant ;

              // intraband part
              if (_compute_intraband_terms && (band1 == band2))
              {
                libMesh::Complex val(0.0, -JxW[qp] * phi[p1][qp]);
                value[0] += val * dphi[p2][qp](0);
                value[1] += val * dphi[p2][qp](1);
                value[2] += val * dphi[p2][qp](2);
              }


              //----            ------------------------------------

              //--------linear left-----------------------------------------
              for (short i = 0; i < dim; i++)
              {
                temp = JxW[qp]* dphi[p1][qp](i) * phi[p2][qp];
                for (short pol = 0; pol < 3; pol++)
                  value[pol] -= temp * P[pol][band1][band2].linear_left[i]*libMesh::Complex(0.0, -1.0) ;
              }
              //--------           ----------------------------------------


              //--------linear right---------------------------------------
              for (short i = 0; i < dim; i++)
              {
                temp = JxW[qp]* dphi[p2][qp](i) * phi[p1][qp];
                for (short pol = 0; pol < 3; pol++)
                  value[pol] += temp * P[pol][band1][band2].linear_right[i] * libMesh::Complex(0.0, -1.0) ;
              }
              //-----------------------------------------------------------

              // for (short pol = 0; pol < 3; pol++) value[pol] *= my_Jacobian;

              Px_real_sub(p1,p2) += value[0].real();
              Py_real_sub(p1,p2) += value[1].real();
              Pz_real_sub(p1,p2) += value[2].real();

              Px_imag_sub(p1,p2) += value[0].imag();
              Py_imag_sub(p1,p2) += value[1].imag();
              Pz_imag_sub(p1,p2) += value[2].imag();


            }


          }


        }

      }

    }
    Px_matr_real->add_matrix(Px_real,dof_indices );
    Py_matr_real->add_matrix(Py_real,dof_indices );
    Pz_matr_real->add_matrix(Pz_real,dof_indices );

    Px_matr_imag->add_matrix(Px_imag,dof_indices );
    Py_matr_imag->add_matrix(Py_imag,dof_indices );
    Pz_matr_imag->add_matrix(Pz_imag,dof_indices );

  }
  Px_matr_imag->close();
  Px_matr_real->close();
  Py_matr_imag->close();
  Py_matr_real->close();
  Pz_matr_imag->close();
  Pz_matr_real->close();

  //Px_matr_imag->print_matlab("Px_i.m");
  //Px_matr_real->print_matlab("Px_r.m");

}


//----------------------------------------------------------------------------------------//

std::vector<libMesh::Complex> OpticsKP::calculate_matrix_element(unsigned int i, unsigned int j)
{
  std::vector<libMesh::Complex> result(3,libMesh::Complex(0.0,0.0));

  vector<  PetscMatrix<libMesh::Number>* > P_real_p(3);
  vector<  PetscMatrix<libMesh::Number>* > P_imag_p(3);

  const MeshBase& mesh = system->get_mesh();


  libMesh::DofMap& dof_map = system->get_dof_map();


  // ?!! compute nuber of nodes of the active elements: is it safe this way ??
  // n_variables should count the total number of bands ?  
  unsigned int number_of_nodes  = dof_map.n_dofs()/dof_map.n_variables(); 



  P_real_p[0]  = static_cast< PetscMatrix<libMesh::Number>* > (Px_matr_real); P_real_p[0]->close();
  P_real_p[1]  = static_cast< PetscMatrix<libMesh::Number>* > (Py_matr_real); P_real_p[1]->close();
  P_real_p[2]  = static_cast< PetscMatrix<libMesh::Number>* > (Pz_matr_real); P_real_p[2]->close();
  P_imag_p[0]  = static_cast< PetscMatrix<libMesh::Number>* > (Px_matr_imag); P_imag_p[0]->close();
  P_imag_p[1]  = static_cast< PetscMatrix<libMesh::Number>* > (Py_matr_imag); P_imag_p[1]->close();
  P_imag_p[2]  = static_cast< PetscMatrix<libMesh::Number>* > (Pz_matr_imag); P_imag_p[2]->close();

  const vector< libMesh::Complex >&  eigen_vector_i =  _i_states[i].eigen_vector;
  const vector< libMesh::Complex >&  eigen_vector_f =  _f_states[j].eigen_vector;


  int size_matrix = Px_matr_real->n();


  //number of bands in initial state
  short num_bands_initial = initial_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_in = initial_state_model->get_band_map();

  assert( kp_bands_map_in.size() > 0);

  //number of bands in final state
  short num_bands_final   = final_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_fi = final_state_model->get_band_map();

  assert( kp_bands_map_fi.size() > 0);

  map<short, short>::const_iterator  band_it;

  //AutoPtr<NumericVector<double>> eigvec_i_re = NumericVector<double>::build();
  //eigvec_i_re->init(size_matrix);

  // create vectors with size of the P-matrix and put the solution vectors'
  // entries into the right place
  numeric_index_type local_size = Px_matr_real->row_stop() - Px_matr_real->row_start();
  libMesh::PetscVector<Real> vec_i_r(this->get_solver_communicator(), size_matrix, local_size, PARALLEL);
  libMesh::PetscVector<Real> vec_i_i(this->get_solver_communicator(), size_matrix, local_size, PARALLEL);
  libMesh::PetscVector<Real> vec_f_r(this->get_solver_communicator(), size_matrix, local_size, PARALLEL);
  libMesh::PetscVector<Real> vec_f_i(this->get_solver_communicator(), size_matrix, local_size, PARALLEL);

  for (size_t i = 0; i < size_matrix; ++i)
  {
    // calculate the dof in the optical matrix
    // assumes band-major ordering
    //short band_number = i / number_of_nodes;
    // NOTE 15/3/2016 in new libmesh default is node major
    short band_number = i % 8;

    band_it = kp_bands_map_in.find(band_number);

    if (band_it != kp_bands_map_in.end())
    {
      //band exists in kp model of the initial state
      unsigned int dof_in_initial_eigenvector =
          band_it->second + (i / 8) * num_bands_initial;
          //band_it->second * number_of_nodes + i % number_of_nodes;

      libMesh::Complex value = eigen_vector_i[dof_in_initial_eigenvector];
      vec_i_r.set(i, std::real(value));
      vec_i_i.set(i, std::imag(value));
    }


    band_it = kp_bands_map_fi.find(band_number);

    if (band_it != kp_bands_map_fi.end())
    {
      //band exists in kp model of the initial state
      unsigned int dof_in_final_eigenvector =
          band_it->second + (i / 8) * num_bands_final;
          //band_it->second * number_of_nodes + i % number_of_nodes;

      libMesh::Complex value = eigen_vector_f[dof_in_final_eigenvector];
      vec_f_r.set(i, std::real(value));
      vec_f_i.set(i, std::imag(value));
    }
  }
  vec_i_r.close();
  vec_i_i.close();
  vec_f_r.close();
  vec_f_i.close();


  libMesh::PetscVector<Real> tmp_r(this->get_solver_communicator(), size_matrix, local_size, PARALLEL);
  libMesh::PetscVector<Real> tmp_i(this->get_solver_communicator(), size_matrix, local_size, PARALLEL);
  // now we have the initial and final eigenvector expanded to the bands
  // of the optical matrix

  for (short pol = 0; pol < 3; pol++)
  {
    // y_r = P_r * v_r - P_i * v_i
    // y_i = P_r * v_i + P_i * v_r

    tmp_r.zero();
    tmp_i.zero();
    P_imag_p[pol]->vector_mult(tmp_r, vec_i_i);
    tmp_r.scale(-1.0);
    P_real_p[pol]->vector_mult_add(tmp_r, vec_i_r);

    P_imag_p[pol]->vector_mult(tmp_i, vec_i_r);
    P_real_p[pol]->vector_mult_add(tmp_i, vec_i_i);

    double mme_r = vec_f_r.dot(tmp_r) + vec_f_i.dot(tmp_i);
    double mme_i = vec_f_r.dot(tmp_i) - vec_f_i.dot(tmp_r);

    result[pol] = libMesh::Complex(mme_r, mme_i);
  }

  return(result);
}





//----------------------------------------------------------------------------------------//


