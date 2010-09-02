// $Id$

#include "OpticsKP.h"
#include "SimulationInterface.h"
#include "KPbulkHamiltonian.h"
#include "EnvelopFunctionApprox.h"
#include "EFAbulkModel.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "DataOutput.h"
#include "SimulationOptions.h"

#include <stdlib.h>

using namespace std;
using namespace Constants;

Device*   OpticsKP:: _device;

OpticsKP::~OpticsKP()
{

  delete(_energy_mesh);

}



//===============================================//
PhysicalModel* OpticsKP::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{
  ModelOptions kp8x8options = options;

  kp8x8options["model"] = "kp";
  kp8x8options["kp_model"] = "8x8";

  EFAbulkModel* model = dynamic_cast<EFAbulkModel*>
    ( PhysicalModelInterface::create("EFAmodel", kp8x8options ) );

  if (model == NULL)
    throw ModelErrorException("OpticsKP: cannot create EFAbulkModel");

  return(model);

}

//===============================================//
BoundaryProperties* OpticsKP::create_boundary_model(const ModelOptions& options) const
 throw (ModelErrorException)
{

  return NULL;
}

//===============================================//
OpticsKP::OpticsKP(const ModelOptions& options)
 : Optics(options)
{
  initial_state_model = NULL;
  final_state_model = NULL;
  _energy_mesh = NULL;

}

//==============================================//

void OpticsKP::parse_options()
{
  const ModelOptions& mod_opt = get_options();

  _initial_eigen_state_numbers.clear();
  const std::vector<EnvelopFunctionApprox::eigen_propblem_solution>& in_solution 
                                               = initial_state_model->get_solution();
  {
    std::vector<unsigned int> temp;
    mod_opt.get_option("initial_eigenstates", temp);

    if (temp.size() == 2)
      if (temp[0] <= temp[1])
        if (temp[0] >= 0 && temp[1]-temp[0] <=in_solution.size())
        {
          _initial_eigen_state_numbers.resize(temp[1] - temp[0] + 1 );
          unsigned int j = 0;
          for (unsigned i = temp[0]; i <= temp[1]; i++)
          {
            _initial_eigen_state_numbers[j] = i;
            j++;
          }
        }

    if (_initial_eigen_state_numbers.size() == 0)
      throw InitFailedException("OpticsKP: initial states are not acceptable\n");
  }

  _final_eigen_state_numbers.clear();
  const std::vector<EnvelopFunctionApprox::eigen_propblem_solution>& fin_solution =
    final_state_model->get_solution();

  {
    std::vector<unsigned int> temp;
    mod_opt.get_option("final_eigenstates", temp);
    if (temp.size() == 2)
      if (temp[0] <= temp[1])
        if (temp[0] >= 0 && temp[1]-temp[0] <=fin_solution.size())
        {
          _final_eigen_state_numbers.resize(temp[1] - temp[0] + 1 );
          unsigned int j = 0;
          for (unsigned i = temp[0]; i <= temp[1]; i++)
          {
            _final_eigen_state_numbers[j] = i;
            j++;
          }
        }

    if (_final_eigen_state_numbers.size() == 0)
      throw InitFailedException("OpticsKP: final states are not acceptable\n");
  }




  //k-vector
  std::vector<double> k_vec(3, 0.0);
  mod_opt.get_option("k_vector",k_vec);
  if (k_vec.size() == 3)
  {

    for (short i = 0; i < 3; i++) k_vector[i] = k_vec[i];


  }
  else
    throw InitFailedException( "OpticsKP: k_vector size must be equal to 3 instead of " + k_vec.size());







  std::string  job_name = mod_opt.get_option("job","matrix_elements");

  if (job_name == "matrix_elements")
    job = MATREL;
  else if (job_name == "bulk_matrix_elements")
    job = BULKMATREL;
  else
    throw InitFailedException( "OpticsKP: Incorrect job: " + job_name);


}

//==============================================//
void OpticsKP::do_init()
{

  SimulationEnvironment& si = get_environment();

  _device = &( si.get_device() );



  //initial state----------------
  const ModelOptions& mod_opt = get_options();
  if  (mod_opt.find_option("initial_state_model"))
  {
    std::string quantum_model;
    quantum_model = mod_opt.get_option("initial_state_model" , "");
    initial_state_model = dynamic_cast<EnvelopFunctionApprox*> ( find_simulation ( quantum_model));
    if (initial_state_model == NULL)
      throw InitFailedException("OpticsKP: initial_state_model " + quantum_model + " does not exist\n");
  }
  else
    throw InitFailedException("OpticsKP: initial_state_model must be defined\n");
  //--------------------------

  //final state----------------
  if  (mod_opt.find_option("final_state_model"))
  {
    std::string quantum_model;
    quantum_model = mod_opt.get_option("final_state_model" , "");
    final_state_model = dynamic_cast<EnvelopFunctionApprox*> ( find_simulation ( quantum_model));
    if (final_state_model == NULL)
      throw InitFailedException("OpticsKP: final_state_model " + quantum_model + " does not exist\n");
  }
  else
    throw InitFailedException("OpticsKP: final_state_model must be defined\n");
  //---------------------------
  system_name = get_equation_system_name();
  es = &(get_equation_systems());

  //-------------------------------
  //-------------------------------------------------------------------------------------------------
  //My Jacobian
  double mesh_units = get_environment().get_device().get_mesh_units();

  Scaling& scaling = get_scaling();

  scaling.set_length_scaling(Constants::bohr_radius);

  scaling.set_calc_mesh_units(mesh_units);

  unsigned int dim = (es->get_mesh()).mesh_dimension();



  //--------------------------------------------------------------------------------------------------
  const MeshBase* mesh = &(es->get_mesh());


  es->add_system<LinearImplicitSystem> (system_name);

  system = &( es->get_system<LinearImplicitSystem>(system_name));


  //-------------------------------------------------------------------------------------------------
  //add variables

  psivar.resize(8);

  vector<string> psi_name;

  psi_name.clear();

  for (short i = 0; i < 8; i++)
  {
    std::ostringstream var_str;
    var_str << "psi" << i ;
    string name = var_str.str();
    psi_name.push_back(name);

    system->add_variable(name,FIRST);
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


//==============================================//






//==============================================//
void OpticsKP::do_solve()
{
  int verbose = SimulationOptions::verbose();

  if (verbose > 0)
    cout << "calculation of  matrix elements for dipole optical transition..." << flush;

  parse_options();


  if (job == BULKMATREL)
    calculate_matrix_bulk(); //only for bulk
  else
  {
    calculate_matrix(); //normal calculation
    calculate_P_matrix_elements();
  }

  unsigned int n1 =  _initial_eigen_state_numbers.size();
  unsigned int n2 =  _final_eigen_state_numbers.size();



  // calculate_P_matrix_elements();


  if (verbose > 0)
    cout << "done\n" << flush;



  if (verbose > 2)
  {
    for (int i = 0; i < n1; i++)
      for (int j = 0; j < n2; j++)
      {

        for (int p = 0; p < 3; p++)
        {
          cout << "polarization = " <<  p << "  " << "state i = " << i  <<"   " << "state j =   "
            << j <<"      "  << P_matrix[p][i][j] << "\n" << flush;


        }

      }
  }

}



//=====================================================================================================


void OpticsKP::calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz,
    std::map<const Elem*, double>& spectrum )
{


  spectrum.clear();




  std::vector<double> fs_eigen_values;
  std::vector<double> is_eigen_values;

  std::vector<double> fs_occupations;
  std::vector<double> is_occupations;




  double     trans_energy, f1, f2;


  unsigned int n1 =  _initial_eigen_state_numbers.size();
  unsigned int n2 =  _final_eigen_state_numbers.size();


  initial_state_model->get_eigenenergies(is_eigen_values);  

  final_state_model->get_eigenenergies(fs_eigen_values);


  initial_state_model->get_occupations(is_occupations);

  final_state_model->get_occupations(fs_occupations);



  // loop on  eigenstates

  for (unsigned i = 0; i < n1; i++)  // "upper" states
  {

    for (unsigned j = 0; j < n2; j++)  // "lower" states
    {

      trans_energy =  is_eigen_values[_initial_eigen_state_numbers[i]]
	             - fs_eigen_values[ _final_eigen_state_numbers[j]];


      f1 = is_occupations[_initial_eigen_state_numbers[i]];   // occupation for  electron

      f2 = fs_occupations[_final_eigen_state_numbers[j]]; // occupation for  holes





      Complex Me = P_matrix[0][i][j] * polariz(1) +
	           P_matrix[1][i][j] * polariz(2) +
	           P_matrix[2][i][j] * polariz(3);



      MeshBase::const_element_iterator       el     = Energy.active_elements_begin();
      const MeshBase::const_element_iterator end_el = Energy.active_elements_end();

      for ( ; el!= end_el ; ++el)
      {

        const Elem* elem = *el;

        double En = elem->centroid()(0);

        double Lorenzian =  0.5*Gamma/( ( trans_energy - En) *  ( trans_energy - En)
                            + (0.5*Gamma)*(0.5*Gamma)) * Hartree;
	// Note(alex): the division by Hartree seems wrong. Lorenzian is in 1/eV, so transformation
	// should be "Lorenzian * Hartree"

        double c = 1.0/Constants::fine_structure_constant;

        double omega = trans_energy/Hartree;

        //This is the right formula, as f1 is electron occupation probability and f2 is hole occupation probability.
        //Note that it differs from usual literature where usually f1 and f2 states initial state and final state
	//occupation probability, so it's related to electrons and it becomes f1*(1-f2)

        spectrum[elem] += 1 / (2 * M_PI ) * (omega * omega) /(c*c*c)  * Lorenzian * abs (Me) * abs (Me) * f1 * f2;

	//Note(alex): This factor 1/(2*PI*PI) was changed to 1/(2*PI). nr still missing
        //
	//According to Chuang's book the recombination rate should contain a pre-factor
	//(including 2 for spin sum and 2 for polarization and 4 Pi for angle integration)
	//
	//    nr^2 w^2           pi e^2      2       8 nr w e^2 c
	// ---------------- --------------- --- = ----------------------
	// pi^2 hbar c^2     nr c m^2 e0 w   V    4 pi e0 hbar (m^2 c^4)
	//
	// Extracting the prefactor m/hbar from the P-matrix, and multipling by a factor (hbar w) to get power emitted
	// we get the following prefactor (in which V has been removed to get total power emitted):
	//
	//     8 nr (hbar w)^2 e^2/(4 pi e0)
	//  = ------------------------------
	//             (hbar c)^3  hbar
	//
	//  Expressed in atomic units, hbar=1, e=1, 4 pi e0=1, c=1/fine_struct.
	//
	//  So the formula above, multiplied by 4*Pi for angle integration, and a factor of 4 for spin/pol deg
	//  agrees to Chuang's only if the factor is nr/(2*PI) rather than 1/(2*PI*PI).
	//

      }


    }


  }





}





//=====================================================================================================




//=====================================================================================================
void OpticsKP::calculate_P_matrix_elements ( )
{
  unsigned int n_i =  _initial_eigen_state_numbers.size();
  unsigned int n_f =  _final_eigen_state_numbers.size();


  P_matrix.clear();
  P_matrix.resize(3);
  for (unsigned i = 0; i < 3; i++)
  {
    P_matrix[i].resize(n_i);
    for (unsigned j = 0; j < n_i; j++)   P_matrix[i][j].resize(n_f);
  }

  for (unsigned int i1 = 0; i1 < n_i; i1++)
    for (unsigned int i2 = 0; i2 < n_f; i2++)
    {
      vector<Complex> matr_elements =  calculate_matrix_element(i1, i2);
      for (unsigned i = 0; i < 3; i++)  P_matrix[i][i1][i2] = matr_elements[i];
    }



}

void OpticsKP::calculate_matrix_bulk(void)
{

  const Elem* mat_elem = initial_state_model->return_bulk_element();

  Point qp = mat_elem->centroid();


  EFAbulkHamiltonian* element_hamiltonian;

  const ID subdomain = mat_elem->subdomain_id();

  const Material* mat = _device->get_material(subdomain);

  element_hamiltonian =
    (  dynamic_cast<EFAbulkModel*> (  mat ->get_model(get_id()) ))->get_Hamiltonian_model();

  element_hamiltonian->set_k_vector(k_vector);





  KPbulkHamiltonian* element_kp_hamiltonian;
  element_kp_hamiltonian = dynamic_cast<KPbulkHamiltonian*>  (element_hamiltonian);

  const vector < vector <vector <EFAbulkHamiltonian::MatrixElement> > >&
    P = element_kp_hamiltonian->get_optical_operator() ;







  unsigned int n_i =  _initial_eigen_state_numbers.size();
  unsigned int n_f =  _final_eigen_state_numbers.size();




  P_matrix.clear();
  P_matrix.resize(3);
  for (unsigned i = 0; i < 3; i++)
  {
    P_matrix[i].resize(n_i);
    for (unsigned j = 0; j < n_i; j++)   P_matrix[i][j].resize(n_f);
  }



  //!number of bands in initial state
  short    num_bands_initial = initial_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_in = initial_state_model->get_kp_bands();

  assert( kp_bands_map_in.size() > 0);

  //number of bands in final state
  short    num_bands_final   = final_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_fi = final_state_model->get_kp_bands();

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
	    for (unsigned int i2 = 0; i2 < n_f; i2++)
	    {

	      for (short pol = 0; pol < 3; pol++)
	      {
		P_matrix[pol][i1][i2] +=  P[pol][band1][band2].constant *
		  conj( (initial_state_model->get_solution())[i1].eigen_vector[number1] )*
		  ( (final_state_model->get_solution())[i2].eigen_vector[number2] );



	      }
	    }
	}
      }

    }











}

//=========================================================================//
void OpticsKP::calculate_matrix(void)
{

  unsigned int dim = (es->get_mesh()).mesh_dimension();


  //--------------------------------------------------------------------------------------------------------//
  const MeshBase* mesh = &(es->get_mesh());


  DofMap& dof_map = system->get_dof_map();



  Px_matr_real->zero();


  Py_matr_real->zero();


  Pz_matr_real->zero();


  Px_matr_imag->zero();


  Py_matr_imag->zero();


  Pz_matr_imag-> zero();

  system->reinit();




  FEType fe_type = dof_map.variable_type(psivar[0]); //all the variable have the same FE representation



  AutoPtr<FEBase> fe (  build_finite_element(dim, fe_type, true)  );

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe -> attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW = fe->get_JxW();

  // properties at the quadrature points.
  const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();


  // The element shape function gradients evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();


  //------------------------------------------------------------
  std::vector<unsigned int> dof_indices_component;

  std::vector<unsigned int> dof_indices;

  //-------------------------------------------------------------

  DenseMatrix<Number> Px_real;
  DenseMatrix<Number> Px_imag;

  DenseMatrix<Number> Py_real;
  DenseMatrix<Number> Py_imag;

  DenseMatrix<Number> Pz_real;
  DenseMatrix<Number> Pz_imag;


  DenseSubMatrix<Number> Px_real_sub(Px_real);
  DenseSubMatrix<Number> Px_imag_sub(Px_imag);

  DenseSubMatrix<Number> Py_real_sub(Py_real);
  DenseSubMatrix<Number> Py_imag_sub(Py_imag);

  DenseSubMatrix<Number> Pz_real_sub(Pz_real);
  DenseSubMatrix<Number> Pz_imag_sub(Pz_imag);


  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();
  //                                                              //
  // we do not apply neither electric potential nor strain        //
  // because in our model they do not affect optical properties   //
  //                                                              //

  EFAbulkHamiltonian* element_hamiltonian;
  KPbulkHamiltonian* element_kp_hamiltonian;

  unsigned int el_number = 0;

  double temp;
  for ( ; el != end_el ; ++el)
  {//el
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem* elem = *el;

    const ID subdomain = elem->subdomain_id();
    const Material* mat = _device->get_material(subdomain);

    element_hamiltonian = (  dynamic_cast<EFAbulkModel*> (  mat ->get_model(get_id()) )  )->get_Hamiltonian_model();

    element_kp_hamiltonian = dynamic_cast<KPbulkHamiltonian*>  (element_hamiltonian);

    // ******************************************************* for k integration **********
    element_kp_hamiltonian->set_k_vector(k_vector);

    // *********************************************************************************


    dof_map.dof_indices (elem, dof_indices);
    const unsigned int n_dofs   = dof_indices.size();
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


              vector<Complex> value(3, Complex(0.0,0.0));

              //----constant part-----------------------------------
              temp = JxW[qp] * phi[p1][qp] * phi[p2][qp];

              for (short pol = 0; pol < 3; pol++)
                value[pol] += temp * P[pol][band1][band2].constant ;



              //----            ------------------------------------

              //--------linear left-----------------------------------------
              for (short i = 0; i < dim; i++)
              {
                temp = JxW[qp]* dphi[p1][qp](i) * phi[p2][qp];
                for (short pol = 0; pol < 3; pol++)
                  value[pol] -= temp * P[pol][band1][band2].linear_left[i]*Complex(0.0, -1.0) ;
              }
              //--------           ----------------------------------------


              //--------linear right---------------------------------------
              for (short i = 0; i < dim; i++)
              {
                temp = JxW[qp]* dphi[p2][qp](i) * phi[p1][qp];
                for (short pol = 0; pol < 3; pol++)
                  value[pol] += temp * P[pol][band1][band2].linear_right[i] * Complex(0.0, -1.0) ;
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


    el_number++;

  }

}


//----------------------------------------------------------------------------------------//

std::vector<Complex> OpticsKP::calculate_matrix_element(unsigned int i, unsigned int j)
{
  vector<Complex> result(3,Complex(0.0, 0.0));

  vector<  PetscMatrix<Number>* > P_real_p(3);
  vector<  PetscMatrix<Number>* > P_imag_p(3);

  const MeshBase& mesh = system->get_mesh();





  DofMap& dof_map = system->get_dof_map();

  unsigned int number_of_nodes  = dof_map.n_dofs()/dof_map.n_variables(); //only nodes that belong to active elements



  P_real_p[0]  = static_cast< PetscMatrix<Number>* > (Px_matr_real); P_real_p[0]->close();
  P_real_p[1]  = static_cast< PetscMatrix<Number>* > (Py_matr_real); P_real_p[1]->close();
  P_real_p[2]  = static_cast< PetscMatrix<Number>* > (Pz_matr_real); P_real_p[2]->close();
  P_imag_p[0]  = static_cast< PetscMatrix<Number>* > (Px_matr_imag); P_imag_p[0]->close();
  P_imag_p[1]  = static_cast< PetscMatrix<Number>* > (Py_matr_imag); P_imag_p[1]->close();
  P_imag_p[2]  = static_cast< PetscMatrix<Number>* > (Pz_matr_imag); P_imag_p[2]->close();

  const vector< Complex >   eigen_vector_i =  (initial_state_model->get_solution())[i].eigen_vector;
  const vector< Complex >   eigen_vector_f =  (final_state_model->get_solution())[j].eigen_vector;



  int size_matrix = Px_matr_real->n();



  //!number of bands in initial state
  short    num_bands_initial = initial_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_in = initial_state_model->get_kp_bands();

  assert( kp_bands_map_in.size() > 0);

  //number of bands in final state
  short    num_bands_final   = final_state_model->get_number_of_bands();
  const map<short, short>&  kp_bands_map_fi = final_state_model->get_kp_bands();

  assert( kp_bands_map_fi.size() > 0);

  map<short, short>::const_iterator  band_it;



  for (short pol = 0; pol < 3; pol++)
  {//polarization
    for (int row = 0 ; row < size_matrix; row++)
    {//rows of P matrix

      // short band_number1 = row%8;
      short band_number1 = row/number_of_nodes;

      band_it = kp_bands_map_in.find( band_number1 );


      if (band_it != kp_bands_map_in.end())
      {//band exists in kp model of the initial state

        unsigned int dof_in_initial_eigenvector = band_it->second * number_of_nodes + row%number_of_nodes;
        int ierr = 0;
        const  PetscScalar *petsc_row_vals_real;
        const  PetscInt *petsc_cols_real;
        int n_cols_real = 0;

        const  PetscScalar *petsc_row_vals_imag;
        const  PetscInt *petsc_cols_imag;
        int n_cols_imag = 0;


        map<int,double> real_values, imag_values;
        map<int, double>::iterator  position;


        set<int> real_column, imag_column, complex_column;
        set<int>::iterator com_col_it;

        insert_iterator<set<int> >  com_ins(complex_column,complex_column.begin() );


        ierr = MatGetRow(P_real_p[pol]->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
        CHKERRABORT(libMesh::COMM_WORLD,ierr);
        real_column.clear();

        for (int i = 0; i < n_cols_real; i++)
          if (petsc_row_vals_real[i] != 0.0)
          {
            real_column.insert(petsc_cols_real[i]);

            real_values.insert(make_pair(petsc_cols_real[i],petsc_row_vals_real[i] ));
          }


        ierr = MatGetRow(P_imag_p[pol]->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
        CHKERRABORT(libMesh::COMM_WORLD,ierr);
        imag_column.clear();


        for (int i = 0; i < n_cols_real; i++)
          if (petsc_row_vals_imag[i] != 0.0)
          {
            imag_column.insert(petsc_cols_imag[i]);
            imag_values.insert(make_pair(petsc_cols_imag[i],petsc_row_vals_imag[i] ));
          }


        set_union(real_column.begin(), real_column.end(), imag_column.begin(), imag_column.end(), com_ins);

        for (com_col_it = complex_column.begin(); com_col_it != complex_column.end(); com_col_it++)
        {
          int n1 = *com_col_it;

          short band_number2  = n1/number_of_nodes;

          band_it = kp_bands_map_fi.find(band_number2 );

          if (band_it != kp_bands_map_fi.end())
          {//band exists in kp model of the final state
            double value_real;
            double value_imag;

            //real part------
            position = real_values.find(n1);
            if (position != real_values.end())
              value_real = position->second;
            else
              value_real = 0.0;

            //imag part
            position = imag_values.find(n1);
            if (position != imag_values.end())
              value_imag = position->second;
            else
              value_imag = 0.0;

            Complex value_complex = Complex(value_real, value_imag);


            unsigned int dof_in_final_eigenvector = band_it->second * number_of_nodes + n1%number_of_nodes;


            result[pol] += value_complex * conj(eigen_vector_i[dof_in_initial_eigenvector]) * eigen_vector_f[dof_in_final_eigenvector];


          }
        }


        ierr = MatRestoreRow(P_real_p[pol]->mat(), row ,&n_cols_real, &petsc_cols_real,&petsc_row_vals_real);
        CHKERRABORT(libMesh::COMM_WORLD,ierr);

        ierr = MatRestoreRow(P_imag_p[pol]->mat(), row ,&n_cols_imag, &petsc_cols_imag,&petsc_row_vals_imag);
        CHKERRABORT(libMesh::COMM_WORLD,ierr);
      }
    }
  }

  return(result);
}





//----------------------------------------------------------------------------------------//



void OpticsKP::do_plot()
{
  // get  spectrum calculation  options from  opticsKP model  section
  // for  calculation of  spectrum for a single k-point

  const ModelOptions& mod_spectrum = get_options();

  double  Gamma = mod_spectrum.get_option("broadening", 0.007);

  double   Emin,  Emax,dE;




  if (mod_spectrum.find_option("Emin"))
    Emin = mod_spectrum.get_option("Emin", 0.0);
  else
    throw InitFailedException("OpticsKP: Emin must be defined\n");



  if (mod_spectrum.find_option("Emax"))
    Emax = mod_spectrum.get_option("Emax", 0.0);
  else
    throw InitFailedException("OpticsKP: Emin must be defined\n");



  if (Emax < Emin)  throw InitFailedException("OpticsKP: Emax < Emin");

  if (mod_spectrum.find_option("dE"))
    dE = mod_spectrum.get_option("dE", 0.0);
  else
    throw InitFailedException("OpticsKP: dE must be defined\n");


  if (dE <= 0)  throw InitFailedException("OpticsKP: dE <= 0");

  unsigned int num_nodes = (int)((Emax - Emin)/dE) + 1;


  // do  energy_mesh
  _energy_mesh = new Mesh(1);



  MeshTools::Generation::build_cube (*_energy_mesh,
      num_nodes, 0, 0,
      Emin, Emax,
      0, 0,
      0, 0,
      EDGE2);


  //  calculate_spectrum             ************  only  for case  k_0
  //calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz,
  //                         std::map<const Elem*, double>& spectrum )

  //std::map<const Elem*, double> spectrum;

  std::map<const Elem*, double> spectrum_x;
  std::map<const Elem*, double> spectrum_y;
  std::map<const Elem*, double> spectrum_z;

  Tensor1 polariz_x(0); polariz_x(1) = 1.0;
  Tensor1 polariz_y(0); polariz_y(2) = 1.0;
  Tensor1 polariz_z(0); polariz_z(3) = 1.0;






  // calculate_spectrum( *_energy_mesh, Gamma,polariz, spectrum ) ;

  calculate_spectrum( *_energy_mesh, Gamma, polariz_x, spectrum_x ) ;
  calculate_spectrum( *_energy_mesh, Gamma, polariz_y, spectrum_y ) ;
  calculate_spectrum( *_energy_mesh, Gamma, polariz_z, spectrum_z ) ;





  //  Emin max , Gamma, polariz,  =  given in  optics model ???

  string dimension;
  double area_dim_factor = 1;



  vector<string> names(3);

  names[0] = "power_density_k0_Px[W/eV]";  //As it's a power density (W) per photon (eV)
  names[1] = "power_density_k0_Py[W/eV]";  //Depending on simulation dimension it will be
  names[2] = "power_density_k0_Pz[W/eV]";  // W/(eV) (3D), W/(eV*cm) (2D), W/(eV*cm^2) (1D)

  vector<double> results;

  {
    MeshBase::const_element_iterator       elem_it  = _energy_mesh->active_elements_begin();
    const MeshBase::const_element_iterator elem_end = _energy_mesh->active_elements_end();
    int n = 0;
    for(;elem_it != elem_end; ++elem_it)
      n++;

    results.resize(n*3);
  }


  MeshBase::const_element_iterator       elem_it  = _energy_mesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end = _energy_mesh->active_elements_end();
  int point = 0;

  for(;elem_it != elem_end; ++elem_it)
  {
    const Elem* el = *elem_it;
    double value;


    //--x - polarization
    value = spectrum_x[el];
    value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
    value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
    value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
    results[3*point + 0] = value;

    //--y - polarization
    value = spectrum_y[el];
    value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
    value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
    value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
    results[3*point + 1] = value;

    //--z - polarization
    value = spectrum_z[el];
    value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
    value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
    value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
    results[3*point + 2] = value;



    point++;
  }


  if (plot_solution("optical_spectrum_k_0"))
  {
    string filename(get_name() +
        "_spectrum_k_0" + TiberCad::get_filename_suffix());

    string format = get_options().get_option("output_format", "grace");

    DataOutput data_output(*_energy_mesh, format);
    data_output.set_output_directory(get_output_directory());

    data_output.write_cell_data(filename, results, names);
  }
}
