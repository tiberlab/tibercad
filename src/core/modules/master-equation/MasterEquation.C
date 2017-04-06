
// modules includes
#include "MasterEquation.h"
#include "Material.h"
#include "Constants.h"
#include "TiberNonlinearSystem.h"
#include "TiberLinearSystem.h"
#include "SolveFailedException.h"
#include "Variable.h"



// libmesh includes
#include "node.h"
#include "mesh.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_interface.h"
#include "equation_systems.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "numeric_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "petsc_matrix.h"



#include "DataOutput.h"
#include "Messages.h"



// C++ includes
#include <tr1/random>
#include <fstream>
#include <vector>



#include "TiberModule.h"
#include "Database.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "SimulationInterface.h"
#include "SimulationOptions.h"


// namespaces
using namespace std;
using namespace libMesh;
using namespace MasterEquationDefs;



// pointer _this initialize
MasterEquation*
MasterEquation::_this = NULL;



MasterEquation::Options::Options(void)
  : coupling(BOTH),
    solver_method(NEWTON)
{
}



MasterEquation::MasterEquation(const ModelOptions& options)
  : SimulationInterface(options),
    _rebuild_eq_system(true)
    //_reference_potential(0.0);
{
  // there's nothing to be done
}



MasterEquation::~MasterEquation(void)
{
  cleanup_solver();
}



const Database&
 MasterEquation::get_database(void)
 {
   return _owner->get_database();
 }



void
 MasterEquation::read_database(void)
 {
    //const Database& db = get_database();
    //db.set_section("valenceband");
    //double prova;
    //db.get("E_v_1", prova, true);

 }





void
MasterEquation::do_init(void)
{
  cerr<<"do_init_begin"<<endl;
  
  parse_options();
  
  //! last added
  _env = &get_environment();
  _device = &(_env->get_device());
  
  

  // Create non linear equation system
  create_equation_system("nonlinear");



  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();
  


  // add variables "Fermi leves" for electrons and holes
  system.add_variable("Ef_e", libMeshEnums::CONSTANT, libMeshEnums::MONOMIAL);
  system.add_variable("Ef_h", libMeshEnums::CONSTANT, libMeshEnums::MONOMIAL);



  // assembly of the system
  system.attach_assembly_routine(assemble);



  // attach transformation
  //system.attach_transformation_routine(transformation);



  /// initialize the system
  system.init();

  //cerr<<"do_init attach assembly done"<<endl;


  // Solutions vector (non-linear system)
  NumericVector<Number>& solution = get_solution_vector();



  // define number of variable which belong to the system
  const unsigned int Ef_e_var = system.variable_number("Ef_e");
  const unsigned int Ef_h_var = system.variable_number("Ef_h");



  // linear system of support to set values of energies (Ec and Ev) perturbated by gaussians
  create_equation_system("linear");
  TiberLinearSystem* _sys_EcEv = &get_equation_system<TiberLinearSystem>(1);



  // set band edges to initial values
  _sys_EcEv->add_variable("LUMO",CONSTANT,MONOMIAL);
  _sys_EcEv->add_variable("HOMO",CONSTANT,MONOMIAL);



  // initialize
  _sys_EcEv->init();



  // get dofs map for Ec and Ev
  const libMesh::DofMap& dof_map_E = _sys_EcEv->get_dof_map();



  //Solutions vector (linear system)
  NumericVector<Number>& solution_E = _sys_EcEv->get_solution_vector();



  // define number of variable which belong to the system
  const unsigned int LUMO_var = _sys_EcEv->variable_number("LUMO");
  const unsigned int HOMO_var = _sys_EcEv->variable_number("HOMO");

/*
  read_database();


  const Database& db = get_database();
  db.set_section("valenceband");
  //double prova;
  double data;
  db.get("E_v_1", data, true);

  //prova = data;

  //cerr << prova << endl;
  */



  // get mesh
  MeshBase& mesh = get_mesh();


  const int n_el = solution_E.size()/2;


  // define iterators which point to the active elements
  MeshBase::element_iterator el = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();


  //cerr<<"do_init set level done"<<endl;

  int seed = 110101010;//= std::tr1::random_device(); //time(NULL) *
  //const double log_sigma = log(1 / sigma * sqrt(2 * M_PI));
  std::tr1::mt19937 generator(seed);

  double kbT = SimulationOptions::temperature * Constants::kb; 

  const double dE_min = log(1e-4) * kbT;

  //cerr<<"do_init loop over elements"<<endl;

  for (; el != end; ++el)
  {
    const Elem* elem = *el;


    unsigned const id_Ef_e = elem->dof_number(system.number(), Ef_e_var, 0 );
    unsigned const id_Ef_h = elem->dof_number(system.number(), Ef_h_var, 0 );

    unsigned const id_LUMO = elem->dof_number(_sys_EcEv->number(), LUMO_var, 0 );
    unsigned const id_HOMO = elem->dof_number(_sys_EcEv->number(), HOMO_var, 0 );


    const double eps = std::numeric_limits<double>::min();

    double r_n;
    double rn_1;
    double rn_2;
    double rn_3;
    double rn_4;

    double r_n_c;
    double rn_1_c;
    double rn_2_c;
    double rn_3_c;
    double rn_4_c;


    for (unsigned int i = 0; i < 4; i++)
    {
      double ru1;
      do
      {
        ru1 = static_cast<double>(generator()) / generator.max();
      } while (ru1 <= eps);

      double ru2 = static_cast<double>(generator()) / generator.max();

      r_n = sigma * sqrt(-2.0 * log(ru1)) * cos(2.0 * M_PI * ru2);

      if (i == 0)
      {
        rn_1 = r_n;
      }
      if (i == 1)
      {
        rn_2 = r_n;
      }
      if (i == 2)
      {
        rn_3 = r_n;
      }
      else
      {
        rn_4 = r_n;
      }

    }


    for (unsigned int i = 0; i < 4; i++)
    {
      double ru1_c;
      do
      {
        ru1_c = static_cast<double>(generator()) / generator.max();
      } while (ru1_c <= eps);

      double ru2_c = static_cast<double>(generator()) / generator.max();

      r_n_c = sigma * sqrt(-2.0 * log(ru1_c)) * cos(2.0 * M_PI * ru2_c);

      if (i == 0)
      {
        rn_1_c = r_n_c;
      }
      if (i == 1)
      {
        rn_2_c = r_n_c;
      }
      if (i == 2)
      {
        rn_3_c = r_n_c;
      }
      else
      {
        rn_4_c = r_n_c;
      }

    }


    // Set energetic gaussian disorder for HOMO and LUMO
    solution_E.set(id_LUMO, level_LUMO + rn_1);
    solution_E.set(id_HOMO, level_HOMO + rn_2); //


    // Set initial guess for Fermi level
    solution.set(id_Ef_e, Ef); //
    solution.set(id_Ef_h, Ef); //


  }

  //cerr<<"do_init loop over elements done"<<endl;

  //cerr << solution << endl;
  //cerr << solution_E << endl;

  ofstream outFile;
  string outputpath = TiberCad::get_output_dir();
  outFile.open(outputpath+"/Energies_1.txt");
  outFile << solution_E << endl;
  outFile.close();

  ofstream outFile1;
  outFile1.open(outputpath+"/Densities.txt");
  outFile1 << solution << endl;
  outFile1.close();

  // cerr<<"do_init linear"<<endl;

  //cerr<<"do_init linear init"<<endl;

  //cerr <<"do_init: done"<<endl;
}


void
MasterEquation::do_newton(void)
{
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);

  system.set_options(get_solver_options());
  system.solve();

}


void
MasterEquation::reset_solver(void)
{
  if (!_rebuild_eq_system)
  {
    //get_equation_systems().delete_system(get_equation_system_name());
   _rebuild_eq_system = true;
  }
}


void
MasterEquation::cleanup_solver(void)
{
  reset_solver();
}



void
MasterEquation::do_solve(void)
{

  // this is dirty, but at the moment we have to provide a static method
  // for assembly. This should change in future
  _this = this;
  //cerr<<"do_solve"<<endl;


  //TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);
  //get_solution_vector().close();
  //system.get_vector("x") = get_solution_vector();

  //int coupling = get_my_options().coupling;

  Messages::info("Solving Newton");

  do_newton();

  Messages::info("Newton done");

  //system.set_options(get_solver_options());
  //cerr<<"solve"<<endl;
  //system.solve();
  //cerr<<"done"<<endl;

}



void
MasterEquation::rebuild_equation_system(void)
{
  if (!_rebuild_eq_system) return;


  //EquationSystems& equation_systems = get_equation_systems();
  clear_systems();


  ModelOptions& solveropts = get_solver_options();
  if (solveropts.get_option("absolute_tolerance", -1.0) < 0)
    solveropts["absolute_tolerance"] = "1e-3";


  // the coupled DD system
  create_equation_system("nonlinear");
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);

  system.attach_assembly_routine(assemble);


  system.add_variable("Ef_e", CONSTANT, MONOMIAL);
  system.add_variable("Ef_h", CONSTANT, MONOMIAL);

  //system.add_vector("me_old_sol");


  // finally initialize the newly created system
  system.init();


  _rebuild_eq_system = false;

}




void
MasterEquation::do_print_info(void)
{
  Messages::info("Module implementation example: simple Master Equations");

  Options& myopts = get_my_options();

  Messages::newline();

  ostringstream os;
  os << "solving for : ";
  if (myopts.coupling & ELECTRONS)
    os << "electrons";
  if (myopts.coupling & HOLES)
    os << "holes";
  if (myopts.coupling & BOTH)
    os << " electrons and holes";

  os << endl;

}



void
MasterEquation::parse_options(void)
{

  sigma = get_option("sigma",0.0);

  level_LUMO = get_option("LUMO",0.0); 

  level_HOMO = get_option("HOMO",0.0); 

  Ef = get_option("Ef", (level_HOMO+level_LUMO)/2.0 );

  alpha = get_option("alpha",1e-9);

  v_0 = get_option("v0",1e12); 
 
  string coupling = get_option("coupling", "");

  Options& myopts = get_options();
  if (coupling == "electrons")
    myopts.coupling = ELECTRONS;
  else if (coupling == "holes")
    myopts.coupling = HOLES;
  else
    myopts.coupling = BOTH;


}



// Coupling

void
MasterEquation::assemble(const libMesh::NumericVector<Number>& x,
    libMesh::NumericVector<Number>* residual,
      libMesh::SparseMatrix<Number>* jacobian,
        libMesh::NonlinearImplicitSystem&)
{
  //cerr << "assemble" << endl;

  switch (_this->_options.coupling)
  {
    case (ELECTRONS):
      _this->do_assembly<ELECTRONS>(x, residual, jacobian);
      break;
    case (HOLES):
      _this->do_assembly<HOLES>(x, residual, jacobian);
      break;
    default:
      _this->do_assembly<BOTH>(x, residual, jacobian);
      break;
  }

}

void
MasterEquation::transformation(NumericVector<Number>& u,
     NumericVector<Number>& T,
       NumericVector<Number>& TX, bool transf,
        libMesh::NonlinearImplicitSystem&)
{
  //cerr << "transformation" << endl;

  _this->do_transformation(u, T, TX, transf);
}



template <int coupling>
void
MasterEquation::do_assembly(const libMesh::NumericVector<Number>& x, libMesh::NumericVector<Number>* residual,
    libMesh::SparseMatrix<Number>* jacobian)
{

  //cerr << "do_assembly: init" << endl;
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  libMesh::NumericVector<Number>& solution = get_solution_vector();

  PetscMatrix<Number>* petsc_matrix = dynamic_cast<PetscMatrix<Number>* > (system.matrix);
  MatSetOption(petsc_matrix->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  unsigned int it_do_ass = 0;

  it_do_ass +=1;


  const unsigned int Ef_e_var = system.variable_number("Ef_e");
  const unsigned int Ef_h_var = system.variable_number("Ef_h");


  libMesh::DofMap& dof_map =  system.get_dof_map();


  TiberLinearSystem* _sys_EcEv = &get_equation_system<TiberLinearSystem>(1);


  const unsigned int LUMO_var = _sys_EcEv->variable_number("LUMO");
  const unsigned int HOMO_var = _sys_EcEv->variable_number("HOMO");

  libMesh::DofMap& dof_map_E = _sys_EcEv->get_dof_map();

  libMesh::NumericVector<Number>& solution_E = _sys_EcEv->get_solution_vector();

  //NumericVector<Number>& oldx = system.get_vector("me_old_sol");


  libMesh::FEType fe_type = system.variable_type(Ef_e_var); //

  // the volume finite element
  libMesh::UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));


  std::cout<<"sigma: "<<sigma<<endl;
  std::cout<<"Ef: "<<Ef<<endl;
  std::cout<<"lumo: "<<level_LUMO<<endl;


  // Jacobian * quadrature weight at each integration point.
  //const vector<Real>& JxW = fe->get_JxW();

  // element shape functions
  //const vector<vector<Real> >& phi = fe->get_phi();

  // element shape function gradients
  //const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  //cerr<<"do_assembly get centroids"<<endl;

  // physical coordinates of the quadrature points
  //const vector<Point>& center = fe->centroid();
  const vector<Point>& center = fe->get_xyz();


  // the system matrix (will hold also element jacobian contribution)
  libMesh::DenseMatrix<Number> Ke;
  // the system rhs (will hold also element rhs contribution)
  libMesh::DenseVector<Number> Fe;
  // the local solution
  libMesh::DenseVector<Number> X;
  
  
  // Dense matrix for BC
  //libMesh::DenseMatrix<Number> BC;


  libMesh::DenseSubMatrix<Number>
    Knn(Ke), Knp(Ke),
    Kpn(Ke), Kpp(Ke);

  libMesh::DenseSubVector<Number>
    Fn(Fe),
    Fp(Fe);

  libMesh::DenseSubVector<Number>
    Xn(X),
    Xp(X);

 // cerr<<"do_assembly matrix and vector init done"<<endl;

  vector<unsigned int> dof_indices_tot;
  vector<unsigned int> dof_indices_tot_neigh;

  vector<unsigned int> dof_indices_n;
  vector<unsigned int> dof_indices_p;

  vector<unsigned int> dof_indices_tot_E;
  vector<unsigned int> dof_indices_Ec;
  vector<unsigned int> dof_indices_Ev;
  //vector<unsigned int> dof_indices_Ev;



  // zero out residual and jacobian !! IMPORTANT !!
  if (residual != NULL)
    residual->zero();
  if (jacobian != NULL)
    jacobian->zero();



  MeshBase::const_element_iterator el =
                                    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end =
                                    mesh.active_local_elements_end();


  //cerr<<"do_assembly loop over elements begin"<<endl;

  int unsigned iter = 0;
  int unsigned iteration ;


  const int n_el = solution_E.size();

  const double n_elem = n_el/2;

  iteration = iteration + 1;

  double kbT = SimulationOptions::temperature * Constants::kb; 

  ofstream outFile_log;
  string outputpath = TiberCad::get_output_dir();
  outFile_log.open(outputpath+"/log.txt");

  //cerr<<"do_assembly: begin loop over elements"<<endl;

  // loop over all active elements
  for (; el != end; ++el)
  {

    //cerr << "element loop" << endl;

    const Elem* elem = *el;
    //cout << el << endl;
    iter += 1;

    //cerr<<"do_assembly bug1"<<endl;

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_tot);
    dof_map.dof_indices(elem, dof_indices_n, Ef_e_var);
    dof_map.dof_indices(elem, dof_indices_p, Ef_h_var);

    dof_map_E.dof_indices(elem, dof_indices_Ec, LUMO_var);
    dof_map_E.dof_indices(elem, dof_indices_Ev, HOMO_var);


    unsigned int n_dofs     = dof_indices_n.size(); // number of dofs of element for single variable
    unsigned int n_dofs_tot = dof_indices_tot.size();   // number of dofs of element for all variables


    //cerr << n_dofs << endl;
    //cerr << n_dofs_tot<< endl;


    //dof_indices_tot.reserve(30);
    //dof_indices_tot_neigh.reserve(30);
    //dof_indices_n.reserve(30);
    //dof_indices_p.reserve(30);
    //dof_indices_Ec.reserve(30);
    //dof_indices_Ev.reserve(30);



    //cerr << dof_indices_p[0] << endl;
    //cerr << dof_indices_tot[1] << endl;

    //cerr<<"do_assembly: add neighbors dof indices begin"<<endl;

    for (unsigned int k=0; k < elem->n_sides(); ++k)
    {



      const Elem* elem_neig = elem->neighbor(k); // pointer who points to neighbors

      if (elem_neig != NULL) // out of boundary
      {

        std::vector<unsigned int> neig_dof_indices_tot;
        std::vector<unsigned int> neig_dof_indices_n;
        std::vector<unsigned int> neig_dof_indices_p;

        //std::vector<unsigned int> neig_dof_indices_Ec;
        //std::vector<unsigned int> neig_dof_indices_Ev;

        //dof_map.dof_indices(elem_neig, neig_dof_indices_tot);
        dof_map.dof_indices(elem_neig, neig_dof_indices_n, Ef_e_var);
        dof_map.dof_indices(elem_neig, neig_dof_indices_p, Ef_h_var);

        //dof_map_E.dof_indices(elem_neig, neig_dof_indices_Ec, Ec_var);
        //dof_map_E.dof_indices(elem_neig, neig_dof_indices_Ev, Ev_var);
        //dof_map_E.dof_indices(elem_neig, neig_dof_indices_Ec);
        //dof_map_E.dof_indices(elem_neig, neig_dof_indices_Ev);

        //cerr << "neig_dof_indices_p = " << neig_dof_indices_p[0] << endl;
        //cerr << "neig_dof_indices_tot = " << neig_dof_indices_tot[1] << endl;


        //cerr<<"do_assembly bug7"<<endl;
        //dof_indices_tot.insert( dof_indices_tot.end(), neig_dof_indices_tot[0] );
        //dof_indices_n.insert( dof_indices_n.end(), neig_dof_indices_n[0] );
        //dof_indices_p.insert( dof_indices_p.end(), neig_dof_indices_p[0] );
        //dof_indices_E.insert( dof_indices_E.end(), neig_dof_indices_E[0] );

        //dof_indices_tot.insert( dof_indices_tot.end(), neig_dof_indices_tot.begin(), neig_dof_indices_tot.end() );
        //dof_indices_n.insert( dof_indices_n.end(), neig_dof_indices_n.begin(), neig_dof_indices_n.end() );
        //dof_indices_p.insert( dof_indices_p.end(), neig_dof_indices_p.begin(), neig_dof_indices_p.end() );
        //dof_indices_E.insert( dof_indices_E.end(), neig_dof_indices_E.begin(), neig_dof_indices_E.end() );

        // push_back

          //cerr << neig_dof_indices_tot[0] << endl;
          //cerr << neig_dof_indices_tot[1] << endl;
          //cerr << neig_dof_indices_tot[2] << endl;

          //dof_indices_tot.push_back(neig_dof_indices_tot[0]);
          //dof_indices_tot.push_back(neig_dof_indices_tot[1]);// add dofs of neighbors soon after the dofs of element

          dof_indices_n.push_back(neig_dof_indices_n[0]);
          dof_indices_p.push_back(neig_dof_indices_p[0]);

          //dof_indices_Ec.push_back(neig_dof_indices_n[0]);
          //dof_indices_Ev.push_back(neig_dof_indices_p[0]);
          //dof_indices_Ec.push_back(neig_dof_indices_Ec[0]);
          //dof_indices_Ev.push_back(neig_dof_indices_Ev[0]);
        //}
      }
    }




    //std::vector<unsigned int> dof_indices_tot_neigh;
    dof_indices_tot_neigh.clear();
    dof_indices_tot_neigh.insert(dof_indices_tot_neigh.begin(), dof_indices_n.begin(), dof_indices_n.end());
    dof_indices_tot_neigh.insert(dof_indices_tot_neigh.end(), dof_indices_p.begin(), dof_indices_p.end());

    dof_indices_Ec = dof_indices_n;
    dof_indices_Ev = dof_indices_p;

    //dof_indices_tot_E.clear();
    //dof_indices_tot_E.insert(dof_indices_tot_E.begin(), dof_indices_n.begin(), dof_indices_n.end());
    //dof_indices_tot_E.insert(dof_indices_tot_E.end(), dof_indices_p.begin(), dof_indices_p.end());

    //dof_indices_tot_neigh.insert(dof_indices_tot_neigh.end(), dof_indices_p.begin()+1, dof_indices_p.end());


    //cerr<<"do_assembly: add neighbors dof indices done"<<endl;

    unsigned int n_dofs_neigh     = dof_indices_n.size(); // number of dofs of element+neighbors for single variable
    unsigned int n_dofs_neigh_tot = dof_indices_tot_neigh.size();  // number of dofs of element+neighbors for all variables

    unsigned int n_dofs_neigh_E = dof_indices_Ec.size();  // number of dofs of element+neighbors for all variables

   //cerr<< "single variable: number dofs indices " <<endl;
    //cerr<< n_dofs_neigh <<endl;

    //cerr<< "both variables: number dofs indices " <<endl;
    //cerr<< n_dofs_neigh_tot <<endl;
    //for (unsigned int i =0; i < n_dofs_neigh_tot; i++ )
      //cerr << "dof_indices_tot = " <<  dof_indices_tot_neigh[i]<< endl;
    //for (unsigned int i =0; i < n_dofs_neigh; i++ )
    //  cerr << "dof_indices_n = " << dof_indices_n[i]<< endl;
    //for (unsigned int i =0; i < n_dofs_neigh; i++ )
    //  cerr << "dof_indices_p = " <<  dof_indices_p[i]<< endl;
    //for (unsigned int i =0; i < n_dofs_neigh; i++ )
    //  cerr << "dof_indices_Ec = " <<  dof_indices_Ec[i]<< endl;
    //for (unsigned int i =0; i < n_dofs_neigh; i++ )
    //  cerr << "dof_indices_Ev = " <<  dof_indices_Ev[i]<< endl;

    //unsigned int n_dofs_neigh_tot_max = elem->n_sides();



    Ke.resize(n_dofs_tot, n_dofs_neigh_tot);

    Fe.resize(n_dofs_tot); // n_dofs_tot

    //cerr << "bug Fe" << endl;

    X.resize(n_dofs_neigh_tot);
	
	//BC.resize(1, n_el-1);

    dof_map.extract_local_vector(x, dof_indices_tot_neigh, X);


    // Reposition the submatrices according to this scheme:
    //
    //         -           -          -  -
    //   Ke = | Knn Knp |   Fe = | Fn |
    //        | Kpn Kpp |;       | Fp |
    //

    Knn.reposition(0, 0, n_dofs, n_dofs_neigh);
    Knp.reposition(0, n_dofs_neigh, n_dofs, n_dofs_neigh); // nulli n_dofs_neigh


    //
    Kpn.reposition(n_dofs, 0, n_dofs, n_dofs_neigh); // nulli
    Kpp.reposition(n_dofs, n_dofs_neigh, n_dofs, n_dofs_neigh);

    //
    Fn.reposition(0, n_dofs);
    Fp.reposition(n_dofs, n_dofs);


    //
    Xn.reposition(0, n_dofs_neigh); // n_dofs_neigh
    Xp.reposition(n_dofs, n_dofs_neigh ); // n_dofs_neigh

    Knn.zero();
    Knp.zero();
    Kpn.zero();
    Kpp.zero();

    Fn.zero();
    Fp.zero();

    Ke.zero();

    //cerr << "Knn" << Knn << endl;

    //cerr << "Knp" << Knp << endl;

    //cerr << "Kpn" << Kpn << endl;

    //cerr << "Kpp" << Kpp << endl;

    //cerr << n_dofs_neigh << endl;


    //cerr << "iteration =" << iteration << endl;



    //const int sides = elem->n_sides();


    //cerr << "n_avg" << n_avg << endl;
    //cerr << "p_avg" << p_avg << endl;

    //
    // Each sub-matrix Ke_ij looks like this:
    //
    //      Ke_ij = summation on j != i of W_ji
    //
    // Each sub-vector Fe_i looks like this:
    //
    //      Fe_i = X_i * (summation on j != i of W_ij)
    //
    // The residual looks like this:
    //
    //      r_i = Ke_ij*X_j - Fe_i
    //
    // The jacobian looks like this:
    //
    //      J_ij = dr_i/dX_j
    //
    //           = Ke_ij - dFe_i/dX_j
    //



    double LUMO_i = solution_E(dof_indices_Ec[0]); // Ec of element
    double HOMO_i = solution_E(dof_indices_Ev[0]); // Ev of element

    //cerr << "LUMO_i = " << LUMO_i << endl;
    //cerr << "HOMO_i = " << HOMO_i << endl;
    //cerr << "Ef_e_i = " << Ef_e_i << endl;
    //cerr << "Ef_h_i = " << Ef_h_i << endl;
    //cerr << "Pe_i = " << Pe_i << endl;
    //cerr << "Ph_i = " << Ph_i << endl;




    // First we will build the system matrix Ke_ij
    if (jacobian != NULL)
    {

      //cerr << "prova" << endl;

      double Ef_e_i = solution(dof_indices_n[0]); // Ec of element
      double Ef_h_i = solution(dof_indices_p[0]); // Ev of element

      //cerr << "Ef_e_i = " << solution(dof_indices_n[0]) << endl;
      //cerr << "Ef_h_i = " << solution(dof_indices_p[0]) << endl;


      double Pe_i = 1/( 1 + exp( -( solution(dof_indices_n[0]) - solution_E(dof_indices_Ec[0]) )/kbT ) );
      double Ph_i = 1/( 1 + exp( -( solution_E(dof_indices_Ev[0]) - solution(dof_indices_p[0]) )/kbT ) );

      //cerr << "Pe_i = " << Pe_i << endl;
      //cerr << "Ph_i = " << Ph_i << endl;

      double LUMO_i = solution_E(dof_indices_Ec[0]); // Ec of element
      double HOMO_i = solution_E(dof_indices_Ev[0]); // Ev of element


      //cerr << "LUMO_i = " << PLUMO_i << endl;
      //cerr << "HOMO_i = " << HOMO_i << endl;


      for (unsigned int i = 0; i < n_dofs ; i++) // loop over i-elements
      {

        unsigned int ind = 0;  // index initialization

        // hopping rate calculation
        for (unsigned int j = 0; j < elem->n_sides(); j++)
        {
          const Elem* elem_neig = elem->neighbor(j); // pointer who points to neighbors

          if (elem_neig != NULL) // neighbors
          {

            ind += 1;

            Point d = elem->centroid() - elem_neig->centroid();
            double dist = d.size();
            double abs_dist = abs(dist);
            double exp_d = exp(-2 * alpha * abs_dist);

            // get conduction band edge for the calculation of Miller-Abrhams energetic term
            double LUMO_j = solution_E(dof_indices_Ec[ind]);

            // get valence band edge for the calculation of Miller-Abrhams energetic term
            double HOMO_j = solution_E(dof_indices_Ev[ind]);


            double Pe_j = 1/( 1 + exp( -( solution(dof_indices_n[ind]) - solution_E(dof_indices_Ec[ind]) )/kbT ) );

            double Ph_j = 1/( 1 + exp( -( solution_E(dof_indices_Ev[ind]) - solution(dof_indices_p[ind]) )/kbT ) );

            //cerr << "Pe_j = " << Pe_j << endl;

            // from site j to site i
            double delta_LUMO_ji = LUMO_j - LUMO_i;
            double delta_HOMO_ji = HOMO_j - HOMO_i;

            // from site i to site j
            double delta_LUMO_ij = LUMO_i - LUMO_j;
            double delta_HOMO_ij = HOMO_i - HOMO_j;

            // from site j to site i
            double exp_LUMO_ji = 0.0;
            double exp_HOMO_ji = 0.0;

            // from site i to site j
            double exp_LUMO_ij = 0.0;
            double exp_HOMO_ij = 0.0;


            //double d_pi_d_fi = -kbT * Ph_i * (1 - Ph_i);
            //double d_ni_d_fi = kbT * Pe_i * (1 - Pe_i);

            //double d_pj_d_fj = -kbT * Ph_j * (1 - Ph_j);
            //double d_nj_d_fj = kbT * Pe_j * (1 - Pe_j);



            double d_pi_d_fi = -( exp( -( solution(dof_indices_n[0]) - solution_E(dof_indices_Ec[0]) )/kbT) )/( kbT * ( 1 + exp( -( solution(dof_indices_n[0]) - solution_E(dof_indices_Ec[0]) )/kbT ) ) * ( 1 + exp( -( solution(dof_indices_n[0]) - solution_E(dof_indices_Ec[0]) )/kbT ) ));
            double d_ni_d_fi = exp( -( solution_E(dof_indices_Ev[0]) - solution(dof_indices_p[0]) )/kbT )/ ( kbT * ( 1 + exp( -( solution_E(dof_indices_Ev[0]) - solution(dof_indices_p[0]) )/kbT ) ) * ( 1 + exp( -( solution_E(dof_indices_Ev[0]) - solution(dof_indices_p[0]) )/kbT ) ));

            double d_pj_d_fj = -( exp( -( solution(dof_indices_n[ind]) - solution_E(dof_indices_Ec[ind]) )/kbT) )/( kbT * ( 1 + exp( -( solution(dof_indices_n[ind]) - solution_E(dof_indices_Ec[ind]) )/kbT ) ) * ( 1 + exp( -( solution(dof_indices_n[ind]) - solution_E(dof_indices_Ec[ind]) )/kbT ) ));
            double d_nj_d_fj = exp( -( solution_E(dof_indices_Ev[ind]) - solution(dof_indices_p[ind]) )/kbT  )/ ( kbT * ( 1 + exp( -( solution_E(dof_indices_Ev[ind]) - solution(dof_indices_p[ind]) )/kbT ) ) * ( 1 + exp( -( solution_E(dof_indices_Ev[ind]) - solution(dof_indices_p[ind]) )/kbT ) ) );


            if (delta_LUMO_ji > 0)
            {
              exp_LUMO_ji = 1;

              exp_LUMO_ij = exp( delta_LUMO_ij / kbT);
            }
            if (delta_LUMO_ji < 0)
            {
              exp_LUMO_ji = exp( -delta_LUMO_ij / kbT);

              exp_LUMO_ij = 1;
            }
            if (delta_LUMO_ji == 0)
            {
              exp_LUMO_ji = 1;

              exp_LUMO_ij = 1;
            }
            if (delta_HOMO_ji < 0)
            {
              exp_HOMO_ji = 1;

              exp_HOMO_ij = exp( -delta_HOMO_ij / kbT);
            }
            if (delta_HOMO_ji > 0)
            {
              exp_HOMO_ji = exp( delta_HOMO_ij / kbT);

              exp_HOMO_ij = 1;
            }
            if (delta_HOMO_ji == 0)
            {
              exp_HOMO_ji = 1;

              exp_HOMO_ij = 1;
            }

            // HO INVERTITO i E j DELLE DERIVATE E SEMBRA CONVERGERE
            if (coupling &BOTH)
            {
				
				//if (i == n_dofs - 1)
					
					//{
						//Knn(i,i) = v_0/n_elem*1 * d_nj_d_fj ;//*1* d_ni_d_fi;
						
						//Knn(i,ind) = v_0/n_elem*1 * d_ni_d_fi;
						
						//Kpp(i,i) = v_0/n_elem*1 * d_pj_d_fj;//*1* d_pi_d_fi;
						
						//Kpp(i,ind) = v_0/n_elem * d_pi_d_fi;
					//}
					
				//if (i < n_dofs - 1)
					
					//{
						Knn(i,i) +=  v_0 * exp_d * ( Pe_j * (exp_LUMO_ij - exp_LUMO_ji) - exp_LUMO_ij ) * d_nj_d_fj;  // diagonal  Jp - Jw

						Knn(i,ind) += v_0 * exp_d * ( exp_LUMO_ji - Pe_i * ( exp_LUMO_ij - exp_LUMO_ji ) ) * d_ni_d_fi; //off diagonal -> Wji - Pi * (  Wij - Wji )

						Kpp(i,i) +=  v_0 * exp_d * ( Ph_j * (exp_HOMO_ij - exp_HOMO_ji) - exp_HOMO_ij ) * d_pj_d_fj;  // diagonal Jp - Jw

						Kpp(i,ind) += v_0 * exp_d * ( exp_HOMO_ji - Ph_i * ( exp_HOMO_ij - exp_HOMO_ji ) ) * d_pi_d_fi; //off diagonal -> Wji - Pi * (  Wij - Wji )
					//}

            }

          }
        }
      }


      //cerr << "VECTOR x" << x << endl;
      //cerr << "solution" << solution << endl;
    }

    //cerr << "jacobian" << endl;

    if (jacobian != NULL)
    {
      //cerr << "Knn" << endl << Knn << endl;
      //cerr << "iter =" << iter << endl;
      //cerr << "Kpp" << endl << Kpp << endl;
    }

    int iter = 0;

    outFile_log << "Newton iteration " << iter << endl;



    //cerr << "test" << endl;



    if (residual != NULL)
    {
      //cerr << "test_1" << endl;

      int iter = iter +1;

      //cerr << "test = " << solution(dof_indices_n[0]) << endl;

      for (unsigned int i = 0; i < n_dofs ; i++) // trovare funzione per calc num di dofs di un elem
      {


        double Pe_i = 1/( 1 + exp( -( x(dof_indices_n[0]) - solution_E(dof_indices_Ec[0]) )/kbT  ) );
        double Ph_i = 1/( 1 + exp( -( solution_E(dof_indices_Ev[0]) - x(dof_indices_p[0]) )/kbT  ) );

        //outFile_log << "Pe_i_it = " << Pe_i << endl;
        //outFile_log << "Ph_i_it = " << Ph_i << endl;

        unsigned int ind = 0;

        for (unsigned int j = 0; j < elem->n_sides(); j++)
        {
          const Elem* elem_neig = elem->neighbor(j); // pointer who points to neighbors

          if (elem_neig != NULL)
          {
            ind += 1;

            Point d = elem->centroid() - elem_neig->centroid();
            double dist = d.size();
            double abs_dist = abs(dist);
            double exp_d = exp(-2 * alpha * abs_dist);

            // get conduction band edge for the calculation of Miller-Abrhams energetic term
            double LUMO_j = solution_E(dof_indices_Ec[ind]);

            // get valence band edge for the calculation of Miller-Abrhams energetic term
            double HOMO_j = solution_E(dof_indices_Ev[ind]);


            double Pe_j = 1/( 1 + exp( -( x(dof_indices_n[ind]) - solution_E(dof_indices_Ec[ind]) )/kbT  ) );

            double Ph_j = 1/( 1 + exp( -( solution_E(dof_indices_Ev[ind]) - x(dof_indices_p[ind]) )/kbT  ) );


            outFile_log << "Pe_j_it = " << Pe_i << endl;
            outFile_log << "Ph_j_it = " << Ph_j << endl;


            // from site j to site i
            double delta_LUMO_ji = LUMO_j - LUMO_i;
            double delta_HOMO_ji = HOMO_j - HOMO_i;

            // from site i to site j
            double delta_LUMO_ij = LUMO_i - LUMO_j;
            double delta_HOMO_ij = HOMO_i - HOMO_j;

            // from site j to site i
            double exp_LUMO_ji;
            double exp_HOMO_ji;

            // from site i to site j
            double exp_LUMO_ij;
            double exp_HOMO_ij;


            if (delta_LUMO_ji > 0)
            {
              exp_LUMO_ji = 1;

              exp_LUMO_ij = exp( delta_LUMO_ij / kbT);
            }
            if (delta_LUMO_ji < 0)
            {
              exp_LUMO_ji = exp( -delta_LUMO_ij / kbT);

              exp_LUMO_ij = 1;
            }
            if (delta_LUMO_ji == 0)
            {
              exp_LUMO_ji = 1;

              exp_LUMO_ij = 1;
            }
            if (delta_HOMO_ji < 0)
            {
              exp_HOMO_ji = 1;

              exp_HOMO_ij = exp( -delta_HOMO_ij / kbT);
            }
            if (delta_HOMO_ji > 0)
            {
              exp_HOMO_ji = exp( delta_HOMO_ij / kbT);

              exp_HOMO_ij = 1;
            }
            if (delta_HOMO_ji == 0)
            {
              exp_HOMO_ji = 1;

              exp_HOMO_ij = 1;
            }


            if (coupling & BOTH)
            {
              Fn(i) += v_0 * exp_d * Pe_j * ( exp_LUMO_ji + Pe_i * ( exp_LUMO_ij - exp_LUMO_ji ) ); // + Kij * Pj = Pj * ( Wji + Pi * ( Wij - Wji ) )

              Fn(i) -= v_0 * exp_d * exp_LUMO_ij * Pe_i; // - Fi = - Pi * sum (Wij)

              Fp(i) += v_0 * exp_d * Ph_j * ( exp_HOMO_ji + Ph_i * ( exp_HOMO_ij - exp_HOMO_ji ) ); // Kij * Pj = Wji + Pi * ( Wij - Wji )

              Fp(i) -= v_0 * exp_d * exp_HOMO_ij * Ph_i; // - Fi = - Pi * sum (Wij)


            }

          }

        }

      }


      residual->add_vector(Fe, dof_indices_tot);

      //cerr << "Fn = " << Fn << endl;
      //cerr << "Fp = " << Fp << endl;

    }
    else
    {
      jacobian->add_matrix(Ke, dof_indices_tot, dof_indices_tot_neigh);
    }

	//cerr<<"do_assembly_first: done"<<endl;
	
  }// end loop over all active elements

cerr<<"do_assembly: done"<<endl;

// boundary conditions jacobian
if (jacobian != NULL)
{
  const int n_el = solution_E.size();
  
  const double n_elem = n_el/2;
  
  libMesh::DenseMatrix<Number> BC;
  
  BC.resize(1, n_el);
 
  std::vector<unsigned int> BC_cols(n_el);
  std::vector<unsigned int> BC_rows_n(1);
  std::vector<unsigned int> BC_rows_p(1);
  

  BC_rows_n[0] = n_elem-1;
  BC_rows_p[0] = n_el-1;
 
 
  for (unsigned int l = 0; l < n_el  ; l++)
  {
     //cerr << 1/ n_elem << endl;
 	
     //jacobian->add(0, l , 1 / n_elem);
     //jacobian->set(n_elem - 1, l , 1 / n_elem);
     //jacobian->set(n_el - 1, l + n_elem, 1 / n_elem);
     //jacobian->set(n_elem - 1, l , 0);
     //jacobian->set(n_el - 1, l + n_elem, 0);
 	
 	   //jacobian->add(l ,n_elem - 1,  1 / n_elem);
 	
     BC_cols[l] = l;
   	 BC(0,l) = 1.0;
   }

   jacobian->add_matrix(BC, BC_rows_n, BC_cols);
   jacobian->add_matrix(BC, BC_rows_p, BC_cols);
 
    //jacobian->add_matrix(BC, BC_rows_p, BC_cols);
 	  //jacobian->set(n_elem - 1, n_elem-1,  v_0/n_elem);
    //jacobian->set(n_el - 1, n_el-1, v_0/n_elem);
}


// boundary conditions residual
if (residual != NULL)
{
  const int n_el = solution_E.size();

  const double n_elem = n_el/2;
  
  const double c = 1; 

  double BCn = 0.0;
  double BCp = 0.0;


  for (unsigned int j = 0; j < n_elem; j++)
  {

    BCn += 1/( 1 + exp( -( x(j) - solution_E(j) )/kbT  ) );
    BCp += 1/( 1 + exp( -(  solution_E(n_elem + j) - x(n_elem + j) )/kbT  ) );


    //BCn -= 1/ ( 1 + exp ( - ( Ef - solution_E(j) )/kbT ) );
    //BCp -= 1/ ( 1 + exp ( - ( solution_E(n_elem + j ) - Ef )/kbT ) );



  }

  //BCn -= 1/ ( 1 + exp ( - ( Ef - level_LUMO )/kbT ) );
  //BCp -= 1/ ( 1 + exp ( - ( level_HOMO - Ef )/kbT ) );

  BCn *= 1/n_elem;
  BCp *= 1/n_elem;
  
  
  BCn -= c;
  BCp -= c;
  
  //cerr <<  "BCn = " << BCn << endl;
  //cerr <<  "BCp = " << BCp << endl;
  
  //BCn -= 1;
  //BCp -= 1;
  
  //cerr <<  "BCn = " << BCn << endl;
  //cerr <<  "BCp = " << BCp << endl;


  //residual->set(0, BCn_0);
  residual->set(n_elem - 1, 0.0);
  residual->set(n_el - 1, 0.0);

}

if (residual != NULL)
{
  //cerr << "Fn = " << Fn << endl;
  //cerr << "Fp = " << Fp << endl;
}

    //cerr << "VECTOR X" << endl << X << endl;
    //cerr << "MATRIX Ke" << endl << Ke << endl;
    //cerr << "VECTOR Fe" << endl << Fe << endl;

  //cerr <<  "BCn_tot = " << BCn_0 << endl;
  //cerr <<  "BCp_tot = " << BCp_0 << endl;

   //cerr << "VECTOR x" << x << endl;

   if (jacobian != NULL)
   {
     //jacobian->print_matlab("jacobian.m");
   }
   if (residual != NULL)
   {
     //residual->print_matlab("residual.m");
   }


    if (jacobian != NULL)
    {
      jacobian->close();
    }
    else
    {
      residual->close();
    }




  string num;
  string it_el; // string which will contain the result

  ostringstream convert;   // stream used for the conversion
  ostringstream convert1;

  convert << it_do_ass;
  convert1 << iter;

  num = convert.str();
  it_el = convert1.str();


  //string outputpath = TiberCad::get_output_dir();


  ofstream outFile_x;
  outFile_x.open(outputpath+"/solution.txt");
  outFile_x << x << endl;
  outFile_x.close();


  ofstream outFile_res;
  outFile_res.open(outputpath+"/Jac_"+ num + ".txt");
  outFile_res << Ke << endl;
  outFile_res.close();


  outFile_log.close();
  //bool transf = 1;

  //cerr<<"do_assembly: done"<<endl;



}



void
MasterEquation::do_transformation(libMesh::NumericVector<Number>& u,
    libMesh::NumericVector<Number>& T,
    libMesh::NumericVector<Number>& TX, bool transf)
{


  //cerr << "do_transformation" << endl;
  //T.zero();


  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  const unsigned int Ef_e_var = system.variable_number("Ef_e");
  const unsigned int Ef_h_var = system.variable_number("Ef_h");

  libMesh::DofMap& dof_map =  system.get_dof_map();

  libMesh::NumericVector<Number>& solution = get_solution_vector();

  TiberLinearSystem* _sys_EcEv = &get_equation_system<TiberLinearSystem>(1);

  const unsigned int LUMO_var = _sys_EcEv->variable_number("LUMO");
  const unsigned int HOMO_var = _sys_EcEv->variable_number("HOMO");

  libMesh::DofMap& dof_map_E = _sys_EcEv->get_dof_map();

  libMesh::NumericVector<Number>& solution_E = _sys_EcEv->get_solution_vector();


  libMesh::FEType fe_type = system.variable_type(Ef_e_var); //

  // the volume finite element
  libMesh::UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));

  //AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));


  MeshBase::const_element_iterator el =
                                    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end =
                                    mesh.active_local_elements_end();

  const int n_el = solution_E.size();

  const int n_elem = n_el/2;

  double kbT = SimulationOptions::temperature * Constants::kb; 

  vector<unsigned int> dof_indices_tot;
  vector<unsigned int> dof_indices_tot_neigh;

  vector<unsigned int> dof_indices_n;
  vector<unsigned int> dof_indices_p;

  vector<unsigned int> dof_indices_tot_E;
  vector<unsigned int> dof_indices_Ec;
  vector<unsigned int> dof_indices_Ev;

  unsigned int ind = 0.0;
  unsigned int n = 1;



  for (; el != end; ++el)
  {
    double Ec_i = 0.0;
    double Ev_i = 0.0;
    double n_i = 0.0;
    double p_i = 0.0;

    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices_tot);
    dof_map.dof_indices(elem, dof_indices_n, Ef_e_var);
    dof_map.dof_indices(elem, dof_indices_p, Ef_h_var);

    dof_map_E.dof_indices(elem, dof_indices_Ec, LUMO_var);
    dof_map_E.dof_indices(elem, dof_indices_Ev, HOMO_var);

    Ec_i = solution_E(dof_indices_Ec[0]);
    Ev_i = solution_E(dof_indices_Ev[0]);


    n_i = u(dof_indices_n[0]);
    p_i = u(dof_indices_p[0]);


    if (transf == 1 )
    {

      TX.set(ind, Ec_i - kbT * log( ( 1 - n_i )/ n_i ));
      TX.set(ind + n_elem, Ev_i + kbT * log( ( 1 - p_i )/ p_i ));

      T.set(ind, kbT / ( n_i * ( 1 - n_i ) ));
      T.set(ind + n_elem, -kbT / ( p_i * ( 1 - p_i )));


    }
    if (transf == 0)
    {

      if (TX(ind) >= Ec_i + n*sigma)
      {
        TX.set(ind,  Ec_i + n*sigma/2);
      }
      if (TX(ind + n_elem) <= Ev_i - n*sigma)
      {
        TX.set(ind + n_elem, Ev_i - n*sigma/2 );
      }
     // if (TX(ind) <= Ev_i - n*sigma)
     // {
     //   TX.set(ind, Ev_i + kbT/2);
      //}
     // if (TX(ind + n_elem) >= Ec_i)
      //{
      //  TX.set(ind + n_elem, Ec_i - kbT/2);
      //}

      //cerr << "Ec_i ="  << Ec_i << endl;
      //cerr << "TX(0)(" << ind << ") ="  << TX(ind) << endl;
      u.set(ind, 1 / ( 1 + exp( ( Ec_i - TX(ind) ) / kbT ) ));
      u.set(ind + n_elem, 1 / ( 1 + exp( ( -Ev_i + TX(ind + n_elem ) ) / kbT ) ));

      //cerr << "u(" << ind << ") ="  << u(ind) << endl;
    }

    ind += 1;
  }



  if (transf == 1 )
  {
    //cerr << "transf = 1" << endl;
    //cerr << "u (n)" << u << endl;
    //cerr << "TX" << TX << endl;
    //cerr << "T" << T << endl;

  }

  if (transf == 0 )
  {
    //cerr << "transf = 0" << endl;
    //cerr << "u (n+1) "<< u << endl;
  }





}


void
MasterEquation::get_solution_secure(const Elem* elem, std::map<ID,
    std::vector<double> >& solutions,
    const std::vector<Point>& points)
{


  unsigned int np = points.size();

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

  const libMesh::NumericVector<Number>& solution = system->get_solution_vector();


  const unsigned int dim = get_mesh().mesh_dimension();

  const libMesh::DofMap& dof_map = system->get_dof_map();

  const unsigned int Ef_e_var = system->variable_number("Ef_e");
  const unsigned int Ef_h_var = system->variable_number("Ef_h");

  TiberLinearSystem* _sys_EcEv = &get_equation_system<TiberLinearSystem>(1);

  libMesh::NumericVector<Number>& solution_E = _sys_EcEv->get_solution_vector();

  libMesh::DofMap& dof_map_E = _sys_EcEv->get_dof_map();

  const unsigned int LUMO_var = _sys_EcEv->variable_number("LUMO");
  const unsigned int HOMO_var = _sys_EcEv->variable_number("HOMO");

  libMesh::FEType fe_type = system->variable_type(Ef_e_var);
  libMesh::UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));

  vector<unsigned int> dof_indices_n;
  vector<unsigned int> dof_indices_p;

  vector<unsigned int> dof_indices_Ec;
  vector<unsigned int> dof_indices_Ev;

  //const int n_elem = solution_E.size()/2;

  fe->reinit(elem, &points);

  dof_map.dof_indices(elem, dof_indices_n, Ef_e_var);
  dof_map.dof_indices(elem, dof_indices_p, Ef_h_var);

  dof_map_E.dof_indices(elem, dof_indices_Ec, LUMO_var);
  dof_map_E.dof_indices(elem, dof_indices_Ev, HOMO_var);

  const unsigned int n_dofs = dof_indices_n.size();

  const double beta = 1.0/(SimulationOptions::temperature * Constants::kb); 

  //const int n_elem = solution.size() /2;


  for (unsigned int e = 0; e < np; e++)
  {
    double Ef_e = 0.0;
    double Ef_h = 0.0;

    double LUMO = 0.0;
    double HOMO = 0.0;

    double Pe = 0.0;
    double Ph = 0.0;

    for (unsigned int i = 0; i < n_dofs; i++) // RIDONDANTE
    {
     Ef_e += solution(dof_indices_n[i]);
     Ef_h += solution(dof_indices_p[i]);


     LUMO += solution_E(dof_indices_Ec[i]);
     HOMO += solution_E(dof_indices_Ev[i]);

     Pe = 1/ ( 1 + exp( -( Ef_e - LUMO )*beta  ) );
     Ph = 1/ ( 1 + exp ( - ( HOMO - Ef_h )*beta ) );

    }


    if (solutions.count(eFermiLevel))
    {
      solutions[eFermiLevel][e] = Ef_e;
    }

    if (solutions.count(hFermiLevel))
    {
      solutions[hFermiLevel][e] = Ef_h;
    }

    if (solutions.count(LUMOLevel))
    {
      solutions[LUMOLevel][e] = LUMO;
    }

    if (solutions.count(HOMOLevel))
    {
      solutions[HOMOLevel][e] = HOMO;
    }

    if (solutions.count(eOccProbability))
    {
      solutions[eOccProbability][e] = Pe;
    }

    if (solutions.count(hOccProbability))
    {
      solutions[hOccProbability][e] = Ph;
    }
    
    if (solutions.count(ChargeDensity))
    {
      double u = get_mesh_units() * 1.0e2;  // units of lenth expressed in cm    
      solutions[ChargeDensity][e] = Constants::elementary_charge * 
	                            (Ph - Pe) / (elem->volume()*u*u*u);
    }
  }
}


void
MasterEquation::do_setup_solution_variables(void)
{
  // declare solutions variables
  declare_solution(eFermiLevel, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "eV");
  declare_solution(hFermiLevel, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "eV");

  declare_solution(LUMOLevel, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "eV");
  declare_solution(HOMOLevel, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "eV");

  declare_solution(eOccProbability, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "");
  declare_solution(hOccProbability, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "");

  declare_solution(ChargeDensity, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "C/cm^3");

}





