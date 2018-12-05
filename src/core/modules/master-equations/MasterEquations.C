
// modules includes
#include "MasterEquations.h"
#include "Material.h"
#include "MasterEquationsProperties.h"
#include "MEBulkModel.h"
#include "Constants.h"
#include "TiberNonlinearSystem.h"
#include "TiberLinearSystem.h"
#include "SolveFailedException.h"
#include "Variable.h"
#include "MasterEquationsModelInterface.h"



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


#include "DataOutput.h"
#include "Messages.h"



// C++ includes
#include <tr1/random>
#include <fstream>
#include <vector>

#include "TiberModule.h"

// namespaces
using namespace std;
using namespace MasterEquationsDefs;


// pointer _this initialize
MasterEquations*
MasterEquations::_this = NULL;



MasterEquations::Options::Options(void)
  : coupling(BOTH),
    solver_method(NEWTON)
{
}



MasterEquations::MasterEquations(const ModelOptions& options)
  : SimulationInterface(options),
    _rebuild_eq_system(true)
{
  // there's nothing to be done
}


MasterEquations::~MasterEquations(void)
{
  //cleanup_solver();
}


PhysicalModel*
MasterEquations::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  //string modelname;

  //modelname = options.get_option("model", "default"); // default


  //MEBulkModel* model =
    //MEBulkModel::create(modelname, mat, options);

  //if (model == NULL)
    //throw ModelErrorException(
        //"MasterEquations: No such physical model: " + modelname);

  //return model;
  return NULL;
}


//MasterEquations*
//MasterEquations::create(const std::string& name, const Material* mat,
//    const ModelOptions& options)
//{
//  return PhysicalModel:: create<MasterEquationsProperties>("me_bulk_" + name, mat, options);
//}


void
MasterEquations::do_init(void)
{
  cerr<<"do_init_begin"<<endl;

  // Create non linear equation system
  create_equation_system("nonlinear");


  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();


  // add variables "electron density" and "hole density"
  system.add_variable("n", CONSTANT, MONOMIAL);
  system.add_variable("p", CONSTANT, MONOMIAL);



  // assembly of the system
  system.attach_assembly_routine(assemble);

  // attach transformation
  system.attach_transformation_routine(transformation);


  // initialize the system
  system.init();

  cerr<<"do_init attach assembly done"<<endl;


  // Solutions vector
  NumericVector<Number>& solution = get_solution_vector();


  // define number of variable which belong to the system
  const unsigned int n_var = system.variable_number("n");
  const unsigned int p_var = system.variable_number("p");

  const double scaling = 1;

  // set initial values of densities
  double level_n = 0.0;//0.1 * scaling;//1e6;
  double level_p = 0.0;//0.1 * scaling;//2e6;

  double Ef_n = 2.8;
  double Ef_p = 0.7;
  double Ef_n0 = 2.7;

  // linear system of support to set values of energies (Ec and Ev) perturbated by gaussians
  create_equation_system("linear");
  TiberLinearSystem* _sys_EcEv = &get_equation_system<TiberLinearSystem>(1);

  // set band edges to initial values
  _sys_EcEv->add_variable("Ec",CONSTANT,MONOMIAL);
  _sys_EcEv->add_variable("Ev",CONSTANT,MONOMIAL);

  // initialize
  _sys_EcEv->init();


  // get dofs map for Ec and Ev
  const DofMap& dof_map_E = _sys_EcEv->get_dof_map();
  //const DofMap& dof_map_Ev = _sys_EcEv->get_dof_map();


  //
  NumericVector<Number>& solution_E = _sys_EcEv->get_solution_vector();

  // define number of variable which belong to the system
  const unsigned int Ec_var = _sys_EcEv->variable_number("Ec");
  const unsigned int Ev_var = _sys_EcEv->variable_number("Ev");

  const double sigma = 0.06; // sigma of Gaussian expression (standard deviation)
  const double mean = 0.00;  // mean value of gaussian

  const double sigma_c = 0.02 * scaling; // sigma of Gaussian expression (standard deviation)


  const double level_Ec = 2.5;//Ec = db->get_conduction_band_edge();

  const double level_Ev = 1;//Ev = db->get_valence_band_edge();


  // get mesh
  MeshBase& mesh = get_mesh();


  const int n_el = solution_E.size()/2;


  // define iterators which point to the active elements
  MeshBase::element_iterator el = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();


  cerr<<"do_init set level done"<<endl;

  int seed = 110101010;//= std::tr1::random_device(); //time(NULL) *
  //const double log_sigma = log(1 / sigma * sqrt(2 * M_PI));
  std::tr1::mt19937 generator(seed);

  const double KbT_elem = 0.026;

  const double dE_min = log(1e-4) * KbT_elem;

  cerr<<"do_init loop over elements"<<endl;

  for (; el != end; ++el)
  {
    const Elem* elem = *el;


    unsigned const id_n = elem->dof_number(system.number(), n_var, 0 );
    unsigned const id_p = elem->dof_number(system.number(), p_var, 0 );

    unsigned const id_Ec = elem->dof_number(_sys_EcEv->number(), Ec_var, 0 );
    unsigned const id_Ev = elem->dof_number(_sys_EcEv->number(), Ev_var, 0 );


    const double eps = std::numeric_limits<double>::min();

    double r_n;
    double rn_1;
    double rn_2;

    double r_n_c;
    double rn_1_c;
    double rn_2_c;


    for (unsigned int i = 0; i < 2; i++)
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
      else
      {
        rn_2 = r_n;
      }

    }


    for (unsigned int i = 0; i < 2; i++)
    {
      double ru1_c;
      do
      {
        ru1_c = static_cast<double>(generator()) / generator.max();
      } while (ru1_c <= eps);

      double ru2_c = static_cast<double>(generator()) / generator.max();

      r_n_c = sigma_c * sqrt(-2.0 * log(ru1_c)) * cos(2.0 * M_PI * ru2_c);

      if (i == 0)
      {
        rn_1_c = r_n_c;
      }
      else
      {
        rn_2_c = r_n_c;
      }

    }

    solution_E.set(id_Ec, level_Ec + rn_1); //
    solution_E.set(id_Ev, level_Ev + rn_2); //


    if (solution_E(id_Ec) - Ef_n < dE_min)
    {
      level_n = 1 / ( 1 + exp (( solution_E(id_Ec) - 0.5 * ( solution_E(id_Ec) + sigma ) )/KbT_elem) );
    }
    else
    {
      level_n = 1 / ( 1 + exp (( solution_E(id_Ec) - Ef_n)/KbT_elem) );
    }

    if (-solution_E(id_Ev) + Ef_p < dE_min)
    {
      level_p = 1 / ( 1 + exp ((-solution_E(id_Ev) + solution_E(id_Ep) - sigma ) )/KbT_elem ) );
    }
    else
    {
      level_p = 1 / ( 1 + exp ((-solution_E(id_Ev) + Ef_p)/KbT_elem ) );
    }

    solution.set(id_n, level_n ); //
    solution.set(id_p, level_p ); //

  }

  cerr<<"do_init loop over elements done"<<endl;

  cerr << solution << endl;
  //cerr << solution_E << endl;

  ofstream outFile;
  outFile.open("/home/drossi/TiberCAD/ME/OLED_2D/Test_1/output/Energies.txt");
  outFile << solution_E << endl;
  outFile.close();

  ofstream outFile1;
  outFile1.open("/home/drossi/TiberCAD/ME/OLED_2D/Test_1/output/Densities.txt");
  outFile1 << solution << endl;
  outFile1.close();

  // cerr<<"do_init linear"<<endl;

  //cerr<<"do_init linear init"<<endl;



  //std::random_device r;
  //cerr<<"do_init el loop"<<endl;

    //double ru1 = static_cast<double>(generator()) / generator.max();
    //double ru2 = static_cast<double>(generator()) / generator.max();

    //rn1 = -sqrt(2 * pow(sigma, 2) * (log_sigma - log(1 - ru1)));
    //rn2 = -sqrt(2 * pow(sigma, 2) * (log_sigma - log(1 - ru2)));


  cerr<<"do_init: done"<<endl;
}


void
MasterEquations::do_newton(void)
{
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>(0);

  system.set_options(get_solver_options());
  system.solve();

}




void
MasterEquations::do_solve(void)
{

  // this is dirty, but at the moment we have to provide a static method
  // for assembly. This should change in future
  _this = this;
  cerr<<"do_solve"<<endl;

  parse_options();

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
MasterEquations::rebuild_equation_system(void)
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


  system.add_variable("n", CONSTANT, MONOMIAL);
  system.add_variable("p", CONSTANT, MONOMIAL);

  //system.add_vector("me_old_sol");


  // finally initialize the newly created system
  system.init();


  _rebuild_eq_system = false;

}




void
MasterEquations::do_print_info(void)
{
  Messages::info("Module implementation example: simple Master Equations");

  parse_options();

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
MasterEquations::parse_options(void)
{

  const ModelOptions& opts = SimulationInterface::get_options();
  Options& myopts = get_options();

  string coupling = opts.get_option("coupling", "");

  if (coupling == "electrons")
    myopts.coupling = ELECTRONS;
  else if (coupling == "holes")
    myopts.coupling = HOLES;
  else
    myopts.coupling = BOTH;

}



// Coupling

void
MasterEquations::assemble(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
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
MasterEquations::transformation(NumericVector<Number>& u,
     NumericVector<Number>& T,
       NumericVector<Number>& TX, bool transf)
{
  //cerr << "transformation" << endl;

  _this->do_transformation(u, T, TX, transf);
}

template <int coupling>
void
MasterEquations::do_assembly(const NumericVector<Number>& x, NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{

  //cerr << "do_assembly: init" << endl;
  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  NumericVector<Number>& solution = get_solution_vector();


  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  unsigned int it_do_ass = 0;

  it_do_ass +=1;


  const unsigned int n_var = system.variable_number("n");
  const unsigned int p_var = system.variable_number("p");


  DofMap& dof_map =  system.get_dof_map();


  TiberLinearSystem* _sys_EcEv = &get_equation_system<TiberLinearSystem>(1);


  const unsigned int Ec_var = _sys_EcEv->variable_number("Ec");
  const unsigned int Ev_var = _sys_EcEv->variable_number("Ev");

  DofMap& dof_map_E = _sys_EcEv->get_dof_map();

  NumericVector<Number>& solution_E = _sys_EcEv->get_solution_vector();

  //NumericVector<Number>& oldx = system.get_vector("me_old_sol");



  //cerr<< solution_E << endl;

  //cerr<<"do_assembly build finite elements"<<endl;

  const double alpha = 1e-9; //1e-9;
  const double v_0 = 1; //1e12;

  FEType fe_type = system.variable_type(n_var); //

  // the volume finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));


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
  DenseMatrix<Number> Ke;
  // the system rhs (will hold also element rhs contribution)
  DenseVector<Number> Fe;
  // the local solution
  DenseVector<Number> X;


  DenseSubMatrix<Number>
    Knn(Ke), Knp(Ke),
    Kpn(Ke), Kpp(Ke);

  DenseSubVector<Number>
    Fn(Fe),
    Fp(Fe);

  DenseSubVector<Number>
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


 // cerr<<"do_assembly loop over elements begin"<<endl;

  int unsigned iter = 0;
  int unsigned iteration ;

  double n_avg = 0.0;
  double p_avg = 0.0;

  const int n_el = solution_E.size();

  const double n_elem = n_el/2;

  iteration = iteration + 1;

  n_avg = 0.1;
  p_avg = 0.1;


  double Ef_n = 2.35;
  double Ef_p = 1.15;
  double Ef_n0 = 2.6;
  double Ef_p0 = 1.15;

  const double KbT_elem = 0.026;


  //cerr<<"do_assembly: begin loop over elements"<<endl;

  // loop over all active elements
  for (; el != end; ++el)
  {

    //cerr << "element loop" << endl;

    const Elem* elem = *el;
    //cout << el << endl;
    iter += 1;

    //cerr<<"do_assembly bug1"<<endl;

    const double KbT_elem = 0.026;// [eV]  sc->get_lattice_temperature();//lattice_vt;

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices_tot);
    dof_map.dof_indices(elem, dof_indices_n, n_var);
    dof_map.dof_indices(elem, dof_indices_p, p_var);

    dof_map_E.dof_indices(elem, dof_indices_Ec, Ec_var);
    dof_map_E.dof_indices(elem, dof_indices_Ev, Ev_var);


    unsigned int n_dofs     = dof_indices_n.size(); // number of dofs of element for single variable
    unsigned int n_dofs_tot = dof_indices_tot.size();   // number of dofs of element for all variables


    //cerr << n_dofs << endl;
    //cerr << n_dofs_tot<< endl;


    dof_indices_tot.reserve(30);
    dof_indices_tot_neigh.reserve(30);
    dof_indices_n.reserve(30);
    dof_indices_p.reserve(30);
    dof_indices_Ec.reserve(30);
    dof_indices_Ev.reserve(30);



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
        dof_map.dof_indices(elem_neig, neig_dof_indices_n, n_var);
        dof_map.dof_indices(elem_neig, neig_dof_indices_p, p_var);

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


    double Ec_i = solution_E(dof_indices_Ec[0]); // Ec of element
    double Ev_i = solution_E(dof_indices_Ev[0]); // Ev of element


    // First we will build the system matrix Ke_ij
    if (jacobian != NULL)
    {


      double Ec_i = solution_E(dof_indices_Ec[0]); // Ec of element
      double Ev_i = solution_E(dof_indices_Ev[0]); // Ev of element

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
            double Ec_j = solution_E(dof_indices_Ec[ind]);

            // get valence band edge for the calculation of Miller-Abrhams energetic term
            double Ev_j = solution_E(dof_indices_Ev[ind]);

            // from site j to site i
            double delta_Ec_ji = Ec_j - Ec_i;
            double delta_Ev_ji = Ev_j - Ev_i;

            // from site i to site j
            double delta_Ec_ij = Ec_i - Ec_j;
            double delta_Ev_ij = Ev_i - Ev_j;

            // from site j to site i
            double exp_Ec_ji = 0.0;
            double exp_Ev_ji = 0.0;

            // from site i to site j
            double exp_Ec_ij = 0.0;
            double exp_Ev_ij = 0.0;

            if (delta_Ec_ji > 0)
            {
              exp_Ec_ji = 1;

              exp_Ec_ij = exp( delta_Ec_ij / KbT_elem);
            }
            if (delta_Ec_ji < 0)
            {
              exp_Ec_ji = exp( -delta_Ec_ij / KbT_elem);

              exp_Ec_ij = 1;
            }
            if (delta_Ec_ji == 0)
            {
              exp_Ec_ji = 1;

              exp_Ec_ij = 1;
            }
            if (delta_Ev_ji < 0)
            {
              exp_Ev_ji = 1;

              exp_Ev_ij = exp( -delta_Ev_ij / KbT_elem);
            }
            if (delta_Ev_ji > 0)
            {
              exp_Ev_ji = exp( delta_Ev_ij / KbT_elem);

              exp_Ev_ij = 1;
            }
            if (delta_Ev_ji == 0)
            {
              exp_Ev_ji = 1;

              exp_Ev_ij = 1;
            }


            if (coupling &BOTH)
            {
              Knn(i,i) +=  v_0 * exp_d * ( x(dof_indices_n[ind]) * (exp_Ec_ij - exp_Ec_ji) - exp_Ec_ij );  // diagonal  Jp - Jw

              Knn(i,ind) += v_0 * exp_d * ( exp_Ec_ji - x(dof_indices_n[i]) * ( exp_Ec_ij - exp_Ec_ji ) ); //off diagonal -> Wji - Pi * (  Wij - Wji )

              Kpp(i,i) +=  v_0 * exp_d * ( x(dof_indices_p[ind]) * (exp_Ev_ij - exp_Ev_ji) - exp_Ev_ij );  // diagonal Jp - Jw

              Kpp(i,ind) += v_0 * exp_d * ( exp_Ev_ji - x(dof_indices_p[i]) * ( exp_Ev_ij - exp_Ev_ji ) ); //off diagonal -> Wji - Pi * (  Wij - Wji )
            }

          }
        }
      }
    }


    if (jacobian != NULL)
    {
      //cerr << "Knn" << endl << Knn << endl;
      //cerr << "iter =" << iter << endl;
      //cerr << "Kpp" << endl << Kpp << endl;
    }

    if (residual != NULL)
    {
      for (unsigned int i = 0; i < n_dofs ; i++) // trovare funzione per calc num di dofs di un elem
      {

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
            double Ec_j = solution_E(dof_indices_Ec[ind]);

            // get valence band edge for the calculation of Miller-Abrhams energetic term
            double Ev_j = solution_E(dof_indices_Ev[ind]);

            // from site j to site i
            double delta_Ec_ji = Ec_j - Ec_i;
            double delta_Ev_ji = Ev_j - Ev_i;

            // from site i to site j
            double delta_Ec_ij = Ec_i - Ec_j;
            double delta_Ev_ij = Ev_i - Ev_j;

            // from site j to site i
            double exp_Ec_ji;
            double exp_Ev_ji;

            // from site i to site j
            double exp_Ec_ij;
            double exp_Ev_ij;

            if (delta_Ec_ji > 0)
            {
              exp_Ec_ji = 1;

              exp_Ec_ij = exp( delta_Ec_ij / KbT_elem);
            }
            if (delta_Ec_ji < 0)
            {
              exp_Ec_ji = exp( -delta_Ec_ij / KbT_elem);

              exp_Ec_ij = 1;
            }
            if (delta_Ec_ji == 0)
            {
              exp_Ec_ji = 1;

              exp_Ec_ij = 1;
            }
            if (delta_Ev_ji < 0)
            {
              exp_Ev_ji = 1;

              exp_Ev_ij = exp( -delta_Ev_ij / KbT_elem);
            }
            if (delta_Ev_ji > 0)
            {
              exp_Ev_ji = exp( delta_Ev_ij / KbT_elem);

              exp_Ev_ij = 1;
            }
            if (delta_Ev_ji == 0)
            {
              exp_Ev_ji = 1;

              exp_Ev_ij = 1;
            }


            if (coupling & BOTH)
            {
              Fn(i) += v_0 * exp_d * x(dof_indices_n[ind]) * ( exp_Ec_ji + x(dof_indices_n[i]) * ( exp_Ec_ij - exp_Ec_ji ) ); // + Kij * Pj = Pj * ( Wji - Pi * ( Wij - Wji ) )

              Fn(i) -= v_0 * exp_d * exp_Ec_ij * x(dof_indices_n[i]); // - Fi = - Pi * sum (Wij)

              Fp(i) += v_0 * exp_d * x(dof_indices_p[ind]) * ( exp_Ev_ji + x(dof_indices_p[i]) * ( exp_Ev_ij - exp_Ev_ji ) ); // Kij * Pj = Wji - Pi * ( Wij - Wji )

              Fp(i) -= v_0 * exp_d * exp_Ev_ij * x(dof_indices_p[i]); // - Fi = - Pi * sum (Wij)
            }

          }

        }
      }

      residual->add_vector(Fe, dof_indices_tot);

    }
    else
    {
      jacobian->add_matrix(Ke, dof_indices_tot, dof_indices_tot_neigh);
    }


  }// end loop over all active elements


// boundary conditions jacobian
if (jacobian != NULL)
{
  const int n_el = solution_E.size();

  const double n_elem = n_el/2;
      for (unsigned int l = 0; l < n_elem  ; l++)
      {
        //cerr << n_elem << endl;
        //jacobian->set(0, l , 1 / n_elem);
        jacobian->set(n_elem - 1, l , 1 / n_elem);
        jacobian->set(n_el - 1, l + n_elem, 1 / n_elem);
        //jacobian->set(n_elem - 1, l , 0);
        //jacobian->set(n_el - 1, l + n_elem, 0);
      }
      //jacobian->set(n_elem - 1, n_elem-1,  1);
      //jacobian->set(n_el - 1, n_el-1, 1);
  }


// boundary conditions residual
if (residual != NULL)
{
  const int n_el = solution_E.size();

  const double n_elem = n_el/2;

  double BCn = 0.0;
  double BCp = 0.0;
  double BCn_0 = 0.0;
  double BCp_0 = 0.0;

  for (unsigned int j = 0; j < n_elem; j++)
  {
    //BCn_0 += ( 1 / n_elem )* x(j);
    BCn += ( 1 / n_elem )* x(j);
    BCp += ( 1 / n_elem )* x(j + n_elem);
  }

  for (unsigned int j = 0; j < n_elem; j++)
  {
    if (j == 0)
    {
      //BCn_0 -= 1 / (( 1 + exp (( solution_E(j) - Ef_n0)/KbT_elem) ));
    }


    ///BCn_0 -= 1 / (( 1 + exp (( solution_E(j) - Ef_n)/KbT_elem) ));
    BCn -= 1 / (n_elem * ( 1 + exp (( solution_E(j) - Ef_n)/KbT_elem) )); //n_avg;
    BCp -= 1 / (n_elem * ( 1 + exp (( -solution_E(n_elem + j ) + Ef_p)/KbT_elem) ));
  }
  //BCn -= 1 / ( 1 + exp (( solution_E(n_elem - 1) - Ef_n)/KbT_elem) ); //n_avg;
  //BCp -= 1 / ( 1 + exp (( -solution_E(n_el - 1) + Ef_p)/KbT_elem) );


  //residual->set(0, BCn_0);
  residual->set(n_elem - 1, BCn);
  residual->set(n_el - 1, BCn);


}

if (residual != NULL)
{
  //cerr << "Fe = " << Fe << endl;
  //cerr << "Fp" << Fn << endl;
}

    //cerr << "VECTOR X" << endl << X << endl;
    //cerr << "MATRIX Ke" << endl << Ke << endl;
    //cerr << "VECTOR Fe" << endl << Fe << endl;

  //cerr <<  "n_tot = " << n_tot << endl;
  //cerr <<  "p_tot = " << p_tot << endl;

   //cerr << "VECTOR x" << x << endl;

   if (jacobian != NULL)
   {
     jacobian->print_matlab("/home/drossi/TiberCAD/ME/OLED_2D/Test_1/output/jacobian.m");
   }
   if (residual != NULL)
   {
     residual->print_matlab("/home/drossi/TiberCAD/ME/OLED_2D/Test_1/output/residual.m");
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




  ofstream outFile_x;
  outFile_x.open("/home/drossi/TiberCAD/ME/OLED_2D/Test_1/output/solution.txt");
  outFile_x << x << endl;
  outFile_x.close();


  ofstream outFile_res;
  outFile_res.open("/home/drossi/TiberCAD/ME/OLED_2D/Test_1/output/Jac_"+ num + ".txt");
  outFile_res << Ke << endl;
  outFile_res.close();

  //bool transf = 1;

  //cerr<<"do_assembly: done"<<endl;



}



void
MasterEquations::do_transformation(NumericVector<Number>& u,
    NumericVector<Number>& T,
       NumericVector<Number>& TX, bool transf)
{


  //cerr << "do_transformation" << endl;
  //T.zero();


  TiberNonlinearSystem& system = get_equation_system<TiberNonlinearSystem>();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  const unsigned int n_var = system.variable_number("n");
  const unsigned int p_var = system.variable_number("p");

  DofMap& dof_map =  system.get_dof_map();

  NumericVector<Number>& solution = get_solution_vector();

  TiberLinearSystem* _sys_EcEv = &get_equation_system<TiberLinearSystem>(1);

  const unsigned int Ec_var = _sys_EcEv->variable_number("Ec");
  const unsigned int Ev_var = _sys_EcEv->variable_number("Ev");

  DofMap& dof_map_E = _sys_EcEv->get_dof_map();

  NumericVector<Number>& solution_E = _sys_EcEv->get_solution_vector();


  FEType fe_type = system.variable_type(n_var); //

  // the volume finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));


  MeshBase::const_element_iterator el =
                                    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end =
                                    mesh.active_local_elements_end();

  const int n_el = solution_E.size();

  const int n_elem = n_el/2;

  const double KbT_elem = 0.026;

  vector<unsigned int> dof_indices_tot;
  vector<unsigned int> dof_indices_tot_neigh;

  vector<unsigned int> dof_indices_n;
  vector<unsigned int> dof_indices_p;

  vector<unsigned int> dof_indices_tot_E;
  vector<unsigned int> dof_indices_Ec;
  vector<unsigned int> dof_indices_Ev;

  unsigned int ind = 0.0;
  unsigned int n = 1;
  double sigma = 0.06;



  for (; el != end; ++el)
  {
    double Ec_i = 0.0;
    double Ev_i = 0.0;
    double n_i = 0.0;
    double p_i = 0.0;

    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices_tot);
    dof_map.dof_indices(elem, dof_indices_n, n_var);
    dof_map.dof_indices(elem, dof_indices_p, p_var);

    dof_map_E.dof_indices(elem, dof_indices_Ec, Ec_var);
    dof_map_E.dof_indices(elem, dof_indices_Ev, Ev_var);

    Ec_i = solution_E(dof_indices_Ec[0]);
    Ev_i = solution_E(dof_indices_Ev[0]);


    n_i = u(dof_indices_n[0]);
    p_i = u(dof_indices_p[0]);


    if (transf == 1 )
    {

      TX.set(ind, Ec_i - KbT_elem * log( ( 1 - n_i )/ n_i ));
      TX.set(ind + n_elem, Ev_i + KbT_elem * log( ( 1 - p_i )/ p_i ));

      T.set(ind, KbT_elem / ( n_i * ( 1 - n_i ) ));
      T.set(ind + n_elem, -KbT_elem / ( p_i * ( 1 - p_i )));


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
     //   TX.set(ind, Ev_i + KbT_elem/2);
      //}
     // if (TX(ind + n_elem) >= Ec_i)
      //{
      //  TX.set(ind + n_elem, Ec_i - KbT_elem/2);
      //}

      //cerr << "Ec_i ="  << Ec_i << endl;
      //cerr << "TX(0)(" << ind << ") ="  << TX(ind) << endl;
      u.set(ind, 1 / ( 1 + exp( ( Ec_i - TX(ind) ) / KbT_elem ) ));
      u.set(ind + n_elem, 1 / ( 1 + exp( ( -Ev_i + TX(ind + n_elem ) ) / KbT_elem ) ));

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
MasterEquations::get_solution_secure(const Elem* elem, std::map<ID,
    std::vector<double> >& solutions,
    const std::vector<Point>& points)
{


  unsigned int np = points.size();

  TiberNonlinearSystem* system = &get_equation_system<TiberNonlinearSystem>();

  const NumericVector<Number>& solution = system->get_solution_vector();


  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system->get_dof_map();

  const unsigned int n_var = system->variable_number("n");
  const unsigned int p_var = system->variable_number("p");

  TiberLinearSystem* _sys_EcEv = &get_equation_system<TiberLinearSystem>(1);

  NumericVector<Number>& solution_E = _sys_EcEv->get_solution_vector();

  DofMap& dof_map_E = _sys_EcEv->get_dof_map();

  const unsigned int Ec_var = _sys_EcEv->variable_number("Ec");
  const unsigned int Ev_var = _sys_EcEv->variable_number("Ev");

  FEType fe_type = system->variable_type(n_var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices_n;
  vector<unsigned int> dof_indices_p;

  vector<unsigned int> dof_indices_Ec;
  vector<unsigned int> dof_indices_Ev;

  const int n_elem = solution_E.size()/2;

  fe->reinit(elem, &points);

  dof_map.dof_indices(elem, dof_indices_n, n_var);
  dof_map.dof_indices(elem, dof_indices_p, p_var);

  dof_map_E.dof_indices(elem, dof_indices_Ec, Ec_var);
  dof_map_E.dof_indices(elem, dof_indices_Ev, Ev_var);

  const unsigned int n_dofs = dof_indices_n.size();


  for (unsigned int e = 0; e < np; e++)
  {
    double n = 0.0;
    double p = 0.0;

    double Ec = 0.0;
    double Ev = 0.0;

    for (unsigned int i = 0; i < n_dofs; i++) // RIDONDANTE
    {
      n += solution(dof_indices_n[i]);
      p += solution(dof_indices_p[i]);

      Ec += solution_E(dof_indices_Ec[i]);
      Ev += solution_E(dof_indices_Ev[i]);
    }


    if (solutions.count(eDensity))
    {
      solutions[eDensity][e] = n;
    }

    if (solutions.count(hDensity))
    {
      solutions[hDensity][e] = p;
    }

    if (solutions.count(Ec_edge))
    {
      solutions[Ec_edge][e] = Ec;
    }

    if (solutions.count(Ev_edge))
    {
      solutions[Ev_edge][e] = Ev;
    }
  }
}


void
MasterEquations::do_setup_solution_variables(void)
{
  // declare solutions variables
  declare_solution(eDensity, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "cm^-3");
  declare_solution(hDensity, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "cm^-3");

  declare_solution(Ec_edge, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "eV");
  declare_solution(Ev_edge, SolutionDescriptor::REAL, SolutionDescriptor::CELL, "eV");

}


void MasterEquations::get_solution_secure(std::map<ID, std::vector<double> >& solutions)
{

}







