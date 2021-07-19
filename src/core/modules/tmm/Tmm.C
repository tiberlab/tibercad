// $Id$

#include "Tmm.h"
#include "TiberLinearSystem.h"
#include "TmmBulkModel.h"
#include "Messages.h"


#include "libmesh/dof_map.h"


// This is needed in order to create the shared module library
#include "TiberModule.h"

#include "function.h"

using namespace libMesh;
using std::vector;


Tmm::Tmm(const ModelOptions& options) :
  SimulationInterface(options)
{
}


Tmm::~Tmm(void)
{
}


Tmm*
Tmm::create(const ModelOptions& options)
{
  return new Tmm(options);
}



void
Tmm::do_init(void)
{
  parse_options();

  // create a linear equation system 
  create_equation_system("linear");

  // get the reference to it
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // add variables and attach the assemble function
  system.add_variable("E", libMeshEnums::CONSTANT, MONOMIAL, &get_region_ids());
 
  system.init();
}


void
Tmm::parse_options(void)
{
  // read wavelengths from input
  get_option("wavelengths", _wavelengths);
  if (_wavelengths.empty())
  {
    Messages::warning("You did not provide any wavelengths for TMM.");
  }
  get_option("incident_angle", _incident_angle);
  if (_incident_angle.empty())
  {
    Messages::warning("You did not provide any incident_angle for TMM.");
  }


}


void
Tmm::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(EField, REAL, CELL, "V/cm");
  declare_solution(HField, VECTOR, CELL, "A/cm");
  //declare_solution(Displacement, VECTOR, CELL, "C/cm^2");

  // we can define aliases but the same name cannot refer to
  // several IDs
  //add_alias("ElectricField", EField);
}


void
Tmm::do_solve(void)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  NumericVector<libMesh::Number>& solution = system.get_local_solution_vector();
  solution.close();
  solution.zero();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  // NB: Tibercad by default uses length-scale in meters 
  //     This means that FEM derivatives d/dx are in 1/m
  //     To change this behavior it is necessary to define a different 'scaling'
  //     For instance if we want to use mesh_units in the assembly we need to:
  //     1. set the scaling to mesh units:
  //        get_scaling().set_length_scaling(get_mesh_units());
  //     2. use  build_finite_element(dim, fe_type, true)  
  //                                                ^ false is the default  
  //                                                
  // Now 2nd derivatives will be 1/mesh_units^2
  // We need a factor to transform rho/eps0 into V/mesh_units^2
  // Charge density is cm^-3, and Constants::e is in Coulomb, 
  // Constant::e0 is in C/Vm
  // The factor Lambda is such that rho*Lambda is in V/mesh_units^2
  // BUT (BUT) 
  // This is not that clever! Since Displacement and Polarization are already in C/m^2
  // it is easier to work with the derivatives in 1/m and rho/eps0 in V/m^2 
  // The factor 1e6 is for cm^3 -> m^3 in rho 
  get_scaling().set_length_scaling(1.0);

  DofMap& dof_map =  system.get_dof_map();
  vector<unsigned int> dof_indices;

  for (unsigned int i = 0; i < _wavelengths.size(); ++i)   //loop over wavelength
  {

  for (unsigned int j = 0; j < _incident_angle.size(); ++j)  //loop over incident angle
  {

    double lambda = _wavelengths[i];


    const unsigned int uvar = system.variable_number("E");

    MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

    double E = 1.0;
    double incoming_angle=_incident_angle[j];

    vector<double> n_real;
    vector<double> n_imag;
    vector<double> l_length;
    vector<double> l;


    for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;

      dof_map.dof_indices(elem, dof_indices, uvar);
      const unsigned int n_dofs = dof_indices.size();

      TmmBulkModel& mod = *get_bulk_model<TmmBulkModel>(elem);

      mod.reinit(elem);


      
      /*      getting refractive index and length of layers from model
      n_imag.push_back(sqrt((abs(mod.get_permittivity(lambda))-real(mod.get_permittivity(lambda)))/2));
      n_real.push_back(sqrt((abs(mod.get_permittivity(lambda))+real(mod.get_permittivity(lambda)))/2));
      l.push_back(elem->volume()); 
      */

    }
    //******************************************************
    // manully defining and printing refractive index and length of layers 
    for(int nm=0;nm<1;nm++)
       n_real.push_back(1);
    for(int nm=1;nm<2;nm++)
       n_real.push_back(1.3);
    for(int nm=2;nm<3;nm++)
       n_real.push_back(1.8);
    for(int nm=3;nm<4;nm++)
       n_real.push_back(1.9);
    
    for(int nm=0;nm<n_real.size();nm++)
       std::cout<<n_real[nm]<<" n "<<endl;
    //****
    for(int nm=0;nm<1;nm++)
       l.push_back(1e-6);
    for(int nm=1;nm<2;nm++)
       l.push_back(1e-6);
    for(int nm=2;nm<3;nm++)
       l.push_back(1e-6);
    for(int nm=3;nm<4;nm++)
       l.push_back(1e-6);

    for(int nm=0;nm<l.size();nm++)
       std::cout<<l[nm]<<" l "<<endl;
    //****
    for(int nm=0;nm<1;nm++)
       n_imag.push_back(0.0);
    for(int nm=1;nm<2;nm++)
       n_imag.push_back(0.0);
    for(int nm=2;nm<3;nm++)
       n_imag.push_back(0.0);
    for(int nm=3;nm<4;nm++)
       n_imag.push_back(0.0);

    for(int nm=0;nm<n_imag.size();nm++)
       std::cout<<n_imag[nm]<<" k "<<endl;
    //***************************************************
    // snell's law
    vector<double> theta(n_real.size());
    theta=theta_cal(n_real,incoming_angle);

    // defining Vectors 
    vector<vector<complex<double>>> I {{1,0},{0,1}};
    vector<vector<complex<double>>> M {{0,0},{0,0}};
    vector<vector<complex<double>>> T_load {{0,0},{0,0}};
    vector<vector<complex<double>>> T {{1,0},{0,1}};

    vector<vector<complex<double>>> E_N {{1,0},{0,0}};
    vector<vector<complex<double>>> E_I {{0,0},{0,0}};

    vector<complex<double>> E_F(n_real.size());
    vector<complex<double>> E_B(n_real.size());
    E_F[n_real.size()-1]=E_N[0][0];
    E_B[n_real.size()-1]=E_N[1][0];

    vector<complex<double>> E_F_NORM(n_real.size());
    vector<complex<double>> E_B_NORM(n_real.size());

    // main loop over layer, calculating matrixs
    for (int k=n_real.size()-1 ; k>=0 ; k--){

      if(k<(n_real.size()-1)){      
         I=get_D(n_real[k],n_imag[k],n_real[k+1],n_imag[k+1],theta[k],theta[k+1]);
         cout<<"D matrix "<<k+1<<"&"<<k+2<<endl;
         show_matrix(I);
      }
      M=get_M(n_real[k],n_imag[k],l[k],lambda,theta[k]);  
      cout<<"M matrix "<<k+1<<endl;
      show_matrix(M);  
      T_load=matrix_product(I,M);
      T=matrix_product(T,T_load);

      E_I=matrix_product(T,E_N);
      E_F[k]+=E_I[0][0];
      E_B[k]+=E_I[1][0];
    }
    cout<<"T matrix "<<endl;
    show_matrix(T);


    // reflection and transmission calculation     
    complex<double> Reflection,Transmission;

    Reflection=pow(abs(T[1][0]/T[0][0]),2);
    complex<double> nc_first (n_real[0],n_imag[0]);
    complex<double> nc_last (n_real[n_real.size()-1],n_imag[n_imag.size()-1]);
    complex<double> ratio_complex;
    ratio_complex=((nc_last)*cosd(theta[theta.size()-1]))/(nc_first*cosd(theta[0]));
    Transmission=ratio_complex*pow(abs(1.0/T[0][0]),2);
    
    // printing tansmission adnd reflection
    cout<<"trasmision is :"<< Transmission << "reflection is :"<<Reflection<<endl;  

    
    // normalizing electric field matrix
    for(double nm=0; nm < E_F.size() ; nm++){	
	E_F_NORM[nm]=E_F[nm]/E_F[0];			
	E_B_NORM[nm]=E_B[nm]/E_F[0];	
    } 

    for(int mm=0;mm<E_F_NORM.size();mm++)
        cout<<"E is "<< E_F_NORM[mm]<<endl;

   // printing electric field
    for(double nm=0; nm <= l_length.size() ; nm++){
	//E=E_F_NORM[nm]+E_B_NORM[nm];	
	E=l_length[nm];                   //just for test
	solution.add(l_length[nm], E);
    } 

  }

  }

  solution.close();
  system.update();
  //solution.print_matlab("sol.m");

}



void
Tmm::do_print_info(void)
{
  Messages::info("1D TMM");
}


PhysicalModel*
Tmm::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  return TmmBulkModel::create(mat, options);
}



PhysicalModel*
Tmm::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{
  // for now there is no boundary model 
  return(nullptr);
  //return(TmmBoundaryModel::create(boundary, options));
}




void
Tmm::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int u_var = system.variable_number("E");

  FEType fe_type = system.variable_type(u_var);
  UniquePtr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();
  const vector<Point>& real_pts = fe->get_xyz();

  ID subdomain = elem->subdomain_id();

  fe->reinit(elem, &p);

  dof_map.dof_indices(elem, dof_indices, u_var);
  const unsigned int n_dofs = dof_indices.size();

  // cell data variable
  RealGradient field(0);
  TmmBulkModel& mod = *get_bulk_model<TmmBulkModel>(elem);


  for (unsigned int n = 0; n < np; n++)
  {
    //double efield  = 0.0;

    //if (values.count(EField))
    //  values[EField][n] = efield;

  }


  if (values.count(EField))
  {
    values[EField][0] = solution(dof_indices[0]);
  //  values[EField][1] = field(1) / np;
  //  values[EField][2] = field(2) / np;
  }
}

