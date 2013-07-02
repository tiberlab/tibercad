// $Id$

#include "Elasticity.h"
#include "ElasticityModel.h"
#include "ElasticityBoundaryModel.h"
#include "TiberLinearSystem.h"
#include "Messages.h"
#include "equation_systems.h"
#include "dof_map.h"
#include "quadrature_gauss.h"
#include "quadrature_trap.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "fe_interface.h"
#include "SimulationOptions.h"
#include "SimulationEnvironment.h"
#include "TensorOperators.h"
#include "RotatedCrystal.h"
#include "AtomisticStructure.h"
#include "QuantumContact.h"
#include "DataOutput.h"

#include "TiberModule.h"


using namespace std;

Elasticity*
Elasticity::_this = NULL;



Elasticity::Elasticity(const ModelOptions& options) :
  SimulationInterface(options)
{
  // there's nothing to be done
}


Elasticity::~Elasticity(void)
{
  // there's nothing to be done
}


Elasticity*
Elasticity::create(const ModelOptions& options)
{
  // we could use the options to create different implementations
  // or something like that.
  return new Elasticity(options);
}




void
Elasticity::do_init(void)
{

  parse_options();
 
  create_equation_system("linear");
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();

  Device*  _device = &get_environment().get_device();
  
  system.add_variable("ux", FIRST);
  system.add_variable("uy", FIRST);
  system.add_variable("uz", FIRST);

  uvar.resize(3);
  uvar[0] =  system.variable_number("ux");
  uvar[1] =  system.variable_number("uy");
  uvar[2] =  system.variable_number("uz");


  system.attach_assemble_function(assemble);
  system.init();

  int n_elem = 0;
  
  //Create node connection
  const unsigned int nn  = get_mesh().n_nodes();
  node_conn.resize(nn);
  {
    vector<unsigned short int> node_conn_local(node_conn.size());
    
    
    MeshBase::const_element_iterator       el     = get_mesh().active_elements_begin();
    const MeshBase::const_element_iterator end_el = get_mesh().active_elements_end();
    
    for ( ; el != end_el; ++el, ++n_elem)
      for (unsigned int n = 0; n < (*el)->n_nodes(); n++)
	node_conn_local[(*el)->node(n)]++;

    node_conn = node_conn_local;
  }
  
  // prepare the vector for the accumulated strain
  // TODO this would not work very well with mesh refinement and needs to be improved,
  // e.g. by using a helper System
  _accumulated_strain.clear();


}


void
Elasticity::parse_options(void)
{

  const ModelOptions& opt = get_options();

  myopt.shape_error = opt.get_option("shape_error",1e-2);
  myopt.non_linear_strain = opt.get_option("non_linear_strain",false);
  myopt.non_linear_strain = opt.get_option("nonlinear_strain", myopt.non_linear_strain);
  myopt.shape_iterations = opt.get_option("shape_iterations",20);
  //myopt.deformation = opt.get_option("do_deformation",false);
  myopt.magnification = opt.get_option("magnification",1);
  myopt.structure_to_be_strained = opt.get_option("strain_atomistic_structure", "all");

}


void
Elasticity::do_setup_solution_variables(void)
{
  // we declare our solution variables

  declare_solution(Strain, TENSOR, NODES, "");
  declare_solution(StrainCell, TENSOR, CELL, "");
  declare_solution(StrainCrystal, TENSOR, NODES, "");
  //declare_solution(Energy, REAL, NODES, "Joule");
  declare_solution(Stress, TENSOR, NODES, "GPa");
  declare_solution(StressCrystal, TENSOR, NODES, "");
  declare_solution(Displacement, VECTOR, NODES, "m");
  declare_solution(StrainSource, TENSOR, NODES, "");
  declare_solution(StressSource, TENSOR, NODES, "GPa");
  declare_solution(ForceSource, TENSOR, NODES, "N/m3 ");
  declare_solution(EnergyDensity, REAL, NODES, "J/m^3");

}


void
Elasticity::do_solve(void)
{
  _this = this;

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // NOTE: this probably breaks shape deformation with lattice mismatch
  //       across different simulations. Accumulated strain should persist,
  //       and also sol should accumulate the deformation across several
  //       solves.
  sol =  (system.solution)->clone();
  sol->zero();
  _accumulated_strain.clear();
  (system.solution)->zero();

  shape_iteration = 0;
  double error_energy = 0.0;
  double error_u = 0.0;
  double tot_norm = 0.0;
  double energy = 0.0;

  //if ((SimulationOptions::verbose() > 2) &&  myopt.shape_iterations>1)
  // {
  //  cout<<"| Iteration | Elastic energy [J] | norm(u) [m] | Energy error [%] | Displ error [%] |"<<endl;
  //  cout<<"-----------------------------------------------------------------------------------------------"<<endl;
  // }


  int max_iterations = 0;
  if (myopt.non_linear_strain == true)
     max_iterations = myopt.shape_iterations;


  do {

    accumulate_strain();
    apply_shape_deformation();
    system.solution->zero();

    system.solve();
    sol->add(1.0,*(system.solution));
    
    double tot_norm = sol->l2_norm();
    double norm = (system.solution)->l2_norm();
    error_u = norm/tot_norm;

    //The error is based on the elastic energy
    //double elastic_energy = abs(compute_elastic_energy());
    //error_energy = abs((new_energy - energy)/energy) * 100.0;
    //energy = new_energ�y;

    if ((verbose() > 1) && shape_iteration > 0)
    {
      ostringstream os;
      os << "iteration " << shape_iteration << ":  Error =  " << error_u;
      Messages::info(os.str());
    }

    // apply_shape_deformation() leads to negative jacobian for certain structures??
    shape_iteration += 1;

 // } while (error_u > max_error);

  } while ((error_u > myopt.shape_error) && (shape_iteration < max_iterations));

  //double elastic_energy = abs(compute_elastic_energy());
  //if ((SimulationOptions::verbose() > 2) && myopt.shape_iterations > 1)
  //  cout<<"Elastic Energy: "<<elastic_energy<<" J"<<endl;

}


void
Elasticity::do_print_info(void)
{
  //Messages::info("ELASTICITY MODULE");
}


PhysicalModel*
Elasticity::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
 
  return ElasticityModel::create(mat, options);
}



PhysicalModel*
Elasticity::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{
  return ElasticityBoundaryModel::create(boundary, options);
}



void
Elasticity::accumulate_strain(void)
{

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();
  const NumericVector<Number>& solution = system.get_solution_vector();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  DofMap& dof_map =  system.get_dof_map();
  vector<vector<unsigned int> > dof_indices_vec(3);
  vector<unsigned int> dof_indices;

  FEType fe_type = dof_map.variable_type(uvar[0]);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  const vector<Point>& ref_points = qrule.get_points();
  const vector<Real>& JxW = fe->get_JxW();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    fe->reinit(elem);

    dof_map.dof_indices(elem, dof_indices_vec[0], uvar[0]);
    dof_map.dof_indices(elem, dof_indices_vec[1], uvar[1]);
    dof_map.dof_indices(elem, dof_indices_vec[2], uvar[2]);
    const unsigned int n_dofs = dof_indices_vec[0].size();

    RealTensor strain(0);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        for (int l = 0; l < n_dofs; l++)
          strain(i,j) += 0.5 * ((solution)(dof_indices_vec[i][l]) * dphi[l][0](j) +
              (solution)(dof_indices_vec[j][l]) * dphi[l][0](i));

    _accumulated_strain[elem] += strain;
  }
}


void
Elasticity::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
 
   unsigned int np = p.size();
   
 
   TiberLinearSystem* system = &get_equation_system<TiberLinearSystem>();

   const NumericVector<Number>& solution = system->get_solution_vector();
   const NumericVector<Number>& accum_sol = *sol;
   const unsigned int dim = get_mesh().mesh_dimension();

   const DofMap& dof_map = system->get_dof_map();

   std::vector<std::vector<unsigned int> > dof_indices(3);
   for (unsigned int i = 0; i< 3 ; i++)
     dof_map.dof_indices(elem, dof_indices[i],uvar[i]);

   FEType fe_type = system->variable_type(uvar[0]);
   AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
   const std::vector<std::vector<Real> >& phi = fe->get_phi();
   const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
   const std::vector<Point>& real_pts = fe->get_xyz();
  
   fe->reinit(elem, &p);
  
   ElasticityModel& mod = *get_bulk_model<ElasticityModel>(elem);
   
   RealTensor total_strain;
   RealTensor total_stress;
   RealTensor stress;
   RealTensor strain;
   RealTensor summed_strain;
   _accumulated_strain[elem].get_tensor(summed_strain);

  
   for (int n = 0; n < np; n++)
   {
     mod.calculate(elem, p[n]);
     
     const Tensor4DSym& C = mod.get_stiffness();
     const RealGradient& force_source =  mod.get_force_source();
     const RealTensor& strain_source  =  mod.get_strain_source();
     const RealTensor& stress_source  =  mod.get_stress_source();
    
     RealGradient u(0);
     for (unsigned int i = 0;i<3; i ++)
       for (unsigned int alpha = 0; alpha<dof_indices[i].size() ;alpha ++)
         u(i) += (accum_sol)(dof_indices[i][alpha]) * phi[alpha][n];

     //------Strain-------------------------

     for (unsigned int i = 0; i < 3; i++)
     {
       double der1 = 0;
       double der2 = 0;
       double der3 = 0;
       for (unsigned int alpha = 0; alpha < dof_indices[i].size(); alpha++)
       {
         double ui = (solution)(dof_indices[i][alpha]);
         der1 += ui * dphi[alpha][n](0);
         der2 += ui * dphi[alpha][n](1);
         der3 += ui * dphi[alpha][n](2);
       }

       strain(i,0) = der1;
       strain(i,1) = der2;
       strain(i,2) = der3;
     }
     strain(0,1) = strain(1,0) = 0.5 * (strain(0,1) + strain(1,0));
     strain(0,2) = strain(2,0) = 0.5 * (strain(0,2) + strain(2,0));
     strain(1,2) = strain(2,1) = 0.5 * (strain(1,2) + strain(2,1));


     strain += summed_strain;
     stress = C*strain;
     total_strain = strain + strain_source;
     total_stress = C * total_strain + stress_source;

     //----Displacemet--------
     if (values.count(Displacement))
     {
       values[Displacement][3*n]   = u(0);
       values[Displacement][3*n+1] = u(1);
       values[Displacement][3*n+2] = u(2);
     }

     //Total Strain
     if (values.count(Strain))
     {
       values[Strain][6*n] =   total_strain(0,0);
       values[Strain][6*n+1] = total_strain(1,1);
       values[Strain][6*n+2] = total_strain(2,2);
       values[Strain][6*n+3] = total_strain(1,0);
       values[Strain][6*n+4] = total_strain(2,1);
       values[Strain][6*n+5] = total_strain(2,0);
     }
     
     //Internal Strain
     if (values.count(InternalStrain))
     {
       values[Strain][6*n] =   strain(0,0);
       values[Strain][6*n+1] = strain(1,1);
       values[Strain][6*n+2] = strain(2,2);
       values[Strain][6*n+3] = strain(1,0);
       values[Strain][6*n+4] = strain(2,1);
       values[Strain][6*n+5] = strain(2,0);
     }

     //Total strain in the crystal system
     if (values.count(StrainCrystal))
     {
       
       const Material* mat = mod.get_material();
       const RotatedCrystal&   cr = mat->get_rotated_crystal ();
       const Tensor2Gen& rotate = cr.RotMatrix;
       RealTensor crystal_strain = rotate.transpose() * (total_strain * rotate);

       values[StrainCrystal][6*n] =   crystal_strain(0,0);
       values[StrainCrystal][6*n+1] = crystal_strain(1,1);
       values[StrainCrystal][6*n+2] = crystal_strain(2,2);
       values[StrainCrystal][6*n+3] = crystal_strain(1,0);
       values[StrainCrystal][6*n+4] = crystal_strain(2,1);
       values[StrainCrystal][6*n+5] = crystal_strain(2,0);
     }


     if (values.count(Stress))
     {
       values[Stress][6*n] = total_stress(0,0);
       values[Stress][6*n+1] = total_stress(1,1);
       values[Stress][6*n+2] = total_stress(2,2);
       values[Stress][6*n+3] = total_stress(1,0);
       values[Stress][6*n+4] = total_stress(2,1);
       values[Stress][6*n+5] = total_stress(2,0);
     }

     //Internal Stress
     if (values.count(InternalStress))
     {

       values[InternalStress][6*n] = stress(0,0);
       values[InternalStress][6*n+1] = stress(1,1);
       values[InternalStress][6*n+2] = stress(2,2);
       values[InternalStress][6*n+3] = stress(1,0);
       values[InternalStress][6*n+4] = stress(2,1);
       values[InternalStress][6*n+5] = stress(2,0);

     }
  
     //----------Stress crystal
     if (values.count(StressCrystal))
     {
       
       const Material* mat = mod.get_material();
       const RotatedCrystal&   cr = mat->get_rotated_crystal ();
       const Tensor2Gen& rotate = cr.RotMatrix;
       RealTensor crystal_stress = rotate.transpose() * (total_stress * rotate);
       
       values[StressCrystal][6*n] =   crystal_stress(0,0);
       values[StressCrystal][6*n+1] = crystal_stress(1,1);
       values[StressCrystal][6*n+2] = crystal_stress(2,2);
       values[StressCrystal][6*n+3] = crystal_stress(1,0);
       values[StressCrystal][6*n+4] = crystal_stress(2,1);
       values[StressCrystal][6*n+5] = crystal_stress(2,0);
     }
     //----Force source--------
     if (values.count(ForceSource))
     { 
       values[ForceSource][3*n]   = force_source(0);
       values[ForceSource][3*n+1] = force_source(1);
       values[ForceSource][3*n+2] = force_source(2);
     }
     
     if (values.count(StressSource))
     {
       values[StressSource][6*n] = stress_source(0,0);
       values[StressSource][6*n+1] = stress_source(1,1);
       values[StressSource][6*n+2] = stress_source(2,2); 
       values[StressSource][6*n+3] = stress_source(1,0);
       values[StressSource][6*n+4] = stress_source(2,1);
       values[StressSource][6*n+5] = stress_source(2,0); 
     }
     
     if (values.count(StrainSource))
     {
       values[StrainSource][6*n] = strain_source(0,0);
       values[StrainSource][6*n+1] = strain_source(1,1);
       values[StrainSource][6*n+2] = strain_source(2,2); 
       values[StrainSource][6*n+3] = strain_source(1,0);
       values[StrainSource][6*n+4] = strain_source(2,1);
       values[StrainSource][6*n+5] = strain_source(2,0); 
     }

     if (values.count(EnergyDensity))
     {
       // E = 1/2*C_ijkl\eps_ij\eps_kl + \sigma_ext_ij\eps_ij
       double energy = 0;
       for (ID i = 0; i < 3; i++)
         for (ID j = 0; j < 3; j++)
           energy += (total_stress(i,j) + stress_source(i,j)) * total_strain(i,j);
       values[EnergyDensity][n] = 0.5e9 * energy;
     }

     //if (values.count(Energy))
     //{
     //  for (ID i = 0; i<dim; i++)
     //    for (ID j = 0; j<dim; j++)
     //      values[Energy][n] += 0.5 * (total_stress(i,j)) * (total_strain(i,j));
     //}

   }


    if (values.count(StrainCell))
    {
      values[StrainCell][0] = total_strain(0,0);
      values[StrainCell][1] = total_strain(1,1);
      values[StrainCell][2] = total_strain(2,2); 
      values[StrainCell][3] = total_strain(1,0);
      values[StrainCell][4] = total_strain(2,1);
      values[StrainCell][5] = total_strain(2,0); 
    }
      

}



Real 
Elasticity::compute_elastic_energy(void)
{


  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  DofMap& dof_map =  system.get_dof_map();
  FEType fe_type = dof_map.variable_type(uvar[0]);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim, FIFTH);
  fe->attach_quadrature_rule(&qrule);

  const vector<Point>& ref_points = qrule.get_points();
  const vector<Real>& JxW = fe->get_JxW();

  //Start assembling
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  Real energy(0);
  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    fe->reinit(elem);  

    std::vector<double> energy_p(qrule.n_points(),0.0);
    std::map<ID, std::vector<double> > values;
    values[Energy] = energy_p;
    get_solution_secure(elem, values, ref_points);
    
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
      energy += JxW[qp] * values[Energy][qp];

  }

  return energy;

}

RealTensor
Elasticity::get_stress(const Elem* elem, const Point& p)
{
  
  RealTensor stress(0);

    std::vector<double> stress_p(6,0.0);
    std::map<ID, std::vector<double> > values;
    values[Stress] = stress_p;
    
    std::vector<Point>  points(1);  points[0] = p;
    get_solution_secure(elem,values,points);

    stress(0,0) = values[Stress][0];
    stress(1,1) = values[Stress][1];
    stress(2,2) = values[Stress][2];
    stress(0,1) = values[Stress][3];
    stress(1,2) = values[Stress][4];
    stress(0,2) = values[Stress][5];
    stress(1,0) = stress(0,1);
    stress(2,0) = stress(0,2); 
    stress(2,1) = stress(1,2);

  return stress;
}


RealTensor
Elasticity::get_internal_stress(const Elem* elem, const Point& p)
{

  RealTensor stress(0);

    std::vector<double> stress_p(6,0.0);
    std::map<ID, std::vector<double> > values;
    values[InternalStress] = stress_p;

    std::vector<Point>  points(1);  points[0] = p;
    get_solution_secure(elem,values,points);

    stress(0,0) = values[InternalStress][0];
    stress(1,1) = values[InternalStress][1];
    stress(2,2) = values[InternalStress][2];
    stress(0,1) = values[InternalStress][3];
    stress(1,2) = values[InternalStress][4];
    stress(0,2) = values[InternalStress][5];
    stress(1,0) = stress(0,1);
    stress(2,0) = stress(0,2);
    stress(2,1) = stress(1,2);

  return stress;
}



void
Elasticity::do_assemble(EquationSystems& es, const std::string& system_name)
{


  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const NumericVector<Number>& solution = *sol;

  const MeshBase& mesh = get_mesh();
  //const unsigned int dim = mesh.mesh_dimension();
  ID dim = mesh.mesh_dimension();

  DofMap& dof_map =  system.get_dof_map();

  FEType fe_type = dof_map.variable_type(uvar[0]);

  SimulationEnvironment& si = get_environment();

  // the volume finite element
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  //AutoPtr<QBase> qrule(QBase::build(QTRAP, dim, FIFTH));

  QGauss qrule(dim, FIFTH);
  //QTrap qrule(dim);
  fe->attach_quadrature_rule(&qrule);

  //const vector<Point>& ref_points = qrule.get_points();

  const vector<Real>& JxW = fe->get_JxW();
  const vector<Point>& q_point = fe->get_xyz();
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();


  // the surface finite element
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  QGauss qface(dim - 1, SIXTH);
  //QTrap qface(dim - 1);
  fe_face->attach_quadrature_rule(&qface);

  const vector<Point>& ref_face_points = qface.get_points();

  const vector<Real>& JxW_face = fe_face->get_JxW();
  const vector<Point>& qface_point = fe_face->get_xyz();
  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  const vector<vector<RealGradient> >& dphi_face = fe->get_dphi();
  const vector<Point>& normal = fe_face->get_normals();


  //std::vector<std::vector<unsigned int> > dof_indices(3);
  //for (unsigned int i = 0; i< 3 ; i++)
  //  dof_map.dof_indices(elem, dof_indices[i],uvar[i]);


  //Initialize-----------------------------------------------
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector< DenseSubVector<Number>* > F(3);
  std::vector<std::vector< DenseSubMatrix<Number>* > > K(3);
  for (unsigned int i= 0;i<3; i++)
  {
    K[i].resize(3);
    for (unsigned int j= 0;j<3; j++)
      K[i][j] = new  DenseSubMatrix<Number> (Ke);

    F[i] = new DenseSubVector<Number> (Fe) ;
  }
  //----------------------------------------------------------
  std::vector< std::vector<unsigned int> > dof_indices_vec(3);
  std::vector<unsigned int> dof_indices;


  //Start assembling
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    
    dof_map.dof_indices (elem, dof_indices);
    const unsigned int n_dofs   = dof_indices.size();
    
    std::vector<unsigned int> n_dofs_vec(3);
    for (unsigned int i = 0;i <3; i++)
    {
      dof_map.dof_indices (elem, dof_indices_vec[i], uvar[i]);
      n_dofs_vec[i] = dof_indices_vec[i].size();
    }
    
    fe->reinit(elem);


    Ke.resize (n_dofs, n_dofs);
    Fe.resize (n_dofs);

    for (unsigned int i = 0;i <3; i++)
    {
      (F[i])->reposition(uvar[i] *n_dofs_vec[i],n_dofs_vec[i]);
      (F[i])->zero();
      for (unsigned int j = 0;j <3; j++)
      {
        (K[i][j])->reposition(uvar[i] * n_dofs_vec[i], uvar[j] * n_dofs_vec[i],  n_dofs_vec[i] , n_dofs_vec[j] );
	(K[i][j])->zero();
      }
    }

    //Bulk
    ElasticityModel& mod = *get_bulk_model<ElasticityModel>(elem);    
  
// commentato perche' centroide ha coordinate globali e produce errore
//    const RealTensor& internal_stress = get_internal_stress(elem, elem->centroid());

    // get the accumulated strain from shape deformation
    RealTensor accum_strain;
    _accumulated_strain[elem].get_tensor(accum_strain);

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
     
       mod.calculate(elem, qrule.qp(qp));
       const Tensor4DSym& C = mod.get_stiffness();
       const RealGradient& force =  mod.get_force_source();
       //const RealTensor& strain_source =  mod.get_strain_source();
       const RealTensor& stress_source =  mod.get_stress_source();

       RealTensor strain(mod.get_strain_source());
       strain += accum_strain;

       for (ID i = 0;i <3; i++)
         for (ID j = 0;j <3; j++)
           for (ID alpha=0; alpha<n_dofs_vec[i]; alpha++)
           {
             //Add internal stress SEEMS TO BE WRONG
             //if (shape_iteration > 0)
             //  strain(i,j) += 0.5 * ((solution)(dof_indices_vec[i][alpha]) * dphi[alpha][qp](j) + (solution)(dof_indices_vec[j][alpha]) * dphi[alpha][qp](i));

             for (ID beta=0; beta<n_dofs_vec[j]; beta++)
               (*(K[i][j]))(alpha,beta) += JxW[qp] * dphi[alpha][qp] * (get_subtensor(C,i,j) * dphi[beta][qp]);
           }

       RealTensor total_stress(stress_source);
       total_stress += C * strain;


        //Right Hand side
	for (unsigned int alpha=0; alpha<n_dofs_vec[0]; alpha++)
	{  

	  RealGradient tmp = total_stress * dphi[alpha][qp] + force *  phi[alpha][qp];
	  for (ID i = 0;i <3; i++)
	    (*(F[i]))(alpha) -= JxW[qp] * tmp(i);

	}
      
    }//End QP



    //Boundary conditions
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      const ElementSide elside(elem->top_parent(),s);

      if (si.is_outer_boundary(elside))
      {
          ElasticityBoundaryModel*  boundary_mod =
              get_interface_model<ElasticityBoundaryModel>(elem, s);

          fe_face->reinit(elem, s);

          RealTensor H(0);
          RealGradient R(0);

          if (boundary_mod != NULL)
          {
            boundary_mod->set_normal(normal[0]);
            boundary_mod->calculate(elem, s,elem->centroid());
            double b(0);
            boundary_mod->get_coefficients(H, b, R);

            //if (boundary_mod->get_type() == "ebnd_extended")
            //{
            //  mod.calculate(elem, qface.qp(0));
            //  const RealTensor& strain =  mod.get_strain_source();
            //  const RealTensor& stress =  mod.get_stress_source();
            //  const Tensor4DSym& C = mod.get_stiffness();
            //  R = ( C * strain) * normal[0];
            //}

            if ((b < 1e-10) && (b >= 0)) b = 1e-20;
            else if ((b > -1e-10) && (b<= 0)) b = -1e-20;
            H /= b;
            R /= b;
          }

          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
            for (unsigned int alpha=0; alpha<n_dofs_vec[0]; alpha++)
              for (unsigned int i =0; i<3; i++)
              {
                (*(F[i]))(alpha) +=  JxW_face[qp] * phi_face[alpha][qp] * R(i);

                for (unsigned int j =0; j<3; j++)
                  for (unsigned int beta=0; beta<n_dofs_vec[0]; beta++)
                    (*(K[i][j]))(alpha,beta) += JxW_face[qp] * H(i,j) * phi_face[alpha][qp] * phi_face[beta][qp];

              }
      }
    }

    //-------------------------
    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices);
  }

  //system.matrix->close();
  //system.matrix->print_matlab("K.m");
  //system.rhs->close();
  //system.rhs->print_matlab("F.m");

}

void
Elasticity::apply_shape_deformation()
{
  
  //Atomistics deformation-----
  vector<AtomisticStructure*> atom_structures;

  Device& device = get_environment().get_device();
  device.get_atomistic_structures(myopt.structure_to_be_strained,atom_structures);
   
  TiberLinearSystem* system = &get_equation_system<TiberLinearSystem>();
  const unsigned int system_number = system->number();
  
  // TODO CHECK THIS! It seems that this was terribly wrong and worked only
  // for one deformation step
  //const NumericVector<Number>& solution = *sol;
  const NumericVector<Number>& solution = system->get_solution_vector();
  const unsigned int dim = get_mesh().mesh_dimension();
  const DofMap& dof_map = system->get_dof_map();
  
  std::vector<std::vector<unsigned int> > dof_indices(3);
  
  FEType fe_type = system->variable_type(uvar[0]);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  const std::vector<std::vector<Real> >& phi = fe->get_phi();


  for (unsigned int ns = 0; ns < atom_structures.size(); ns++)
  {
    std::vector<Atom>& structure = atom_structures[ns]->get_structure_atoms();
    std::vector<std::vector<unsigned int> > bond_map = atom_structures[ns]->get_bond_map();
    double scale = atom_structures[ns]->get_scale();
   
    for (unsigned int na = 0; na < structure.size(); na++)
    {
      Point old_pos(structure[na].get_position() / scale);
       
      // will be the point in the reference coordinates, but is used also to
      // pass the real space coordinates to inverse_map
      vector<Point> p_ref(1, old_pos);

      const Elem* elem = structure[na].get_elem();
      if (elem == NULL)
      {
        // here we may arrive when the atom is a passivation atom
        // we simply take the displacement at the neighbour atom
        if (bond_map[na].size() == 0)
          throw RuntimeException("Found atom with no neighbours!");

        unsigned int neighbor = bond_map[na][0];

        elem = structure[neighbor].get_elem();

        if (elem == NULL)
          throw RuntimeException("Found atom with no neighbour inside the mesh!");

        p_ref[0] = structure[neighbor].get_position() / scale;
      }

      if (!get_environment().contains_element(elem))
      {
        // it must be in a quantum contact, otherwise the problem is badly set up
        Device::quantum_contact_iterator it(device.quantum_contacts_begin());
        Device::quantum_contact_iterator end(device.quantum_contacts_end());
        for ( ; it != end; ++it)
        {
          QuantumContact* qct = it->second;
          pair<const Elem*, Point> proj = qct->project_on_boundary(elem, p_ref[0]);
          if (proj.first != NULL)
          {
            elem = proj.first;
            p_ref[0] = proj.second;
            break;
          }
        }
      }

      p_ref[0] = FEInterface::inverse_map(get_mesh().mesh_dimension(), FEType(), elem, p_ref[0]);

      fe->reinit(elem, &p_ref);

      for (unsigned int i = 0; i< 3 ; i++)
	dof_map.dof_indices(elem, dof_indices[i], uvar[i]);
  
      Point displ(0);
      for (unsigned int i = 0;i<3; i ++)
	for (unsigned int alpha = 0; alpha < dof_indices[i].size(); alpha++)
	  displ(i) += (solution)(dof_indices[i][alpha]) * phi[alpha][0];
        
      displ /= get_scaling().get_calc_mesh_units();

      old_pos += displ;
      old_pos *= scale;
      
      structure[na].set_position(old_pos);
      
    }
    atom_structures[ns]->print_structure("strained.xyz");
  }
  
  //
  // we may have to move some quantum contact regions
  //
  Device::quantum_contact_iterator it(device.quantum_contacts_begin());
  Device::quantum_contact_iterator end(device.quantum_contacts_end());
  for ( ; it != end; ++it)
  {
    QuantumContact* qct = it->second;
    const set<ID>& regions = qct->get_region_ids();
    bool touches = false;
    for (set<ID>::iterator sit(regions.begin()); sit != regions.end(); ++sit)
      touches |= get_environment().contains_region(*sit);

    if (touches)
    {
      vector<Point> proj_point(1);

      set<const Node*> visited;
      QuantumContact::ContactElemIterator cel(qct->contact_elements_begin());
      QuantumContact::ContactElemIterator cend(qct->contact_elements_end());
      for ( ; cel != cend; ++cel)
      {
        const Elem* elem = cel->first;
        for (unsigned int n = 0; n < elem->n_nodes(); ++n)
        {
          Node* node = elem->get_node(n);
          if ((node->n_dofs(system_number, uvar[0]) != 0) ||
              visited.count(node))
          {
            continue;
          }

          visited.insert(node);
          
          pair<const Elem*, Point> proj = qct->project_on_boundary(elem, elem->point(n));
          const Elem* proj_elem = proj.first;
          proj_point[0] = proj.second;

          proj_point[0] = FEInterface::inverse_map(get_mesh().mesh_dimension(),
              FEType(), proj_elem, proj_point[0]);

          fe->reinit(proj_elem, &proj_point);

          for (unsigned int i = 0; i < 3 ; i++)
            dof_map.dof_indices(proj_elem, dof_indices[i], uvar[i]);

          Point displ(0);
          for (unsigned int i = 0; i < 3; i++)
            for (unsigned int alpha = 0; alpha < dof_indices[i].size(); alpha++)
              displ(i) += (solution)(dof_indices[i][alpha]) * phi[alpha][0];

          node->add(displ / get_scaling().get_calc_mesh_units());
        }
      }

    }
  }


  //
  // now move the mesh
  //
  const MeshBase& mesh = get_mesh();
  MeshBase::const_node_iterator  nd  = mesh.active_nodes_begin();
  const MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();

  for ( ;  nd != nd_end ; ++nd)
  {
    Node* node = *nd;

    // If there are no DOFs, it's not a node of the simulation domain
    if (node->n_dofs(system_number, uvar[0]) == 0)
      continue;

    Point pos;
    for(unsigned int i = 0; i<dim; i++)
      pos(i) = (*node)(i);

    for (unsigned int i = 0; i < dim; i++)
    {
      const unsigned int  n_dof = node->dof_number(system_number,uvar[i],0);
      pos(i) += (solution)(n_dof) / get_scaling().get_calc_mesh_units();
    }

    *node = pos;
  }

  /*
  it = device.quantum_contacts_begin();
  for ( ; it != end; ++it)
    it->second->activate_elements();

  DataOutput* dop = DataOutput::create("vtk");
  dop->set_mesh(mesh);
  dop->set_filename("deformed");
  dop->set_output_directory("./");
  dop->write(true);

  it = device.quantum_contacts_begin();
  for ( ; it != end; ++it)
    it->second->inactivate_elements();
  */
}

void
Elasticity::restore_shape()
{


  //MeshDeformation----
  TiberLinearSystem* system = &get_equation_system<TiberLinearSystem>();
  
  const NumericVector<Number>& solution = *sol;

  const unsigned int system_number = system->number();
  const MeshBase& mesh = get_mesh();
  ID  dim = get_mesh().mesh_dimension();
  MeshBase::const_node_iterator  nd  = mesh.active_nodes_begin();
  const MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();
  
  for ( ;  nd != nd_end ; ++nd)
  {
    Node* node = *nd;
    Point pos;
    for(unsigned int i = 0; i<3; i++)
      pos(i) = (*node)(i);
    
    for (unsigned int i = 0; i < 3; i++)
    {
      const unsigned int  n_dof = node->dof_number(system_number,uvar[i],0);
      pos(i) -= (solution)(n_dof)/ get_scaling().get_calc_mesh_units(); 
    }
    
    *node = pos;
  }

}

RealTensor
Elasticity::get_subtensor(const Tensor4DSym& C_calc,unsigned int i,unsigned  int j)
{


  const Tensor4DSym&  C1 = C_calc;



  RealTensor a;
  for (unsigned int k = 0; k<3; k ++)
  {
    for (unsigned int m = 0; m<3; m ++)
    {
      const double  p =  C1(i+1,k+1,j+1,m+1);
      a(k,m) = p;

    }
  }

  return a;

}



//Surface Integration
//    for (unsigned int s = 0; s < elem->n_sides(); s++)
//    {
//      ElasticityBoundaryModel*  boundary_mod =
//        get_interface_model<ElasticityBoundaryModel>(elem, s);
//
//          fe_face->reinit(elem, s);
//
//          if (boundary_mod == NULL)
//          {
//
//            for (unsigned int qp = 0; qp < qface.n_points(); qp++)
//            {
//              for (unsigned int alpha=0; alpha<n_dofs_vec[0]; alpha++)
//              {
//                RealGradient tmp = get_stress(elem,qface_point[qp]) * normal[qp];
//                for (ID i = 0;i <3; i++)
//                  (*(F[i]))(alpha) += JxW_face[qp] * phi_face[alpha][qp] * tmp(i) ;
//              }
//            }
//          }
//    }
//
