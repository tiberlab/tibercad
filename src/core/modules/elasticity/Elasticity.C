// $Id$

#include "Elasticity.h"
#include "ElasticityModel.h"
#include "ElasticityBoundaryModel.h"
#include "TiberLinearSystem.h"
#include "Messages.h"
#include "equation_systems.h"
#include "dof_map.h"
#include "quadrature_gauss.h"
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

// This is needed in order to create the shared module library
// The first string is the class name of the object to be created,
// the second one is the name of the module as it should be referred
// in the input file (the Makefile defines MODULE_NAME, which can be used here).
TIBER_MODULE(Elasticity, MODULE_NAME)


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
  
  // ID  dim = get_mesh().mesh_dimension();
  uvar.resize(3);
  //switch (dim)
  // {
  //case 3:
    system.add_variable("uz", FIRST);
    uvar[2] =  system.variable_number("uz");
    // case 2:
    system.add_variable("uy", FIRST);
    uvar[1] =  system.variable_number("uy");
    //default:
    system.add_variable("ux", FIRST);
    uvar[0] =  system.variable_number("ux");
    // }

  system.attach_assemble_function(assemble);
  system.init();

}


void
Elasticity::parse_options(void)
{


  const ModelOptions& opt = get_options();

  myopt.shape_error = opt.get_option("shape_error",1e-2);
  myopt.shape_iterations = opt.get_option("shape_iterations",1);
  myopt.deformation = opt.get_option("do_deformation",true);
  myopt.magnification = opt.get_option("magnification",1);



}


void
Elasticity::do_setup_solution_variables(void)
{
  // we declare our solution variables

  declare_solution(Strain, TENSOR, NODES, "");
  declare_solution(StrainCell, TENSOR, CELL, "");
  declare_solution(StrainCrystal, TENSOR, NODES, "");
  declare_solution(Energy, REAL, NODES, "Joule");
  declare_solution(Stress, TENSOR, NODES, "GPa");
  declare_solution(Displacement, VECTOR, NODES, "m");
  declare_solution(StrainSource, TENSOR, NODES, "");
  declare_solution(StressSource, TENSOR, NODES, "GPa");
  declare_solution(ForceSource, TENSOR, NODES, "N/m3 ");

}


void
Elasticity::do_solve(void)
{
  _this = this;

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  sol =  (system.solution)->clone();
  sol->zero();

  
  iter = 0;
  double error = 0.0;
  //double tot_norm = 0.0;
  double energy = 0.0;

  //if ((SimulationOptions::verbose() > 2) &&  myopt.shape_iterations>1)
  // {
  //  cout<<"| Iteration | Elastic energy [Joule] | norm(u) [m] | Energy error | Displacement error |"<<endl;
  //  cout<<"-----------------------------------------------------------------------------------------------"<<endl;
  // }

  do {
    iter +=1;
    //system.solution->zero();
    system.solve();
    sol->add(1.0,*(system.solution));
    
    //double norm = (system.solution)->l2_norm();
    //tot_norm +=norm;
    //error = norm/tot_norm;
  
    double new_energy = compute_elastic_energy();
    error = abs(new_energy - energy)/new_energy;
    energy = new_energy;

    if ((SimulationOptions::verbose() > 2))
    {
      if (myopt.shape_iterations > 1)
	cout<<"Iter: "<<iter;
	  
      cout<<" Elastic energy: "<<energy<<" Joule";
    
      if ( iter>1 )
	cout<<" Relative error: "<< error;

      cout<<endl;
    }

    if (myopt.deformation)
      apply_shape_deformation();
    
    //  iter +=1;
  } while (error > myopt.shape_error & iter < myopt.shape_iterations);

  // if (!myopt.deformation)
  // restore_shape();
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
 
  return ElasticityModel::create(options);
}



PhysicalModel*
Elasticity::create_boundary_model(const ModelOptions& options,
    const Material* material_A, const Material* material_B) const
{
  return ElasticityBoundaryModel::create(options);
}



void
Elasticity::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
				const std::vector<Point>& p)
{
 
   unsigned int np = p.size();

   TiberLinearSystem* system = &get_equation_system<TiberLinearSystem>();

   //const NumericVector<Number>& solution = system->get_solution_vector();
   const NumericVector<Number>& solution = *sol;
   const unsigned int dim = get_mesh().mesh_dimension();

   const DofMap& dof_map = system->get_dof_map();

   std::vector<std::vector<unsigned int> > dof_indices(3);
   for (unsigned int i = 0; i< 3 ; i++)
     dof_map.dof_indices(elem, dof_indices[i],uvar[i]);

   FEType fe_type = system->variable_type(uvar[0]);
   AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
   const std::vector<std::vector<Real> >& phi = fe->get_phi();
   const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  
   fe->reinit(elem, &p);
  
   ElasticityModel& mod = *get_bulk_model<ElasticityModel>(elem);
   
   RealTensor total_strain(0);
   RealTensor total_stress(0);

   for (ID n = 0; n < np; n++)
   {
     mod.calculate(elem,p[n]);   
     
     const Tensor4DSym& C = mod.get_stiffness();
     const RealGradient& force_source =  mod.get_force_source();
     const RealTensor& strain_source  =  mod.get_strain_source();
     const RealTensor& stress_source  =  mod.get_stress_source();
    
     //----Displacemet--------
     if (values.count(Displacement))
     {
       //----Displacemet-------- 
       RealGradient u(0);
     for (unsigned int i = 0;i<3; i ++)
	 for (unsigned int alpha = 0; alpha<dof_indices[i].size() ;alpha ++)
	   u(i) +=(solution)(dof_indices[i][alpha]) * phi[alpha][n];
       
       values[Displacement][3*n]   = u(0);
       values[Displacement][3*n+1] = u(1);
       values[Displacement][3*n+2] = u(2);
     }

     RealTensor strain(0);
     if (values.count(Strain) || values.count(StrainCrystal) ||
         values.count(Stress) || values.count(Energy) )
     {
       //------Strain-------------------------
      
       for (unsigned int i = 0;i<3; i ++)
	 for (unsigned int j = 0;j<=i; j ++)
	 {
	   
	   double der1 = 0;
	   for (unsigned int alpha = 0; alpha<dof_indices[i].size() ;alpha ++)
	     der1 += (solution)(dof_indices[i][alpha]) * dphi[alpha][n](j);
	   
	   double der2 = 0;
	   for (unsigned int alpha = 0; alpha<dof_indices[j].size() ;alpha ++)
	     der2 += (solution)(dof_indices[j][alpha]) * dphi[alpha][n](i);
	   
	   strain(i,j) = 0.5 * (der1 + der2);
	   
	 }


       total_strain = strain + strain_source;
       if (values.count(Strain))
       {
	 values[Strain][6*n] =   total_strain(0,0);
	 values[Strain][6*n+1] = total_strain(1,1);
	 values[Strain][6*n+2] = total_strain(2,2); 
	 values[Strain][6*n+3] = total_strain(1,0);
	 values[Strain][6*n+4] = total_strain(2,1);
	 values[Strain][6*n+5] = total_strain(2,0); 
       }

       if (values.count(StrainCrystal))
       {

	 Material* mat = mod.get_material();
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

     }
     
     //Stress
     if (values.count(Stress) || values.count(Energy))
     {
       
       total_stress = C * total_strain + stress_source;

       if (values.count(Stress))
	 {
	   values[Stress][6*n] = total_stress(0,0);
	   values[Stress][6*n+1] = total_stress(1,1);
	   values[Stress][6*n+2] = total_stress(2,2); 
	   values[Stress][6*n+3] = total_stress(1,0);
	   values[Stress][6*n+4] = total_stress(2,1);
	   values[Stress][6*n+5] = total_stress(2,0); 
	 }
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

     if (values.count(Energy))
     {
       for (ID i = 0; i<dim; i++)
	 for (ID j = 0; j<dim; j++)
	   values[Energy][n] += 0.5 * (total_stress(i,j)) * (total_strain(i,j));
     }

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
    get_solution_secure(elem,values,ref_points);
    
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
      energy += JxW[qp] * values[Energy][qp];

  }

  return energy;

}

RealTensor
Elasticity::get_stress(const Elem* elem, const Point& p)
{
  
  RealTensor stress(0);

  if (iter>1)
  {
   
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

   
  }

  return stress;
}




void
Elasticity::do_assemble(EquationSystems& es, const std::string& system_name)
{


  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

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
  fe->attach_quadrature_rule(&qrule);

  const vector<Point>& ref_points = qrule.get_points();

  const vector<Real>& JxW = fe->get_JxW();
  const vector<Point>& q_point = fe->get_xyz();
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();


  // the surface finite element
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  QGauss qface(dim - 1, SIXTH);
  fe_face->attach_quadrature_rule(&qface);

  const vector<Point>& ref_face_points = qface.get_points();

  const vector<Real>& JxW_face = fe_face->get_JxW();
  const vector<Point>& qface_point = fe_face->get_xyz();
  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  const vector<vector<RealGradient> >& dphi_face = fe->get_dphi();
  const vector<Point>& normal = fe_face->get_normals();


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
  
    
    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
     

       mod.calculate(elem,q_point[qp]);   
       const Tensor4DSym& C = mod.get_stiffness();
       const RealGradient& force =  mod.get_force_source();
       const RealTensor& strain =  mod.get_strain_source();
       const RealTensor& stress =  mod.get_stress_source();


      for (ID i = 0;i <3; i++)
	for (ID j = 0;j <3; j++)
	  for (ID alpha=0; alpha<n_dofs_vec[i]; alpha++)
	    for (ID beta=0; beta<n_dofs_vec[j]; beta++)
	      (*(K[i][j]))(alpha,beta) += JxW[qp] * dphi[alpha][qp] * (get_subtensor(C,i,j) * dphi[beta][qp]);
	    


      //Body force of previous iteration
      for (unsigned int alpha=0; alpha<n_dofs_vec[0]; alpha++)
      {  
	RealGradient tmp =   get_stress(elem,ref_points[qp]) * dphi[alpha][qp];
	for (unsigned int i = 0;i <3; i++)
	  (*(F[i]))(alpha) -= JxW[qp] * tmp(i);
      }

      // Total source
      for (unsigned int alpha=0; alpha<n_dofs_vec[0]; alpha++)
      {  
	RealGradient tmp = force *  phi[alpha][qp] + (stress  + (C * strain)) * dphi[alpha][qp];
	for (unsigned int i = 0;i <3; i++)
	  (*(F[i]))(alpha) -= JxW[qp] * tmp(i);
      }

      
    }//End QP



    // bool newb = true;
    //Boundary Conditions
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      ElasticityBoundaryModel*  boundary_mod =
	get_interface_model<ElasticityBoundaryModel>(elem, s);

	if (boundary_mod != NULL)
	{
	  fe_face->reinit(elem, s);
	  
	  for (unsigned int qp = 0; qp < qface.n_points(); qp++)
	  {
	    boundary_mod->calculate(elem, s, qface_point[qp]);
	    
	    RealTensor H(0);
	    RealGradient R(0);
	    double A(0);

            boundary_mod->set_normal(normal[qp]);
	    boundary_mod->get_coefficients(H, R);
	    
	    for (unsigned int alpha=0; alpha<n_dofs_vec[0]; alpha++)
	    {
	      for (unsigned int i =0; i<3; i++)
	      {
		
		(*(F[i]))(alpha) +=  JxW_face[qp] * R(i) * phi_face[alpha][qp];
		
		for (unsigned int j =0; j<3; j++)
		  for (unsigned int beta=0; beta<n_dofs_vec[0]; beta++)
		    (*(K[i][j]))(alpha,beta) += JxW_face[qp] * H(i,j) * phi_face[alpha][qp] * phi_face[beta][qp];
		
	      }
	    }
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

//   as = get_atomistic_structure();

//   Device::StructureMap::iterator str_start = as.begin();
//   const  Device::StructureMap::iterator str_end = as.end();

//   for ( ;  str_start != str_end ; ++str_start)
//   {
//     double scale =  (str_start->second)->get_scale();
//     std::vector< Atom >& structure =  (str_start->second)->get_structure_atoms();

//     for (unsigned int na = 0; na < structure.size(); na++)
//     {

//       std::vector<Point> old_pos(1);
//       old_pos[0](0) = structure[na].get_position()(1) / scale;
//       old_pos[0](1) = structure[na].get_position()(2) / scale;
//       old_pos[0](2) = structure[na].get_position()(3) / scale;

//       std::vector<Tensor1> disp(1);
//       get_displacement(structure[na].get_elem(),old_pos,disp);


//        disp[0] /= get_scaling().get_calc_mesh_units(); //To mesh units

//       Tensor1 new_pos(0);
//       new_pos(1) = disp[0](1) + old_pos[0](0);
//       new_pos(2) = disp[0](2) + old_pos[0](1);
//       new_pos(3) = disp[0](3) + old_pos[0](2);

//       new_pos *=scale;
//       structure[na].set_position(new_pos);

//     }
//     (str_start->second)->print_structure("strained.xyz");
//   }


  //MeshDeformation----
  TiberLinearSystem* system = &get_equation_system<TiberLinearSystem>();
  
  const NumericVector<Number>& solution = system->get_solution_vector();
  
  const unsigned int system_number = system->number();
  const MeshBase& mesh = get_mesh();
  ID  dim = get_mesh().mesh_dimension();
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
      pos(i) += (solution)(n_dof)/ get_scaling().get_calc_mesh_units(); 
    }
    
    *node = pos;
  }

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
//    //Adjusting the point lieing on nodes.
//    vector<Point> global_point(elem->n_nodes());
//    for (ID node = 0; node < elem->n_nodes(); node++)
//      global_point[node] =  elem->point(node);
   
     
//    std::vector<Point> local_point(elem->n_nodes()+1);
//    FEInterface::inverse_map(dim, fe_type, elem, global_point,local_point);
   
//    std::vector<Point> point_new(np);
//    for (ID n = 0; n < np; n++)
//    {
//      bool is_on_node = false;
//      for (ID node = 0; node < elem->n_nodes(); node++)
//      {
//        Point p_node = elem->node(node);
       
//        Tensor1 diff(0);
//        for (ID d = 0;d <3; d++) 
// 	 diff(d+1) = local_point[node](d) - p[n](d);
       
//        double dist = norm(diff);
      
//        if (dist < 1E-8 && dist > -1E-8 )
//          is_on_node = true;
//      }
     
   
//      if (is_on_node)
//        point_new[n] =   p[n] * (1.0);
//      else
//        point_new[n] =   p[n];

//      //cout<< point_new[n]<<endl;
//    }
//    //------------------------------------------------
   //Surface contribution of internal body force
 //    for (unsigned int side = 0; side<elem->n_sides(); side++)
//     {
//       const ElementSide elside(elem->top_parent(), side);
//       if (si.is_on_boundary(elside))
//       {
//         bool do_int = false;

// 	ElasticityBoundaryModel*  boundary_mod =
// 	  get_interface_model<ElasticityBoundaryModel>(elem, side);
	
// 	if (boundary_mod != NULL)
// 	  if  (boundary_mod->get_type() == ElasticityBoundaryModel::SurfaceForce)
// 	    do_int = true;
	
// 	if (do_int)
// 	{
// 	  fe_face->reinit(elem,side);
// 	  for (unsigned int qp=0; qp<qface.n_points(); qp++)
// 	  {
// 	    RealGradient tmp =   get_stress(elem,ref_face_points[qp]) * normal[qp];
// 	    for (unsigned int alpha=0; alpha<n_dofs_vec[0]; alpha++)
// 	      for (ID i = 0;i <dim; i++)
// 		(*(F[i]))(alpha) +=  JxW_face[qp] * phi_face[alpha][qp] * tmp(i);
// 	  }
// 	}
	
//       }
//     }



	//------------------------
	//}
	// else
	// {

	
// 	if (boundary_mod != NULL)
// 	{
	  
// 	  if (boundary_mod->get_type() == ElasticityBoundaryModel::SurfaceForce)
// 	  {
// 	    fe_face->reinit(elem, s);
	    
// 	    RealTensor H(0);
// 	    RealGradient R(0);
// 	    boundary_mod->get_coefficients(H,R);
	    
	    
// 	    for (unsigned int qp=0; qp<qface.n_points(); qp++)
// 	      for (unsigned int k =0; k<3; k++)
// 		for (unsigned int alpha=0; alpha<n_dofs_vec[k]; alpha++)
// 		  (*(F[k]))(alpha) +=  JxW_face[qp] * R(k) * phi_face[alpha][qp];
	    
// 	  }
// 	  else     
// 	  {
	    
// 	    //  //--------------
// 	    RealTensor H(0);
// 	    RealGradient R(0);
// 	    boundary_mod->get_coefficients(H,R);
	    
// 	    for (unsigned int i = 0;i <3; i++)
// 	    {
// 	      for (unsigned int alpha=0; alpha<n_dofs_vec[i]; alpha++)
// 	      {
// 		if (elem->is_node_on_side(alpha,s))
// 		{
// 		  for (unsigned int j = 0;j <3; j++)
// 		    for (unsigned int beta=0; beta<n_dofs_vec[j]; beta++)
// 		      (*(K[i][j]))(alpha,beta) = 0.0;
		  
// 		  for (unsigned int j = 0;j <3; j++)
// 		    (*(K[i][j]))(alpha,alpha) = H(i,j);
		  
// 		  (*(F[i]))(alpha) = R(i);
		  
// 		}// if (elem->is_node_on_side(alpha,side))
// 	      }
// 	    }
// 	  //-----------------
// 	  }
//	}
	
	//     }

	//   }
    
