// $Id$

#include "MicroHeatBalance.h"
#include "BoundaryProperties.h"
#include "TiberLinearSystem.h"
#include "PhysicalModel.h"
#include "mesh.h"
#include "dof_map.h"
#include "equation_systems.h"
#include "fe.h"
#include "fe_base.h"
#include "elem.h"
#include "side.h"
#include "quadrature_gauss.h"

 // Define useful datatypes for finite element
 // matrix and vector components.
 #include "sparse_matrix.h"
 #include "numeric_vector.h"
 #include "dense_matrix.h"
 #include "dense_vector.h"

 // Define the DofMap, which handles degree of freedom
 // indexing.
 #include "dof_map.h"

 #include "fe_interface.h"
#include "dense_submatrix.h"
 #include "dense_subvector.h"

#include "HeatModel.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Boundary.h"
#include "Reservoir.h"
#include "FluxContact.h"
#include "SimulationOptions.h"
#include "SimulationEnvironment.h"
#include "ThermalSurfaceResistance.h"
#include "ThermalSurfaceConductance.h"
#include "MesoThermalContact.h"
#include "FourierBTE.h"
#include "Specular.h"
#include "Diffusive.h"
#include "SimulationOptions.h"
#include<iostream>
#include<fstream>
#include<sstream>


using namespace std;
MicroHeatBalance* MicroHeatBalance::static_this;
Device* MicroHeatBalance::_device;


//-----------------------------------------------------------------//


void MicroHeatBalance::parse_options( )
{

  const ModelOptions& opts = SimulationInterface::get_options();

  //myopts.integration_order = static_cast<libMeshEnums::Order>(
  //   opts.get_option("integration_order", 5));

  myopts.work_units = opts.get_option("Work_length_units", 1e-2);

  // myopts.scale = opts.get_option("scale", 1.0);

  myopts.max_error =  opts.get_option("max_error",1e-3);

  //  myopts.DG =  opts.get_option("DG",0);

  myopts.max_iter =  opts.get_option("max_iter",1);

  myopts.t_0 = opts.get_option("time_scale",1.0);
  myopts.s_0 = opts.get_option("space_scale",1.0);

  //myopts.ref_temp =  opts.get_option("reference_temperature",SimulationOptions::temperature);

  AngInt.theta_slices = opts.get_option("theta_slices",0);
  AngInt.phi_slices =    opts.get_option("phi_slices",0);
  AngInt.dim = dim;

  myopts.diffusive = opts.get_option("diffusive_walls",true); ;


  //Custom Dir
  myopts.alternative = opts.get_option("custom_dir",0);
  int n_dir = opts.get_option("n_dir",1);


  myopts.cd.resize(n_dir);  
  //Custom Direction
  for (ID n = 0; n<n_dir; n++)
  {

    std::string str = static_cast<ostringstream*>( &(ostringstream() << n+1) )->str();
    std::string opt_vec = "dir_" + str;
    
    std::vector<double> cd(3);
    cd[0] = 0;
    cd[1] = 0;
    cd[2] = 1;
    opts.get_option(opt_vec,cd);
    
    myopts.cd[n](0) = cd[0];
    myopts.cd[n](1) = cd[1];
    myopts.cd[n](2) = cd[2];
    
  }  
  


  
  myopts.equilibrium_energy = opts.get_option("equilibrium_energy",SimulationOptions::temperature);
  
  myopts.first_guess = opts.get_option("first_guess","fourier");


}

void MicroHeatBalance::do_init( )
{

  const ModelOptions& sim_opt = get_options();
   SimulationEnvironment& si = get_environment();
 
  _device = &( si.get_device() );
  mesh = &(_device->get_mesh());
  dim = mesh->mesh_dimension();

  parse_options();



  double mesh_units = get_scaling().get_calc_mesh_units() / sim_opt.get_option("Work_length_units", 1e-2 / myopts.s_0);

   get_scaling().set_calc_mesh_units(mesh_units);


  my_system = TiberLinearSystem::create(get_equation_systems(),
      get_equation_system_name(), get_solver_options());

  //if (myopts.DG == 1)
  my_system->add_variable("T",CONSTANT,MONOMIAL);
  //else
    //my_system->add_variable("T",FIRST,LAGRANGE);

   
   // Insert the pointer to function that LibMesh library has to use
    // if (myopts.DG == 1)
    my_system->attach_assemble_function(assemble_heat_matrix);
    //else
    //my_system->attach_assemble_function(assemble_heat_matrix);

  my_system->set_options(get_solver_options());
   // Initialize the data structures for the equation system.
  my_system->init();

  //Inizialize the solution to temperature of simulation options
  gray_sys_number = my_system->number();

  my_system->solution->zero();
  my_system->solution->close();

  heat_legend = "Wq";

  JQ_var.insert(JQX);
  JQ_var.insert(JQY);
  JQ_var.insert(JQZ);

  e0_var.insert(E0);

  if (myopts.alternative)
    AngInt.compute_custom_direction(myopts.cd);
  else
    AngInt.compute_directions();

  if  (SimulationOptions::verbose() > 2)
    AngInt.print_info(); 

  //Initialize thermal flux
  thermal_flux.resize(dim);
  for (ID i = 0; i<dim; i++ )
  { 
    thermal_flux[i] = (my_system->solution)->clone().release();
    thermal_flux[i]->zero();
    thermal_flux[i]->close();
  }

  //Initialize direction solution
  sol_dir.resize(AngInt.n_slices);
  for (ID k = 0; k<AngInt.n_slices ; k++ )
  {
    sol_dir[k] = (my_system->solution)->clone().release();
    sol_dir[k]->zero();
  }
  
  

  //Initialize equilibrium_energy
  equilibrium_energy = (my_system->solution)->clone().release();
  equilibrium_energy->add(SimulationOptions::temperature); //Just put some dummy value
  equilibrium_energy->close();
  energy_norm = equilibrium_energy->l2_norm();


  //DOMAIN PARTITIONING---Prototype for multiscale
  


//  MeshBase::const_element_iterator       it    = mesh->active_elements_begin();
//  const MeshBase::const_element_iterator it_end = mesh->active_elements_end(); 
//  for ( ; it != it_end; ++it)
//  {
//    const Elem* elem = *it;

//    HeatModel& mod = *get_bulk_model<HeatModel>(elem);

      


    
//   }
   





 

  
}
//-------------------------------------------------------------------------------//
void  MicroHeatBalance::do_solve()
{
  
  // std::cout<<"Space scale: "<<myopts.s_0<<std::endl;
  //std::cout<<"Time scale: "<<myopts.t_0<<std::endl;

  static_this = this;
  SimulationEnvironment& si = get_environment();

 
  if (~is_solved())
    compute_fourier_solution();

  //------------------------------------

  //Fill the directional results
  for (ID k = 0; k<AngInt.n_slices ; k++ )
  {
    sol_dir[k]->zero();
    sol_dir[k]->add(*equilibrium_energy);
  }

  //WRITE BOUNDARY DATA
  MeshBase::const_element_iterator       it    = mesh->active_elements_begin();
  const MeshBase::const_element_iterator it_end = mesh->active_elements_end(); 
  for ( ; it != it_end; ++it)
  {
   
    const Elem* elem = *it;
    
    ID dof = elem->dof_number(gray_sys_number,0,0);
    double eq_value = (*equilibrium_energy)(dof);

    for (unsigned int ns = 0; ns < elem->n_sides(); ns++)
    { 
      const ElementSide elside(elem->top_parent(), ns); 
      std::vector<double> sol_d(AngInt.n_slices,eq_value);
      SD[elside] = sol_d;
    }
    
  }
  //-------------------------------------------

 
  
  double err = 0.0;
  iter = 0;
  if  (myopts.max_iter > 0)
  {
    
    do {
   
    for (unsigned int k = 0; k<AngInt.n_slices; k++ )
    { 
      
      vec_spec = AngInt.spec[k];
      d_omega = AngInt.d_omega[k];
      IntDir = AngInt.directions[k];
      dir = AngInt.dir[k];

      if  (SimulationOptions::verbose() > 2)
	AngInt.print_info(k);

      (my_system->solution)->zero();
      my_system->solve();

      sol_dir[k]->zero();
      sol_dir[k]->add(*(my_system->solution));
     
      
    }
    
    //-----UPDATE BOUNDARY VALUES---------
    SideData::iterator it(SD.begin()); 
    SideData::const_iterator it_end(SD.end()); 
    
    for ( ; it != it_end; ++it)
    {
      for (ID k = 0; k<AngInt.n_slices ; k++ ) 
      {
	ID dof = ((it->first).elem())->dof_number(gray_sys_number,0,0);
      	double eq_value =  (*sol_dir[k])(dof);
	SD[it->first][k] = eq_value;
      }
    }
    //-----------------------------


    //----UPDATE EQUILIBRIUM ENERGY---------
    equilibrium_energy->zero();
    for (ID k = 0; k<AngInt.n_slices ; k++ ) 
    {
      sol_dir[k]->scale( AngInt.d_omega[k]);
      equilibrium_energy->add(*sol_dir[k]);
    }
    equilibrium_energy->scale(1.0/AngInt.total_angle);
    old_energy_norm = energy_norm;
    equilibrium_energy->close();
    energy_norm = equilibrium_energy->l2_norm();
    //-------------------------------------

    //----UPDATE THERMAL FLUX--------------
    for (ID i = 0; i<dim ; i++ )
    {
      thermal_flux[i]->close();
      thermal_flux[i]->zero();
      for (ID k = 0; k<AngInt.n_slices ; k++ ) 
      {
	sol_dir[k]->scale( IntDir(i)/AngInt.d_omega[k]);
	thermal_flux[i]->add(*sol_dir[k]);
      }
      thermal_flux[i]->scale(myopts.s_0 * myopts.s_0 * myopts.t_0/AngInt.total_angle);
    }


    {
      MeshBase::const_element_iterator       it    = mesh->active_elements_begin();
      const MeshBase::const_element_iterator it_end = mesh->active_elements_end(); 
    
      for ( ; it != it_end; ++it)
      {
	
	const Elem* elem = *it;
	ID subdomain = elem->subdomain_id();
	const Material* mat = _device->get_material(subdomain);
	HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
	double vg =  heat_model->get_phonon_group_velocity();
	vg *=myopts.s_0;
	vg /=myopts.t_0;
	ID dof = elem->dof_number(gray_sys_number,0,0);		        		        
	
	for (ID i = 0; i<dim ; i++ )    
	{
	  double value = (*thermal_flux[i])(dof);
	  thermal_flux[i]->set(dof, value * vg); 
	}
	
      }
    }

    //---------------------------------------------
    err = abs(energy_norm - old_energy_norm)/max(energy_norm,old_energy_norm);
 
    std::cout<<"Energy error: "<<err<<std::endl;

    iter +=1;
    
    } while (err > myopts.max_error & iter < myopts.max_iter);
  }
    
  energy_conservation_check();

}




//-------------------------------------------------------------------------//
MicroHeatBalance::~MicroHeatBalance()
{
  //Release pointers
  for (ID i = 0; i<dim; i++ )
    delete thermal_flux[i];
  
  for (ID k = 0; k<AngInt.n_slices ; k++ )
    delete sol_dir[k];


  delete  equilibrium_energy;



}
//---------------------------------------------------------------------------------//
MicroHeatBalance::MicroHeatBalance(const ModelOptions& options)
  : SimulationInterface(options)
{
 
 
}


PhysicalModel*
MicroHeatBalance::create_bulk_model(const ModelOptions& options,
			     const Material* mat) const
{
  return HeatModel::create(options);
}


//----------------------------------------------------------------------------------//
PhysicalModel*
MicroHeatBalance::create_physical_model(const ModelOptions &options,
    const Material* mat) const throw (ModelErrorException)
{

  HeatModel* model = dynamic_cast<HeatModel*> ( PhysicalModelInterface::create("thermal",options) );

  if (model == NULL)
    throw ModelErrorException("MicroHeatBalance: Thermal physical model is not created" );

  return model;

}
//----------------------------------------------------------------------------------//

BoundaryProperties* MicroHeatBalance::create_boundary_model (const ModelOptions &options) const
                    throw (ModelErrorException)

{

   const string& modelname = options.get_option("type", "heat_reservoir");


   ThermalContact* model = ThermalContact::create(modelname, options);

   if (model == NULL)
     throw ModelErrorException("MicroHeatBalance: No such boundary model: " + modelname);

  return model;

}
//----------------------------------------------------------------------------------//
MicroHeatBalance*  MicroHeatBalance::create (const ModelOptions& options)
{
  return new MicroHeatBalance(options);
}




//----------------------------------------------------------------------------------//
void MicroHeatBalance::assemble_heat_matrix(EquationSystems& es,
				     const std::string& system_name)
{

  static_this->do_assemble( es, system_name);

}


//----------------------------------------------------------------------------------//
void MicroHeatBalance::do_assemble(EquationSystems& es, const std::string& system_name)
{


  SimulationEnvironment& se = get_environment();

  TiberLinearSystem& system = *my_system;

  const unsigned int uvar = system.variable_number("T");

  DofMap& dof_map = system.get_dof_map();

  FEType fe_type = dof_map.variable_type(uvar);

  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true));

  //AutoPtr<FEBase> fe = FEBase::build(dim, fe_type);


  QGauss qrule (dim, CONSTANT); //may be could be decreased (CHECK!!!)

 
  fe -> attach_quadrature_rule (&qrule);

 // Here we define some references to cell-specific data that
  // will be used to assemble the lin ModelOptions&ear system.
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

  //------------------------------

  //Fe face

  // Declare a special finite element object for
  // boundary integration.
  //AutoPtr<FEBase> fe_face = FEBase::build(dim, fe_type);
  AutoPtr<FEBase>  fe_face(build_finite_element(dim,fe_type,true));

  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // o cout<<"Start loop over lattice thermal conductivity"<<endl;f the element.map
  QGauss qface(dim-1,CONSTANT);

  // Tell the finite element object to use our
  // quadrature rule.

  fe_face->attach_quadrature_rule(&qface);


  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();

  const std::vector<Real>& JxW_face = fe_face->get_JxW();

  const std::vector<Point>& qface_point = fe_face->get_xyz();

  const std::vector<Point>& normal = fe_face->get_normals();

  unsigned int system_number = my_system->number();
  const unsigned int var = my_system->variable_number("T");

  std::vector<unsigned int> dof_indices;

  DenseMatrix<Number>  Ke;

  DenseVector<Number>  Fe;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end(); 

  //Model Variables

  //Tensor2Sym kappa;

  ThermalContact* contact;
  //----------------------------------------------------------LatticeThermalConductivity-------//
 
  //std::cout<<"t_0: "<<t_0<<std::endl;
  //std::cout<<"s_0: "<<s_0<<std::endl;

  std::vector< std::map< ID, double > > temp_sol;

  for ( ; el != end_el ; ++el)   //loop over elements
  {

    const Elem* elem = *el;

    dof_map.dof_indices (elem,dof_indices);

    ID neighbor = elem->n_neighbors();

    for (ID k = 0; k < neighbor; k ++)
    {
      const Elem* elem_n = elem->neighbor(k);
      if (elem_n != NULL)
      {
	std::vector<unsigned int> dof_indices_n;
	dof_map.dof_indices (elem_n,dof_indices_n);
	dof_indices.push_back(dof_indices_n[0]);
	
      }
      else
	dof_indices.push_back(0);
    }



    const unsigned int n_dofs = dof_indices.size();


    e_0 = 0.0;

    fe->reinit(elem);

    Ke.resize(n_dofs,n_dofs);

    Fe.resize(n_dofs);

    Fe.zero();

    Ke.zero();

    const unsigned int num_sides = elem->n_sides();

    ID subdomain = elem->subdomain_id();
    const Material* mat = _device->get_material(subdomain);
    HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
    heat_model->re_init();

     
    std::vector<Point> p(1);
    p[0] = elem->centroid();

    std::vector<double> heat_source;
    heat_model->get_total_heat_source(elem,p,heat_source);

    vg =  heat_model->get_phonon_group_velocity();
    tg =  heat_model->get_phonon_scattering();
    cg =  heat_model->get_lattice_thermal_capacity(); //j/cm3K
   
    vg *=myopts.s_0;
    vg /=myopts.t_0;
    tg *=myopts.t_0;

    e_0 = (*equilibrium_energy)(dof_indices[0]);
  
    int fl = 1;
    
    if (fl==0)
    {
      //Assembly
      Ke(0,0) = 1.0/tg * JxW[0];
      Fe(0)   = (e_0/tg  +  heat_source[0] * d_omega /myopts.t_0) * JxW[0];

      for (ID ns = 0; ns < elem->n_sides(); ns++)
      {
	fe_face->reinit(elem,ns);
	const ElementSide elside(elem->top_parent(),ns);
	double in = dir * normal[0];
	double value = vg * in * JxW_face[0];

	//std::cout<<in<<std::endl;

	if (in<0.0)
	{
	  
	  if (se.is_on_boundary(elside))
	    Fe(0) -= get_boundary_value(elside) * value; 
	  else
	    Ke(0,ns + 1) += value;	  
 	}
	else
	  Ke(0,0) += value;

      }
      
    }
    else
    {
      //Assembly
      Ke(0,0) = 1.0/tg * JxW[0];
      Fe(0)   = (e_0/tg  +  heat_source[0]) * JxW[0];
      
      int OB = 0;
      for (unsigned int ns = 0; ns < elem->n_sides(); ns++)
      {
	fe_face->reinit(elem,ns);
	const ElementSide elside(elem->top_parent(),ns);
	double in = dir * normal[0];
	double value = vg * in * JxW_face[0];
	
	if (in<0.0)
	{
	  
	  if (se.is_on_boundary(elside))
	  {
	    Fe(0) = get_boundary_value(elside); 
	    Ke(0,0) = 1.0;
	    OB = 1; 
	  }
	  else
	    Ke(0,ns + 1) += value;
	  
	}
	else
	  if(OB == 0)
	    Ke(0,0) += value;


      }
      
      
    }
    
    
  
    // END ASSEMBLY
   
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices);


  } //End Loop over elements
  //system.matrix->print_matlab("Matr.m");
  //system.rhs->print();
//  std::cout<<"DONE"<<std::endl;
} //do assembly




double 
MicroHeatBalance::get_boundary_value(ElementSide elside)
{

  double A = 0.0;
  double D = 0.0;
  double R = 0.0;
  double T = 0.0; 
  double E = 0.0;
  double I = 0.0;

  SimulationEnvironment& se = get_environment();
  std::vector< std::map< ID, double > > temp_sol;
  Boundary* bd = se.get_boundary(elside);
  if (bd != NULL)
  {
    //It it's here it means that in the input file is defined  a contact.

    ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );
    
    switch (contact->get_type())
    {
    case  ThermalContact::Meso:
      {
	
	MesoThermalContact* meso_contact = (dynamic_cast<MesoThermalContact*> (contact) );
		
	A = meso_contact->get_absorbivity();
	D = meso_contact->get_diffusivity();
	R = meso_contact->get_reflectivity();
	T = meso_contact->get_temperature(); 

      }
    case  ThermalContact::Reservoir:
      {
	
	Reservoir* res_contact = (dynamic_cast<Reservoir*> (contact) );
	T = res_contact->get_temperature();
        A = 1.0;
      }
    }
   
  }
  else
  {
    if (myopts.diffusive)       
    {
      D = 1.0;
      E = e_0;
    }
    else
    {
      R = 1.0;
      I = SD[elside][vec_spec];
    }
    
  }

  
  double boundary_value = A * T + D * E + R * I;

  return boundary_value; 
}

ID
MicroHeatBalance::convert_variable_name_to_id(const string& variable_name) const
{

  ID id = INVALID_ID;

 
    if (variable_name == "temperature" )
       id  = TEMPERATURE;
    if (variable_name == "e_0" )
       id  = E0;
    if (variable_name == "Jqx" )
       id  = JQX;
    if (variable_name == "Jqy" )
       id  = JQY;
    if (variable_name == "Jqz" )
       id  = JQZ;




  return id;
}



void
MicroHeatBalance::get_solution_secure(const Elem* elem,
    const set<ID>& ids, vector<map<ID, double> >& values)
{


  std::vector<Point> points(elem->n_nodes());

  for (unsigned n = 0 ; n< elem->n_nodes(); ++n)
  {
    points[n] = elem->point(n);
  }

  get_solution_secure(elem,points,ids,values);



}



void
MicroHeatBalance::get_solution_secure(const Elem* elem, const vector<Point>& p,
    const set<ID>& ids, vector<map<ID, double> >& values)
{

 
  unsigned int np = p.size();
  values.resize(np);
  if ((np == 0) || (ids.size() == 0)) return;


  //vg for a given element
  ID subdomain = elem->subdomain_id();
  const Material* mat = _device->get_material(subdomain);
  HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
  double vg =  heat_model->get_phonon_group_velocity();
  //-------

  TiberLinearSystem& system = *my_system;

  DofMap& dof_map = system.get_dof_map();

  //const NumericVector<double>& solution = *(system.solution);
  AutoPtr<NumericVector<Number> >& solution = system.solution;
  const unsigned int var = system.variable_number("T");

  FEType fe_type = dof_map.variable_type(var);


  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));


  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  const std::vector<std::vector<RealGradient> >&  dphi = fe->get_dphi();

  vector<Point> points(np);

  FEInterface::inverse_map(dim, fe_type, elem, p, points);

  fe->reinit(elem, &points);

  vector<unsigned int> dof_indices;

  dof_map.dof_indices(elem, dof_indices);

  const unsigned int n_dofs = dof_indices.size();


 
  for (unsigned int n = 0; n < np; n++)
  {
    RealGradient Jq(0);
    double eq_energy = 0.0;
    for (unsigned int alpha = 0; alpha<dof_indices.size() ;alpha ++)
    {
       eq_energy += (*equilibrium_energy)(dof_indices[alpha]) * phi[alpha][n];
       //     std::cout<<eq_energy<<std::endl;


       for (ID i = 0; i < dim; i++)
	 Jq(i) += vg * (*thermal_flux[i])(dof_indices[alpha]) * phi[alpha][n] /AngInt.total_angle/myopts.t_0 /myopts.s_0/myopts.s_0;
       //       Jq(i) += vg * (*thermal_flux[i])(dof_indices[alpha]) * phi[alpha][n]/t_0 * cg/AngInt.total_angle;

       
    }

    if (ids.count(TEMPERATURE))
      values[n][TEMPERATURE] = eq_energy;
    //	  std::cout<<  Jq(0)<<std::endl;

    if (ids.count(E0))
      values[n][E0] = eq_energy;
    if (ids.count(JQX))
      values[n][JQX] = Jq(0);
    if (ids.count(JQY))
       values[n][JQY] = Jq(1);
    if (ids.count(JQZ))
      values[n][JQZ] = Jq(2);
      
    

  }

}



  //! Order the solution in correct mode
void MicroHeatBalance::build_elemental_results(const std::set<std::string>& variables,
				     std::vector<double>& results,
				     std::vector<std::string>& legend)
{
  
  unsigned int n_vars = 0;
  const unsigned int nn  = mesh->n_active_elem();
  const unsigned int dim = mesh->mesh_dimension();

  legend.resize(variables.size());

  int Temp = -1;
  if (variables.count("LatticeTempElem") ||
      variables.count("thermal"))  
    {
      Temp = n_vars;
      legend.resize(legend.size() + 1);
      legend[n_vars]="LatticeTemp";
      n_vars++;
    }
  
  int TF = -1;
  if (variables.count("ThermalFlux"))
  {
    TF = n_vars;
    legend.resize(legend.size() + dim + 2);
    switch (dim)
    {
    case 3:
	legend[TF + 2] = "thermal_flux_z";
	n_vars++;
      case 2:
	legend[TF + 1] =  "thermal_flux_y";
	n_vars++;
	legend[TF + dim] = "mod_thermal_flux" ;
	n_vars++;
      default:
	legend[TF] = "thermal_flux_x";
	n_vars++;
      }
     
    }



  int par = -1;
  if (variables.count("partial") ||
      variables.count("thermal")   )
  {
    par = n_vars;
    for (ID k = 0;k<  AngInt.n_slices;k++)
    {
      
      std::string label;
      std::ostringstream i_str;
      i_str << "energy_density" << k;
      legend.push_back(i_str.str());
      n_vars++;
    }
  }
  results.resize(nn * n_vars,0.0);
  legend.resize(n_vars);
  const unsigned int  var = my_system->variable_number("T");

  DofMap& dof_map =  my_system->get_dof_map();
  std::vector<unsigned int> dof_indices;

  FEType fe_type = dof_map.variable_type(var);
  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
  
  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;

    dof_map.dof_indices (elem, dof_indices);

    unsigned int id = n_vars * elem_number;

    if (Temp != -1)  
      results[id + Temp] =  (*equilibrium_energy)(dof_indices[0]);
    
    if (par != -1)
      for (ID k = 0;k<  AngInt.n_slices;k++)
	results[id+par+k] = (*sol_dir[k])(dof_indices[0]);

    if (TF != -1)
    {
      unsigned int k = 0;
      ID dof = elem->dof_number(gray_sys_number,0,0);
      
      RealGradient P(0);
      for (ID d = 0; d<dim; d++)
	P(d) = (*thermal_flux[d])(dof);
	switch (dim)
	{
	case 3:
	  results[id + TF  + 2] = P(2);
	case 2:
	  results[id + TF + 1] = P(1);
	  results[id + TF + dim] = 0.0;
	default:
	  results[id + TF ] = P(0);
	}

    }


    elem_number++;
  }

  
 results.resize(elem_number * n_vars);

}

//----------------------------------------------------------------------------------//
void MicroHeatBalance::build_nodal_results (const std::set< std::string > &variables,
				     std::vector< double > &results,
				     std::vector< std::string > &legend)
{
  //----------------Fourier test
  const NumericVector<double>& fourier_solution = *(fourier_system->solution);
  DofMap& dof_map_fourier = fourier_system->get_dof_map();
  std::vector<unsigned int> dof_indices_fourier;
  //-------------------------------

  legend.resize(0);
  unsigned int n_vars = 0;

  int Temp = -1;
  if (variables.count("LatticeTemp")||
      variables.count("thermal")    ) 
  {
    Temp = n_vars;
    legend.push_back("LatticeTemp");
    n_vars++;
  }

  int FourierTemp = -1;
  if (variables.count("FourierTemp")||
      variables.count("thermal")    ) 
  {
    FourierTemp = n_vars;
    legend.push_back("FourierTemp");
    n_vars++;
  }

  int par = -1;
  if (variables.count("partial") ||
      variables.count("thermal")   )
  {
    par = n_vars;
    for (ID k = 0;k<  AngInt.n_slices;k++)
    {
      
      std::string label;
      std::ostringstream i_str;
      i_str << "energy_density" << k;
      legend.push_back(i_str.str());
      n_vars++;
    }
  }

  //Create node connection    
  const unsigned int nn  = mesh->n_nodes();
  vector<unsigned short int> node_conn(nn);
  {
    vector<unsigned short int> node_conn_local(node_conn.size());
    
    MeshBase::const_element_iterator it =
      get_mesh().active_local_elements_begin();
    const MeshBase::const_element_iterator end =
      get_mesh().active_local_elements_end();
    
    for ( ; it != end; ++it)
      for (unsigned int n = 0; n < (*it)->n_nodes(); n++)
	node_conn_local[(*it)->node(n)]++;
    
    node_conn = node_conn_local;
  }
  
  
  results.resize(nn * n_vars,0.0);
  legend.resize(n_vars);
  
  
  std::vector<unsigned int> dof_indices;
  DofMap& dof_map = my_system->get_dof_map();
  
  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
  
  std::vector< std::map< ID, double > > temp_sol;
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    
    dof_map.dof_indices (elem, dof_indices);

    dof_map_fourier.dof_indices(elem, dof_indices_fourier);

    for (int n = 0; n < elem->n_nodes(); n++)
    {
      unsigned int id =  (elem->node(n) * n_vars) ;
      
      if (Temp != -1)
	results[id+Temp]  += (*equilibrium_energy)(dof_indices[0])/node_conn[elem->node(n)];

      if (FourierTemp != -1)
	results[id+FourierTemp] = fourier_solution(dof_indices_fourier[n]);

      if (par != -1)
	for (ID k = 0;k<  AngInt.n_slices;k++)
	  results[id+par+k] += (*sol_dir[k])(dof_indices[0])/node_conn[elem->node(n)];

    }
    
  }

}


void
MicroHeatBalance::build_integrated_quantities_description(
    std::vector<std::string>& legend,
    std::vector<std::string>& description)
{


  if (plot_solution("PowerDissipated"))
  {
    legend.resize(1);

    description.resize(1);

    const unsigned int dim = mesh->mesh_dimension();

    ostringstream s;
    s << "Power Dissipated. Units W";
    switch (dim)
    {
      case 1:
        s << "cm^-2";
        break;
      case 2:
        s << "cm^-1";
        break;
    }
    description[0] = s.str();
  }
}



void
MicroHeatBalance::build_integrated_quantities(
    vector<double>& values)
{


  if (plot_solution("PowerDissipated"))
  {

    // double power = calculate_power_dissipated();

    //  values.resize(1,power);


  }
}



double
MicroHeatBalance::energy_conservation_check()
{



  SimulationEnvironment& se = get_environment();

 //  //------------------------------Gray System---------------
  DofMap& dof_map = my_system->get_dof_map();
  std::vector<unsigned int> dof_indices;
  

//   const unsigned int var = my_system->variable_number("T");
//   FEType fe_type_gray = dof_map.variable_type(var);
//   AutoPtr<FEBase>  fe_gray(build_finite_element(dim,fe_type_gray,true));
//   QGauss qrule_gray(dim, CONSTANT);
//   fe_gray->attach_quadrature_rule(&qrule_gray);
//   const std::vector<std::vector<Real> >&  phi = fe_gray->get_phi();
//   const std::vector<Real>& JxW_gray = fe_gray->get_JxW();
//   const std::vector<Point>& q_point_gray = fe_gray->get_xyz();
//   //------------------------------------------------------------

  FEType fe_type(FIRST,LAGRANGE);
  AutoPtr<FEBase>  fe(build_finite_element(dim,fe_type,true));

  QGauss qrule(dim, CONSTANT);

  fe->attach_quadrature_rule(&qrule);
  const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<RealGradient> >&  dphi_rstf = fe->get_dphi();
  const std::vector<Real>& JxW = fe->get_JxW();

 

  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();

  double check = 0.0;
  //  double PowerEmitted = 0.0;
  for ( ; it != end; ++it)
  {

    const Elem* elem = *it;  
    dof_map.dof_indices(elem, dof_indices);


   //  //-------------------------PowerEmitted------------------------------------------------
//     fe_gray->reinit(elem);
//     ID subdomain = elem->subdomain_id();
//     dof_map.dof_indices(elem, dof_indices,var);
//     const Material* mat = _device->get_material(subdomain);
//     HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
//     heat_model->re_init();
//     std::vector<double> heat_source;
//     heat_model->get_total_heat_source(elem,q_point_gray,heat_source);
//     PowerEmitted  += JxW_gray[0] * heat_source[0]/myopts.s_0 /myopts.s_0 /myopts.s_0;
//     //-------------------------------------------------------------------------------------
 

    //-------------------------PowerDissipated------------------------------------------------
    bool has_node = false;
    const ID num_sides = elem->n_sides();
    for (ID ns = 0; ns<num_sides; ns++)
    {
      const ElementSide elside(elem->top_parent(),ns);
      if (se.is_on_boundary(elside)) //natural boundary
	has_node = true; 
    }

    if (!has_node)
      continue;
    
    fe->reinit(elem);
    
    for (ID d = 0; d<dim; d++)
      for (ID n = 0; n < elem->n_nodes() ;n ++)
      {
	check += JxW[0] * (*thermal_flux[d])(dof_indices[0]) * dphi_rstf[n][0](d);
        cout<<JxW[0] * (*thermal_flux[d])(dof_indices[0]) * dphi_rstf[n][0](d)<<endl;
      }
    //----------------------------------------------------------------------------------------
	    	
  }

  if (SimulationOptions::verbose() > 1)  
    std::cout<<"Energy conservation: "<<check<<std::endl;
   
   



  return check; 
 
}





void
MicroHeatBalance::do_print_info(void)
{

  string space("  ");
  //cout << space << "linear solver is: petsc" <<std::endl;

}



void
MicroHeatBalance::AngularIntegrator::compute_custom_direction(std::vector<Point> custom_dir)
{

  //Omega is not theta dependent but simply uniform. 
  weight = 0.0;
  total_angle = 4.0 * M_PI;
  
  n_slices = custom_dir.size();  

  directions.resize(n_slices);
  d_omega.resize(n_slices);
  dir.resize(n_slices);
  theta_vec.resize(n_slices);
  phi_vec.resize(n_slices);
  spec.resize(1);

  dir = custom_dir;
  // directions = custom_dir * 2*M_PI;

  for (ID k = 0; k<n_slices; k++)
  {

    directions[k] = custom_dir[k];// * M_PI;
    d_omega[k] = 4.0 * M_PI/n_slices;
    directions[k] = custom_dir[k] * d_omega[k];
    theta_vec[k] = 0.0;
    phi_vec[k] = 0.0;
   
  }
  
  //Spec vectors
  for (ID k1 = 0;k1<n_slices; k1++)
    for (ID k2 = 0;k2<n_slices; k2++)
    {
      Point sum = directions[k1] + directions[k2];
      if (sum.size() < 1e-4)
	spec[k1]=k2;
      
    }
 

}
void
MicroHeatBalance::AngularIntegrator::compute_directions()
{

    double min_theta, max_theta, min_phi, max_phi;
   
    switch (dim)
    {
      
    case 1 :
      
      theta_slices = 1;
      phi_slices = 2;
      
      min_theta = 0.0;
      max_theta = M_PI;

      min_phi = M_PI * 0.5;
      max_phi = M_PI * 2.0 + M_PI * 0.5;
      
      weight = 2.0 * M_PI;
      
      n_slices = theta_slices * phi_slices;
      spec.resize(n_slices);
      //  spec[0] = 1;
      //spec[1] = 0;
      
      break;
      
    case 2 :
      
      std::cout<<"2D"<<std::endl;
      
      theta_slices = 1;
      if (phi_slices == 0)
	phi_slices = 4;
      
      
      min_theta = 0.0;
      max_theta = M_PI;
      
      min_phi = M_PI * 0.5;
      max_phi = M_PI * 2.0 + M_PI * 0.5;
      
      weight =  2.0;
      
      n_slices = theta_slices * phi_slices;
      spec.resize(n_slices);
      //  spec[0] = 2;
      //spec[1] = 3;
      //spec[2] = 0;
      //spec[3] = 1;
      
      break;
      
    case 3 :
      
      std::cout<<"3D"<<std::endl;
      
      if (theta_slices == 0)
	theta_slices = 4;
      
      if (phi_slices == 0)
	phi_slices = 2;
      
      
      min_theta = 0.0;
      max_theta = M_PI;
      
      min_phi = 0.0;
      max_phi = M_PI * 2.0;
      
      weight = 1.0;
      
      
      n_slices = theta_slices * phi_slices;

      spec.resize(n_slices);

     
      
      break;
    }
    


    total_angle = 4.0 * M_PI;
    
    directions.resize(n_slices);
    integrate_directions.resize(n_slices);
    d_omega.resize(n_slices);
    dir.resize(n_slices);
    theta_vec.resize(n_slices);
    phi_vec.resize(n_slices);
    
    d_theta =  (max_theta - min_theta) / theta_slices;
    d_phi =  (max_phi - min_phi) / phi_slices;
    double theta, phi;


    ID k = 0;
    for (ID n_phi = 0; n_phi < phi_slices; n_phi++)
    {
      // phi = min_phi + d_phi * 0.5 + d_phi * n_phi;
      phi = min_phi + d_phi * n_phi;
      
      for (ID n_theta = 0; n_theta < theta_slices; n_theta++)
      {
	theta = min_theta + d_theta * 0.5 + d_theta * n_theta;
	
	d_omega[k] =  2.0 * sin(theta) * sin (0.5 * d_theta) * d_phi;
	
	dir[k](0) = sin(theta) * sin(phi);
	dir[k](1) = sin(theta) * cos(phi);
	dir[k](2) = cos(theta);
	
	//directions[k](0) = weight * sin(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
	//directions[k](1) = weight * cos(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
	//directions[k](2) = weight * 0.5 * d_phi * sin(2.0 * theta) * sin(d_theta);
	
	directions[k](0) =  sin(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
	directions[k](1) =  cos(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
	directions[k](2) =  0.5 * d_phi * sin(2.0 * theta) * sin(d_theta);
	
	theta_vec[k] = theta;
	phi_vec[k] = phi;
	
	k++;
      }
    }

    //Spec vectors
    for (ID k1 = 0;k1<n_slices; k1++)
      for (ID k2 = 0;k2<n_slices; k2++)
      {
	Point sum = directions[k1] + directions[k2];
        if (sum.size() < 1e-4)
	  spec[k1]=k2;
	
      }
      
 
}


void
MicroHeatBalance::AngularIntegrator::print_info(void)
{


  std::cout<<"General:  "<<"theta slices: "<<theta_slices<<" phi slices: "<<phi_slices<<std::endl;
  std::cout<<"d_theta:  "<<d_theta / ( M_PI) * 180.0<<"  d_phi:  "<<d_phi / ( M_PI) * 180.0<<std::endl;
  std::cout<<"Weight factor:  "<<weight<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"  "<<std::endl;


  for (ID k =0; k<n_slices; k ++)
  {

    std::cout<<"Direction: "<<k+1<<std::endl;
    std::cout<<"  Theta: "<<theta_vec[k] / ( M_PI) * 180.0 <<"  phi: "<<phi_vec[k] / ( M_PI) * 180.0<<std::endl;

    std::cout<<"d_omega:  "<<  d_omega[k]/(4.0 * M_PI) <<std::endl;


    std::cout<<"sx:  "<<   dir[k](0)  <<"  sy:  "<<   dir[k](1) <<"  sz:  "<<   dir[k](2)<<std::endl;
    std::cout<<"six:  "<<   directions[k](0)  <<"  siy:  "<<   directions[k](1) <<"  siz:  "<<   directions[k](2)<<std::endl;
    std::cout<<"  "<<std::endl;

  }

  std::cout<<"  "<<std::endl;

  //check 1
  double solid_angle = 0.0;
  double s0 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  for (ID n = 0; n < n_slices; n++)
  {
    solid_angle += d_omega[n] ;
    s0 += directions[n](0);
    s1 += directions[n](1);
    s2 += directions[n](2);

  }
  std::cout<<"Check:  "<<std::endl;
  std::cout<<"Solid angle: "<<solid_angle/(4 * M_PI)<<std::endl;


  std::cout<<"Directions:  "<<std::endl;
  std::cout<<s0<<std::endl;
  std::cout<<s1<<std::endl;
  std::cout<<s2<<std::endl;



  std::cout<<"  "<<std::endl;
  std::cout<<"TotalAngle:  "<<total_angle<<std::endl;
  std::cout<<"  "<<std::endl;
  std::cout<<"Specular vector: "<<std::endl;
  for (ID k =0; k<n_slices;k++)
    std::cout<<"Dir "<<k<<" : "<<spec[k]<<std::endl;
  std::cout<<"  "<<std::endl;
}


void
MicroHeatBalance::AngularIntegrator::print_info(ID k)
{



    std::cout<<"Direction: "<<k<<std::endl;
    std::cout<<"  Theta: "<<theta_vec[k] / ( M_PI) * 180.0 <<"  phi: "<<phi_vec[k] / ( M_PI) * 180.0<<std::endl;

    std::cout<<"d_omega:  "<<  d_omega[k]/(4.0 * M_PI) <<std::endl;


    std::cout<<"sx:  "<<   dir[k](0)  <<"  sy:  "<<   dir[k](1) <<"  sz:  "<<   dir[k](2)<<std::endl;
    std::cout<<"six:  "<<   directions[k](0)  <<"  siy:  "<<   directions[k](1) <<"  siz:  "<<   directions[k](2)<<std::endl;
    std::cout<<"  "<<std::endl;
  
  

  std::cout<<"  "<<std::endl;

}


NumericVector< double > & 
MicroHeatBalance::do_get_solution_vector(void)
{
  return  *equilibrium_energy;
}

void
MicroHeatBalance::compute_fourier_solution(void)
{


  std::cout<<"Compute Fourier Solution"<<std::endl;

  fourier_system = TiberLinearSystem::create(get_equation_systems(),
					     "fourier", get_solver_options()); 

  fourier_system->add_variable("T", FIRST);
  
  fourier_system->attach_assemble_function(assemble_macro_heat_matrix);

  fourier_system->init();

  fourier_system->solve();

  //Converting nodal data in elemental data

  //Fourier system INIT
  DofMap& dof_map_fourier = fourier_system->get_dof_map();
  std::vector<unsigned int> dof_indices_fourier;
  
  const NumericVector<double>& solution = *(fourier_system->solution);
  const unsigned int var = fourier_system->variable_number("T");
  FEType fe_type = dof_map_fourier.variable_type(var);

  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  const vector<vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  //Gray System INIT
  DofMap& dof_map = my_system->get_dof_map();
  std::vector<unsigned int> dof_indices;
  
  for (ID j = 0;j<dim; j++)
    thermal_flux[j]->zero();
  
  equilibrium_energy->zero();

  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
  
  for ( ; it != end; ++it)
  {

    const Elem* elem = *it;  

    ID subdomain = elem->subdomain_id();
    const Material* mat = _device->get_material(subdomain);
    HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
    heat_model->re_init();
    vg =  heat_model->get_phonon_group_velocity();
    tg =  heat_model->get_phonon_scattering();
    cg =  heat_model->get_lattice_thermal_capacity();
    
    double cg = 1.0; //j/cm-3
    vg *=myopts.s_0;
    cg /=myopts.s_0;
    cg /=myopts.s_0;
    cg /=myopts.s_0;

    double kappa = cg * vg * vg * tg /3.0;

    dof_map_fourier.dof_indices (elem, dof_indices_fourier);

    dof_map.dof_indices (elem, dof_indices);

    const unsigned int n_dofs = dof_indices.size();

    //Init of the elem at the centroid
    std::vector<Point> p(1);
    p[0]=elem->centroid();
    vector<Point> points(1);
    FEInterface::inverse_map(dim, fe_type, elem, p, points);
    fe->reinit(elem, &points);
    
    //Compute temperature at the centroid
    double value = 0.0;
    RealGradient TF(0);
    for (ID alpha = 0; alpha<dof_indices_fourier.size() ;alpha ++) 
    {
      value += solution(dof_indices_fourier[alpha]) *  phi[alpha][0];
      
      for (ID j = 0;j<dim; j++)
	TF(j) +=  - kappa * solution(dof_indices_fourier[alpha]) * dphi[alpha][0](j) * myopts.s_0 * myopts.s_0;

    }
    
    for (ID j = 0;j<dim; j++)
      thermal_flux[j]->set(dof_indices[0],TF(j));
   
    equilibrium_energy->set(dof_indices[0],value);
        
  }

  old_energy_norm = energy_norm;
  equilibrium_energy->close();
  energy_norm =  equilibrium_energy->l2_norm();
   

}



//----------------------------------------------------------------------------------//

void 
MicroHeatBalance::assemble_macro_heat_matrix(EquationSystems& es,
				     const std::string& system_name)
{

   static_this->do_macro_assemble( es, system_name);

}

//----------------------------------------------------------------------------------//
void 
MicroHeatBalance::do_macro_assemble(EquationSystems& es, const std::string& system_name)
{
  SimulationEnvironment& se = get_environment();
  TiberLinearSystem& f_system = *fourier_system;
  const unsigned int uvar = f_system.variable_number("T");
  DofMap& dof_map = f_system.get_dof_map();
  FEType fe_type = dof_map.variable_type(uvar);
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true));
  QGauss qrule (dim, FIFTH); //may be could be decreased (CHECK!!!)
  fe -> attach_quadrature_rule (&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  
  std::vector<unsigned int> dof_indices;
  DenseMatrix<Number>  Ke;
  DenseVector<Number>  Fe;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  ThermalContact* contact;
  //----------------------------------------------------------LatticeThermalConductivity-------//
  Tensor2Sym kappa2;

  for ( ; el != end_el ; ++el)   //loop over elements
  {
    const Elem* elem = *el;
    dof_map.dof_indices (elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();

    fe->reinit(elem);
    Ke.resize(n_dofs,n_dofs);
    Fe.resize(n_dofs);
    Fe.zero();
    Ke.zero();

    ID subdomain = elem->subdomain_id();
    const Material* mat = _device->get_material(subdomain);
    HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
    heat_model->re_init();

    vg =  heat_model->get_phonon_group_velocity();
    tg =  heat_model->get_phonon_scattering();
    cg =  heat_model->get_lattice_thermal_capacity();

    double cg = 1.0; //j/cm-3
    vg *=myopts.s_0;
    //vg /=myopts.t_0;
    //tg *=myopts.t_0;
    cg /=myopts.s_0;
    cg /=myopts.s_0;
    cg /=myopts.s_0;

    double kappa = cg * vg * vg * tg /3.0;
    //s_0 = 1.0;
    //
    // kappa /=myopts.s_0;
    //kappa /=myopts.t_0;
  					
    std::vector<double> heat_source;
    heat_model->get_total_heat_source(elem,q_point,heat_source);
    
   
    //LHS
    for (unsigned int p1=0; p1<n_dofs; p1++) 
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	for (unsigned int p2 =0; p2<n_dofs; p2++)
	  for (short i = 0; i < dim; i++)
	    Ke(p1,p2) += JxW[qp] * dphi[p1][qp](i) * kappa * dphi[p2][qp](i);
 

    //RHS
    for (unsigned int p1=0; p1<n_dofs; p1++) 
      for (unsigned int qp=0; qp<qrule.n_points(); qp++) 
	Fe(p1) +=JxW[qp] * heat_source[qp] * phi[p1][qp] /myopts.s_0 /myopts.s_0/myopts.s_0;
    

     
    //Boundary Condition
    const unsigned int num_sides = elem->n_sides();
    for (unsigned int side = 0; side<num_sides; side++)
    {
      const ElementSide elside(elem->top_parent(), side);

 
      heat_model->re_init();

      Boundary* bd = se.get_boundary(elside);
      
      if (bd != NULL)
      {
	if (bd->get_boundary_properties( get_id() ) != NULL )
	{

	  ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );

	  switch (contact->get_type())
	  {
	  case  ThermalContact::Reservoir:
	    {
	      for(unsigned int n = 0; n< n_dofs; ++n)
	      {
		if (elem->is_node_on_side(n,side))
		{
		  for (unsigned int nc = 0; nc < n_dofs; nc++)
		    Ke(n,nc) = 0.0;
		  
		  Ke(n,n) = 1.0;
 
                  double value = (dynamic_cast<Reservoir*> (contact) )->get_temperature();
               
		  Fe(n) = value;
		 
		}
	      }

	    }
	    break;

	   case  ThermalContact::Diffusive:
	    {
	      for(unsigned int n = 0; n< n_dofs; ++n)
	      {
		if (elem->is_node_on_side(n,side))
		{
		  for (unsigned int nc = 0; nc < n_dofs; nc++)
		    Ke(n,nc) = 0.0;
		  
		  Ke(n,n) = 1.0;
		  Fe(n) = (dynamic_cast<Diffusive*> (contact) )->get_temperature();
		  
		}
	      }
	      
	    }
	    break;
	  case  ThermalContact::FourierBTE:
	    {
	      for(unsigned int n = 0; n< n_dofs; ++n)
	      {   
		if (elem->is_node_on_side(n,side))
		{
		  const Elem* neighbour = (elside.elem())->neighbor(side);
		  double temp = (dynamic_cast<FourierBTE*> (contact) )->get_temperature(neighbour,elem->point(n));
		  for (unsigned int nc = 0; nc < n_dofs; nc++)
		    Ke(n,nc) = 0.0;
		  
		  Ke(n,n) = 1;
		  Fe(n) = temp; 
		
		}
	      }
	    }
	    
	    break;

	  }//switch
	}
      }
    }

    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    f_system.matrix->add_matrix (Ke, dof_indices);
    f_system.rhs->add_vector    (Fe, dof_indices);


  }//end loop over elements


}





// //----------------------------------------------------------------------------------//
// void MicroHeatBalance::assemble_heat_matrix_DG(EquationSystems& es,
// 				     const std::string& system_name)
// {

//    static_this->do_assemble_DG( es, system_name);

// }
//----------------------------------------------------------------------------------//
// void MicroHeatBalance::do_assemble(EquationSystems& es, const std::string& system_name)
// {
//   std::cout<<" "<<std::endl;
//   std::cout<<"LAGRANGE mode"<<std::endl;
//   std::cout<<" "<<std::endl;


//   SimulationEnvironment& se = get_environment();

//   TiberLinearSystem& system = *my_system;

//   const unsigned int uvar = system.variable_number("T");

//   DofMap& dof_map = system.get_dof_map();

//   FEType fe_type = dof_map.variable_type(uvar);

//   AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true));

//   //AutoPtr<FEBase> fe = FEBase::build(dim, fe_type);


//   QGauss qrule (dim, FIFTH); //may be could be decreased (CHECK!!!)

//   // quadrature rule
//   fe -> attach_quadrature_rule (&qrule);


//  // Here we define some references to cell-specific data that
//   // will be used to assemble the lin ModelOptions&ear system.
//   //
//   // The element Jacobian * quadrature weight at each integration point.
//   const std::vector<Real>& JxW = fe->get_JxW();


//   // The physical XY locations of the quadrature points on the element.
//   // These might be useful for evaluating spatially varying material
//   // properties at the quadrature points.

//   const std::vector<Point>& q_point = fe->get_xyz();

//   // The element shape functions evaluated at the quadrature points.


//   const std::vector<std::vector<Real> >& phi = fe->get_phi();

//   // The element shape function gradients evaluated at the quadrature
//   // points.
//   const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

//   //------------------------------

//   //Fe face

//   // Declare a special finite element object for
//   // boundary integration.
//   //AutoPtr<FEBase> fe_face = FEBase::build(dim, fe_type);
//   AutoPtr<FEBase>  fe_face(build_finite_element(dim,fe_type,true));

//   // Boundary integration requires one quadraure rule,
//   // with dimensionality one less than the dimensionality
//   // o cout<<"Start loop over lattice thermal conductivity"<<endl;f the element.map
//   QGauss qface(dim-1, SIXTH);

//   // Tell the finite element object to use our
//   // quadrature rule.

//   fe_face->attach_quadrature_rule(&qface);


//   const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();

//   const std::vector<Real>& JxW_face = fe_face->get_JxW();

//   const std::vector<Point>& qface_point = fe_face->get_xyz();

//   const std::vector<Point>& normal = fe_face->get_normals();

//   unsigned int system_number = my_system->number();
//   const unsigned int var = my_system->variable_number("T");

//   std::vector<unsigned int> dof_indices;

//   DenseMatrix<Number>  Ke;

//   DenseVector<Number>  Fe;

//   MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
//   const MeshBase::const_element_iterator end_el = mesh->active_elements_end(); 

//   //Model Variables

//   //Tensor2Sym kappa;

//   ThermalContact* contact;
//   //----------------------------------------------------------LatticeThermalConductivity-------//
 
//   // std::cout<<"t_0: "<<t_0<<std::endl;
//   //std::cout<<"s_0: "<<s_0<<std::endl;

//   std::vector< std::map< ID, double > > temp_sol;

//   for ( ; el != end_el ; ++el)   //loop over elements
//   {

//     const Elem* elem = *el;

//     dof_map.dof_indices (elem, dof_indices);


//     const unsigned int n_dofs = dof_indices.size();

//     double e_0 = 0.0;

//     fe->reinit(elem);

//     Ke.resize(n_dofs,n_dofs);

//     Fe.resize(n_dofs);

//     Fe.zero();

//     Ke.zero();

//     const unsigned int num_sides = elem->n_sides();


//     ID subdomain = elem->subdomain_id();

//     const Material* mat = _device->get_material(subdomain);

//     HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
 
//     //HeatModel* heat_model = (  dynamic_cast<HeatModel*> ( get_physical_model(subdomain)  ) );

//     // heat_model->set_element(elem);

//     //heat_model->set_side(-1);

//     heat_model->re_init();

//     std::vector<double> heat_source;
//     heat_model->get_total_heat_source(elem,q_point,heat_source);

//     vg =  heat_model->get_phonon_group_velocity();
//     tg =  heat_model->get_phonon_scattering();
//     cg =  heat_model->get_lattice_thermal_capacity(); //j/cm3K

   
//     vg *=myopts.s_0;
//     vg /=myopts.t_0;
//     tg *=myopts.t_0;
   

//     // std::vector<Point> p(1);
//     //p[0] = elem->centroid();


//       get_solution_secure(elem,q_point,e0_var,temp_sol);

//       //   get_solution_secure(elem,p,e0_var,temp_sol);
//       //   e_0 = temp_sol[0].find(E0)->second;

//     for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//     {//Loop over quadrature points
      
//         e_0 = temp_sol[qp].find(E0)->second;

      
//       //  std::cout<<e_0<<std::endl;
//       for (unsigned int p1=0; p1<n_dofs; p1++) 
//       { // loop over test function
	
// 	for (unsigned int p2=0; p2<n_dofs; p2++)
// 	{//loop over basis functions
	  
// 	  double value1 = 0.0;
// 	  double value2 = 0.0;
	
// 	  value1 = vg * JxW[qp] * (dphi[p2][qp] * IntDir) * phi[p1][qp];

// 	  //value1 = vg * JxW[qp] * (dphi[p1][qp] * IntDir) * phi[p2][qp];
// 	  value2 = 1.0/tg * JxW[qp] *  phi[p1][qp] * phi[p2][qp] * d_omega;

//  	  Ke(p1,p2) += value1;
// 	  Ke(p1,p2) += value2;

// 	  //  std::cout<<"v1: "<<value1<<std::endl;
// 	  // std::cout<<"v2: "<<value2<<std::endl;


// 	} //loop over basis functions

// 	double value1 = 0.0;
// 	double value2 = 0.0;

// 	value1 =  1.0/tg * JxW[qp] * e_0 * phi[p1][qp]  * d_omega;
// 	//std::cout<<"v3: "<<value1<<std::endl;
 
// 	value2 =  JxW[qp]  *  heat_source[qp] * phi[p1][qp] * d_omega /myopts.t_0;

// 	Fe(p1) += value1;
//  	Fe(p1) += value2;

//       }//end Loop over quadrature points
//     } // end loop over test function

   
//     //Second part of assembly (surface)
//     for (unsigned int side = 0; side < elem->n_sides(); side++)
//     {
//       const ElementSide elside(elem->top_parent(), side);
//       if (se.is_outer_boundary(elside))
//       {
// 	fe_face->reinit(elem,side);
// 	for (unsigned int qp=0; qp < qface.n_points(); qp++)
// 	{
// 	  for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
// 	  {
// 	    for (unsigned int p2=0; p2<n_dofs; p2++) //bases test
// 	    {
// 	      double val_plus = 0.0;
// 	       val_plus  =  JxW_face[qp] * vg * phi_face[p1][qp] * phi_face[p2][qp] * (IntDir * normal[qp]);
// 	       //Ke(p1,p2) -= val_plus;

// 	    }// (unsigned int p2=0; p2<n_dofs; p2++)
// 	  }//for (unsigned int p1=0; p1<n_dofs; p1++)
// 	}// for (unsigned int qp=0; qp < qface.n_points(); qp++)
 
//       }
//     }
    

    
//     //Boundary conditions and source
//     //The loop over element is the only loop that is surviving at this point

//     for (unsigned int side = 0; side<num_sides; side++)
//     {

//       const ElementSide elside(elem->top_parent(), side);

//       if (se.is_on_boundary(elside))
//       {
// 	heat_model->re_init();    
// 	fe_face->reinit(elem,side);

// 	{
// 	  //If no boundary is defined
// 	  //  Boundary* bd = se.get_boundary(elside);
// 	  //if (bd != NULL)
// 	  //{std::cout<<"Ciaooo"<<std::endl;
// 	  // if (bd->get_boundary_properties( get_id() ) == NULL )
// 	  // {
// 	  //  std::cout<<"Ciaooo2"<<std::endl;
// 	      //Default Boundary Properties-----Specular

// 	      double in = dir * normal[0];
	   
   
// 	      if (in <0.0 + 1e-6)
// 	      {
// 		for(unsigned int n = 0; n< n_dofs; ++n)
// 		{
// 		  if (elem->is_node_on_side(n,side))
// 		  {
// 		    for (unsigned int nc = 0; nc < n_dofs; nc++)
// 		     Ke(n,nc) = 0.0;
		    
// 		    Ke(n,n) = 1;
		    
// 		    //Ke(n,n) = 1e36;
		    
// 		    const unsigned int  n_dof = (elem->get_node(n))->dof_number(system_number,var,0); 
		    
// 		    double value = 0.0;
		    
// 		    if (myopts.diffusive)
// 		      value = e_0;
// 		    else
// 		      value =  bv[elem->get_node(n)][vec_spec];
		      
                     
		    
// 		    //Fe(n) = value *1e36; //From temp to energy
	
// 		    Fe(n) = value;
		    
// 		  }
// 		}
		
// 	      }
// 	      //	}
// 	      //}
// 	}
// 	    //---------------------------------


// 	Boundary* bd = se.get_boundary(elside);
      
// 	if (bd != NULL)
// 	{
// 	  if (bd->get_boundary_properties( get_id() ) != NULL )

// 	  {
	    
	    
// 	    ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );
	    
// 	    switch (contact->get_type())
// 	    {
// 	    case  ThermalContact::Reservoir:
// 	      {
		
// 		double in = dir * normal[0];
// 		double temp = (dynamic_cast<Reservoir*> (contact) )->get_temperature();
		
		
// 		if (in < 0.0 + 1e-6)
// 		{
// 		 for(unsigned int n = 0; n< n_dofs; ++n)
// 		 {
// 		   if (elem->is_node_on_side(n,side))
// 		   {
		  
// 		     //Ke(n,n) = 1e36;
// 		     // Fe(n) = 1e36* temp;
		       
		   		  
// 		       for (unsigned int nc = 0; nc < n_dofs; nc++)
// 		       Ke(n,nc) = 0.0;
		     
// 		       Ke(n,n) = 1;
// 		       Fe(n) = temp;// * cg /AngInt.total_angle * myopts.scale;
		    

  
// 		   }
//                  }
// 	       }
// 	      }
// 	      break;
	      
// 	    case  ThermalContact::FourierBTE:
// 	      {
		
// 		double in = dir * normal[0];
                     
// 		if (in < 0.0)
// 		{
		  
// 		  for(unsigned int n = 0; n< n_dofs; ++n)
// 		  {
// 		    if (elem->is_node_on_side(n,side))
// 		    {
		      
// 		      const Elem* neighbour = (elside.elem())->neighbor(side);
// 		      double temp = (dynamic_cast<FourierBTE*> (contact) )->get_temperature(neighbour,elem->point(n));
		      
// 		      for (unsigned int nc = 0; nc < n_dofs; nc++)
// 			Ke(n,nc) = 0.0;
		      
// 		      Ke(n,n) = 1;
		     
// 		      Fe(n) =  temp;
		      
// 		    }
// 		  }
// 		}
		
// 		break;
// 	      }
	      
// 	     case  ThermalContact::Specular:	      
// 	      {
		
	
// 		double in = dir * normal[0];
// 				//	std::cout<<in<<std::endl;
// 		if (in < -0.1)
// 		{
// 		  for(unsigned int n = 0; n< n_dofs; ++n)
// 		  {
// 		    if (elem->is_node_on_side(n,side))
// 		    {
// 		      for (unsigned int nc = 0; nc < n_dofs; nc++)
// 			Ke(n,nc) = 0.0;
		      
// 		      Ke(n,n) = 1;
		      
// 		      //   std::cout<< bv[elem->get_node(n)][vec_spec]<<std::endl;
// 		      Fe(n) = bv[elem->get_node(n)][vec_spec]; //From temp to energy
// 		       //Fe(n) = e_0; //From temp to energy
		      

// 		    }
// 		  }
// 		}
		
// 		break;
		
// 	      }
// 	    case  ThermalContact::Diffusive:	      
// 	      {
// 		double temp = (dynamic_cast<Diffusive*> (contact) )->get_temperature();
// 		double emittivity = (dynamic_cast<Diffusive*> (contact) )->get_emittivity();

               
// 		double in = dir * normal[0];
// 		if (in < -0.1)
// 		{

// 		  for(unsigned int n = 0; n< n_dofs; ++n)
// 		  {
// 		    if (elem->is_node_on_side(n,side))
// 		    {
// 		      for (unsigned int nc = 0; nc < n_dofs; nc++)
// 			Ke(n,nc) = 0.0;
		      
// 		      Ke(n,n) = 1;
		      
//                       double en = temp;
// 		      double value = emittivity * en + (1.0 - emittivity) * bv[elem->get_node(n)][vec_spec];
                      
// 		      Fe(n) = value;
		     

// 		    }
// 		  }
// 		}
		
// 		break;
		
// 	      }
 
// 	    }//switch

// 	  }//  if (bd->get_boundary_properties( get_id() ) != NULL )
// 	    //}// if (is_boundary != NULL)
// 	  //	else
// 	  //	{


// 	  //  std::cout<<"ecco"<<std::endl;
// // 	  //Default Boundary Properties
// // 	  double in = dir * normal[0];
	  
// // 	  if (in < -0.1)
// // 	  {
// // 	    for(unsigned int n = 0; n< n_dofs; ++n)
// // 	    {
// // 	      if (elem->is_node_on_side(n,side))
// // 	      {
// // 		for (unsigned int nc = 0; nc < n_dofs; nc++)
// // 		  Ke(n,nc) = 0.0;
		
// // 		Ke(n,n) = 1;
		
		
// // 		Fe(n) = bv[elem->get_node(n)][vec_spec]; //From temp to energy
		
		
// // 	      }
// // 	    }
// // 	  }
			

// //	}// if (is_boundary != NULL)
// 	}
//       } //if side is on boundary
 
//     }// for (unsigned int side = 0; side<num_sides; side++)


//     dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
//     system.matrix->add_matrix (Ke, dof_indices);
//     system.rhs->add_vector    (Fe, dof_indices);


//   } //End Loop over elements
//   //   system.matrix->print_matlab("Matr.m");
//   //   system.rhs->print();
//   std::cout<<"DONE"<<std::endl;
// } //do assembly


// double
// MicroHeatBalance::print_error(void)
// {

//  //Compute Ebergy Error----------------------------
    
//     double energy_error = abs(energy_norm - old_energy_norm)/energy_norm;
    
//     //-----------------------------------------

//     double power_dissipated = 0;
//     double power_emitted = 0;
//     double dissipated_error = 0;
    
//     calculate_power_dissipated(power_dissipated,dissipated_error);
//     power_emitted = calculate_power_emitted();
      
//     if (power_emitted != 0.0)
//       dissipated_error = std::abs((power_dissipated - power_emitted)/(power_dissipated + power_emitted));

//     double  err = energy_error + dissipated_error;
  
//     if  (SimulationOptions::verbose() > 1)
//       std::cout<<"Iter:   "<< iter <<"  Error:"<<err<<std::endl;
    
//     if  (SimulationOptions::verbose() > 1)
//     {
//       std::cout<<"Error: PowerBalance:  "    <<dissipated_error <<std::endl;
//       std::cout<<"Error: Energy Error:  "<<energy_error<<std::endl;
//       std::cout<<"Error: PowerDissipated  "    <<power_dissipated <<" W"<<std::endl;
//       std::cout<<"Error: PowerEmitted:  "<<power_emitted<<" W"<<std::endl;

//       std::cout<<"PowerCheck:  "<<energy_conservation_check()<<std::endl;
//       std::cout<<" "<<std::endl;
      
//     }

//     return err;



// }





// //----------------------------------------------------------------------------------//
// void MicroHeatBalance::build_nodal_results (const std::set< std::string > &variables,
// 				     std::vector< double > &results,
// 				     std::vector< std::string > &legend)
// {


//   legend.resize(0);
//   unsigned int n_vars = 0;

//   int Temp = -1;
//   if (variables.count("LatticeTemp") ||
//       variables.count("thermal")   )
//   {
//     Temp = n_vars;
//     legend.push_back("LatticeTemp");
//     n_vars++;
//   }

//    int eq = -1;
//    if (variables.count("EqEnergy") ||
//        variables.count("thermal")   )
//    {
//      eq = n_vars;a
//      legend.push_back("EqEnergy");
//      n_vars++;
//    }

//    int J_Q = -1;
//    if (variables.count("ThermalFlux") ||
//        variables.count("thermal")   )
//      {
//        J_Q  = n_vars;
//        legend.resize(legend.size() + dim +1);

//        switch (dim)
//        {
//        case 3:
//          legend[J_Q + 2] = "J_z";
//          n_vars++;
//        case 2:
//          legend[J_Q + 1] = "J_y";
//          n_vars++;
//          legend[J_Q + dim] = "modJ";
//          n_vars++;
//        default:
//          legend[J_Q] = "J_x";
//          n_vars++;
//        }
//      }


//    int par = -1;
//    if (variables.count("partial") ||
//        variables.count("thermal")   )
//    {
//      par = n_vars;
//      for (ID k = 0;k<  AngInt.n_slices;k++)
//      {
       
//        std::string label;
//        std::ostringstream i_str;
//        i_str << "energy_density" << k;
//        legend.push_back(i_str.str());
//        n_vars++;
//      }
//    }

//     const unsigned int nn  = mesh->n_nodes();

//     results.resize(nn * n_vars,0.0);
//     legend.resize(n_vars);

//     std::vector<unsigned int> dof_indices;
//     DofMap& dof_map = my_system->get_dof_map();

//     MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
//     const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();

//     std::vector< std::map< ID, double > > temp_sol;
//     for ( ; it != end; ++it)
//     {
//       const Elem* elem = *it;

//       dof_map.dof_indices (elem, dof_indices);

//       ID subdomain = elem->subdomain_id();
//       const Material* mat = _device->get_material(subdomain);
//       HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
//       heat_model->re_init();
//       cg =  heat_model->get_lattice_thermal_capacity();


//       //Get thermal flux
//       std::vector< std::map< ID, double > > jq_solution;
//       get_solution_secure(elem,JQ_var,jq_solution);

//       for (unsigned int n = 0; n < elem->n_nodes(); n++)
//       {

// 	unsigned int id =  (elem->node(n) * n_vars) ;

// 	if (Temp != -1)
// 	{
//           double e_0 = (*equilibrium_energy)(dof_indices[n]);
//           results[id+Temp]  =  e_0;// / cg;
// 	  //results[id+Temp]  = e_0;
// 	}


// 	if (J_Q != -1)
// 	{
// 	  Tensor1 Jq(0);
// 	  Jq(1) = jq_solution[n].find(JQX)->second;
// 	  Jq(2) = jq_solution[n].find(JQY)->second;
// 	  Jq(3) = jq_solution[n].find(JQZ)->second;



// 	  switch (dim)
// 	  {
// 	    case 3:
// 	      results[id + J_Q + 2] = Jq(3);
// 	    case 2:
// 	      results[id + J_Q + 1] =  Jq(2);
// 	      results[id + J_Q + dim] = norm(Jq);
// 	    default:
// 	      results[id + J_Q ] = Jq(1);
// 	  }
// 	}


// 	if (eq != -1)
// 	{
	  
// 	  double val = (*equilibrium_energy)(dof_indices[n]);
	  
// 	  results[id+eq]  = val * cg/ AngInt.total_angle;

// 	}

//          if (par != -1)
//            for (ID k = 0;k<  AngInt.n_slices;k++)
//              results[id+par+k] = (*sol_dir[k])(dof_indices[n]) * cg / AngInt.total_angle;

	   


//       }

//     }

// }



// void
// MicroHeatBalance::calculate_power_dissipated(double& power_dissipated, double& error)
// {

//   TiberLinearSystem& system = *my_system;
//   const unsigned int  var = system.variable_number("T");
//   DofMap& dof_map =  system.get_dof_map();
//   FEType fe_type = dof_map.variable_type(var);
//   const SimulationEnvironment& env = get_environment();
//   //Surface function
//   AutoPtr<FEBase>  fe_face(build_finite_element(dim,fe_type,true));
//   QGauss qface(dim-1, CONSTANT);
//   fe_face->attach_quadrature_rule(&qface);
//   const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
//   const std::vector<Point>& qface_point = fe_face->get_xyz();
//   const std::vector<Real>& JxW_face = fe_face->get_JxW();
//   const std::vector<Point>& face_normals = fe_face->get_normals();


//  //std::vector<unsigned int> dof_indices;

//   MeshBase::const_element_iterator el =    mesh->active_local_elements_begin();
//   const MeshBase::const_element_iterator end_el =     mesh->active_local_elements_end();


//   power_dissipated = 0.0;
//   double  power_dissipated_abs = 0.0;

//   for ( ; el != end_el ; ++el)
//   {

//     const Elem* elem = *el;
//     const Elem* top_parent = (*el)->top_parent();

//     ID dof = elem->dof_number(gray_sys_number,0,0);
    
//     for (unsigned int s = 0; s < elem->n_sides(); s++)
//     {
//       ElementSide side(top_parent, s);
      
//       if (env.is_outer_boundary(side))
//       {
// 	fe_face->reinit(elem,s);

// 	RealGradient P(0);
	
// 	for (ID d = 0; d < dim; d++)
// 	  P(d) = (*thermal_flux[d])(dof);
	 
// 	double power = JxW_face[0] * (P * face_normals[0]); 

// 	power_dissipated_abs  += std::abs(power); 
// 	power_dissipated += power;
	
	
//       }
      
//     } // end loop over elem sides
//   } // end loop over elements
  
//   if (power_dissipated_abs == 0.0)
//     error = 1.0;
//   else
//     error = 2.0 * std::abs(power_dissipated/power_dissipated_abs);

//   //return power_dissipated;
//}

// double
// MicroHeatBalance::calculate_power_emitted_fourier(void)
// {

//   //Fourier system INIT
//   SimulationEnvironment& se = get_environment();

//   const NumericVector<double>& solution = *(fourier_system->solution);
//   const unsigned int var = fourier_system->variable_number("T");

//   DofMap& dof_map = fourier_system->get_dof_map();
//   FEType fe_type = dof_map.variable_type(var);
  
//   AutoPtr<FEBase>  fe(build_finite_element(dim,fe_type,true));
//   QGauss qrule(dim, FIFTH);
//   fe->attach_quadrature_rule(&qrule);

//   const std::vector<std::vector<Real> >&  phi = fe->get_phi();
//   const std::vector<Real>& JxW = fe->get_JxW();
//   const std::vector<Point>& q_point = fe->get_xyz();
 

//   double H = 0.0;

//   MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
//   const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
//   double P = 0.0;
//   for ( ; it != end; ++it)
//   {

//     const Elem* elem = *it;  
    
//     fe->reinit(elem);
    
//     ID subdomain = elem->subdomain_id();
//     const Material* mat = _device->get_material(subdomain);
//     HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
//     heat_model->re_init();
//     std::vector<double> heat_source;
//     heat_model->get_total_heat_source(elem,q_point,heat_source);

//     for (ID qp = 0; qp <  qrule.n_points(); qp++)
//       H += JxW[qp] * heat_source[qp];
    
   
//   }
   
//   H *=myopts.s_0 * myopts.s_0 * myopts.s_0;

//   return H;
 



// }

// double
// MicroHeatBalance::energy_conservation_check_fourier()
// {

//   SimulationEnvironment& se = get_environment();
  
//   DofMap& dof_map = fourier_system->get_dof_map();
//   std::vector<unsigned int> dof_indices;
  
//   const NumericVector<double>& solution = *(fourier_system->solution);
//   const unsigned int var = fourier_system->variable_number("T");
//   FEType fe_type = dof_map.variable_type(var);
  
//   AutoPtr<FEBase>  fe(build_finite_element(dim,fe_type,true));
//   QGauss qrule(dim, FIFTH);
//   fe->attach_quadrature_rule(&qrule);

//   const std::vector<Point>& q_point = fe->get_xyz();
//   const std::vector<std::vector<RealGradient> >&  dphi = fe->get_dphi();
//   const std::vector<Real>& JxW = fe->get_JxW();

//   MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
//   const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();

//   double check = 0.0;
//   for ( ; it != end; ++it)
//   {

//     const Elem* elem = *it;  
//     dof_map.dof_indices(elem, dof_indices);
//     const unsigned int n_dofs = dof_indices.size();
    
//     bool has_node = false;
//     const ID num_sides = elem->n_sides();
//     for (ID ns = 0; ns<num_sides; ns++)
//     {
//       const ElementSide elside(elem->top_parent(),ns);
//       if (se.is_on_boundary(elside)) //natural boundary
// 	has_node = true;

//     }

//     if (!has_node)
//       continue;

//     ID subdomain = elem->subdomain_id();
//     const Material* mat = _device->get_material(subdomain);
//     HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
//     heat_model->re_init();
//     vg =  heat_model->get_phonon_group_velocity();
//     tg =  heat_model->get_phonon_scattering();
//     cg =  heat_model->get_lattice_thermal_capacity();
    
//     double cg = 1.0; //j/cm-3
//     vg *=myopts.s_0;
//     cg /=myopts.s_0;
//     cg /=myopts.s_0;
//     cg /=myopts.s_0;

//     double kappa = cg * vg * vg * tg /3.0;

//     fe->reinit(elem);
    
//     for (ID alpha = 0; alpha<dof_indices.size() ;alpha ++) 
//       for (ID beta = 0; beta<dof_indices.size() ;beta ++) 
// 	for (ID qp = 0; qp <  qrule.n_points(); qp++) 
// 	  check -= kappa *  JxW[qp] * solution(dof_indices[alpha]) * (dphi[alpha][qp] * dphi[beta][qp]) ;
	   
	    
	
//   }

    
//   return check; 
 
// }



// double
// MicroHeatBalance::calculate_power_dissipated_fourier(void)
// {


//     //Fourier system INIT
//   SimulationEnvironment& se = get_environment();

//   DofMap& dof_map = fourier_system->get_dof_map();
//   std::vector<unsigned int> dof_indices;
  
//   const NumericVector<double>& solution = *(fourier_system->solution);
//   const unsigned int var = fourier_system->variable_number("T");
//   FEType fe_type = dof_map.variable_type(var);
  
//   AutoPtr<FEBase>  fe_face(build_finite_element(dim,fe_type,true));
//   QGauss qface(dim-1, FIFTH);
//   fe_face->attach_quadrature_rule(&qface);

//   const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
//   const std::vector<Real>& JxW_face = fe_face->get_JxW();
//   const std::vector<Point>& qface_point = fe_face->get_xyz();
//   const std::vector<Point>& normal = fe_face->get_normals();
//   const std::vector<std::vector<RealGradient> >&  dphi_face = fe_face->get_dphi();

//   MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
//   const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
//   double P = 0.0;
//   for ( ; it != end; ++it)
//   {

//     const Elem* elem = *it;  
//     dof_map.dof_indices(elem, dof_indices);
//     const unsigned int n_dofs = dof_indices.size();

//     ID subdomain = elem->subdomain_id();
//     const Material* mat = _device->get_material(subdomain);
//     HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
//     heat_model->re_init();
//     vg =  heat_model->get_phonon_group_velocity();
//     tg =  heat_model->get_phonon_scattering();
//      cg =  heat_model->get_lattice_thermal_capacity();
    
//     double cg = 1.0; //j/cm-3
//     vg *=myopts.s_0;
//     cg /=myopts.s_0;
//     cg /=myopts.s_0;
//     cg /=myopts.s_0;

//     double kappa = cg * vg * vg * tg /3.0;

//     const unsigned int num_sides = elem->n_sides();

//     for (ID ns = 0; ns<num_sides; ns++)
//     {
//       ElementSide side(elem->top_parent(),ns);
      
//       if (se.is_on_boundary(side))
//       {
//         fe_face->reinit(elem,ns);

// 	for (ID i = 0; i<dim ;i ++) 
// 	  for (ID alpha = 0; alpha<dof_indices.size() ;alpha ++) 
// 	    for (ID qp = 0; qp <  qface.n_points(); qp++) 
// 	    {
// 	      P -= kappa *  JxW_face[qp] * solution(dof_indices[alpha]) * dphi_face[alpha][qp](i) * normal[qp](i);
// 	      // if (normal[qp](2) == 1.0)
// 	      //		std::cout<<kappa *  JxW_face[qp] * solution(dof_indices[alpha]) * dphi_face[alpha][qp](i) * normal[qp](i)<<std::endl;
// 	    }
	
//       }

//     }
//   }
   
   


//   P *=myopts.s_0 * myopts.s_0;

//   return P;

// }

// double
// MicroHeatBalance::calculate_power_emitted(void)
// {

//   // we only do something if we are on processor 0
//   if (libMesh::processor_id() != 0)
//     return 0;

//   //Get dof map
//   TiberLinearSystem& system = *my_system;
//   const unsigned int  var = system.variable_number("T");
//   DofMap& dof_map =  system.get_dof_map();
//   FEType fe_type = dof_map.variable_type(var);

//   //Fe Build
//   AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true));
//   QGauss qrule (dim, CONSTANT);
//   fe -> attach_quadrature_rule (&qrule);
//   const std::vector<Real>& JxW = fe->get_JxW();
//   const std::vector<Point>& q_point = fe->get_xyz();
//   const std::vector<std::vector<Real> >& phi = fe->get_phi();
//   const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  
//   std::vector<unsigned int> dof_indices;
//   MeshBase::const_element_iterator el =    mesh->active_local_elements_begin();
//   const MeshBase::const_element_iterator end_el =     mesh->active_local_elements_end();
  
//   const SimulationEnvironment& env = get_environment();
//   double PowerEmitted = 0.0;

//   for ( ; el != end_el ; ++el)
//   {
    
//     const Elem* elem = *el;
    
//     fe->reinit(elem);

//     ID subdomain = elem->subdomain_id();
//     dof_map.dof_indices(elem, dof_indices,var);
//     const Material* mat = _device->get_material(subdomain);
//     HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
//     heat_model->re_init();
//     std::vector<double> heat_source;
//     heat_model->get_total_heat_source(elem,q_point,heat_source);

//     PowerEmitted  += JxW[0] * heat_source[0]/myopts.s_0 /myopts.s_0 /myopts.s_0   ;

    
//   } // end loop over elements


//   return  PowerEmitted;
// }

