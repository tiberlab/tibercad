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

  myopts.integration_order = static_cast<libMeshEnums::Order>(
      opts.get_option("integration_order", 5));

  myopts.work_units = opts.get_option("Work_length_units", 1e-2);

  myopts.scale = opts.get_option("scale", 1.0);

  myopts.max_error =  opts.get_option("max_error",1e-2);

  myopts.max_iter =  opts.get_option("max_iter",1);

  //myopts.ref_temp =  opts.get_option("reference_temperature",SimulationOptions::temperature);


  AngInt.theta_slices = opts.get_option("theta_slices",1);
  AngInt.phi_slices =    opts.get_option("phi_slices",1);
  AngInt.dim = opts.get_option("dim",dim);
  
  myopts.alternative = opts.get_option("alternative",0);

  //  myopts.eq_temp = opts.get_option("equilibrium_temperature",myopts.ref_temp);

  
  myopts.equilibrium_energy = opts.get_option("equilibrium_energy",0.0);
  
  myopts.first_guess = opts.get_option("first_guess","fourier");


}

void MicroHeatBalance::do_init( )
{

 parse_options();
  

  const ModelOptions& sim_opt = get_options();

  SimulationEnvironment& si = get_environment();

  _device = &( si.get_device() );

  mesh = &(_device->get_mesh());

 //  mesh = & (_device->get_mesh());

  dim = mesh->mesh_dimension();

  double mesh_units = get_scaling().get_calc_mesh_units() / sim_opt.get_option("Work_length_units", 1e-2);
  //double mesh_units = get_scaling().get_calc_mesh_units() / sim_opt.get_option("Work_length_units", 1e-9);

  t_0 =1;//e12; //1e12
  s_0 =1;//e7;


   get_scaling().set_calc_mesh_units(mesh_units);


  my_system = TiberLinearSystem::create(get_equation_systems(),
      get_equation_system_name(), get_solver_options());


  my_system->add_variable("T", FIRST);

   // Insert the pointer to function that LibMesh library has to use
  my_system->attach_assemble_function (assemble_heat_matrix);

 my_system->set_options(get_solver_options());
   // Initialize the data structures for the equation system.
  my_system->init();


  //Inizialize the solution to temperature of simulation options
  

  my_system->solution->zero();
  my_system->solution->close();

  heat_legend = "Wq";

  JQ_var.insert(JQX);
  JQ_var.insert(JQY);
  JQ_var.insert(JQZ);

  e0_var.insert(E0);

  if (myopts.alternative == 1)
    AngInt.compute_very_alternative_directions();
  else
    AngInt.compute_alternative_directions();

  thermal_flux.resize(dim);
  for (ID i = 0; i<dim; i++ )
    thermal_flux[i] = (my_system->solution)->clone().release();

 
  sol_dir.resize(AngInt.n_slices);
  for (ID k = 0; k<AngInt.n_slices ; k++ )
    sol_dir[k] = (my_system->solution)->clone().release();
  

  equilibrium_energy = (my_system->solution)->clone().release();
  equilibrium_energy->add(0.1); //Just put some dummy value
  equilibrium_energy->close();

  equilibrium_energy_new = (my_system->solution)->clone().release();
  equilibrium_energy_new->add(0.1); //Just put some dummy value
  equilibrium_energy_new->close();
 

  //Initialize Boundary Values for Specular Boundary Conditions
  
  {  
    SimulationEnvironment& env = get_environment();
    
    MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh->active_elements_end();
    
    for ( ; el != end_el ; ++el)   //loop over elements
    {
      const Elem* elem = *el;
      const unsigned int num_sides = elem->n_sides();
      for (unsigned int side = 0; side<num_sides; side++)
      {
	const ElementSide elside(elem->top_parent(), side);
	
	if (env.is_on_boundary(elside))
	{
	  for(unsigned int n = 0; n< elem->n_nodes(); ++n)
	    if (elem->is_node_on_side(n,side))
	      bv[elem->get_node(n)].resize(AngInt.n_slices,0.0);
	    

	}
	
      } 
      
    }
  }



  
  //------init is done---------------------------------------------------------------------//
  
}
//-------------------------------------------------------------------------------//
void  MicroHeatBalance::do_solve()
{
  


  std::cout<<"scale: "<<myopts.scale<<std::endl;

  static_this = this;

 

  if (~is_solved())
  { 
        
    equilibrium_energy->zero();
   
  

    if (myopts.first_guess == "fourier")
      compute_fourier_solution();
    else
      equilibrium_energy->add(myopts.equilibrium_energy);
      
    //equilibrium_energy->scale(1/ AngInt.total_angle);
    
   
  }

  //------------------------------------

  //Fill the directional results
  for (ID k = 0; k<AngInt.n_slices ; k++ )
  {
    sol_dir[k]->zero();
    sol_dir[k]->add(*equilibrium_energy);
  }


 //  //Fill the boundary value
  {
    BoundaryData::iterator it(bv.begin());
    const  BoundaryData::iterator end(bv.end());
    
    const unsigned int system_number = my_system->number();
    const unsigned int var = my_system->variable_number("T");
    
    for ( ; it != end; ++it)
    {
      const unsigned int  n_dof = (it->first)->dof_number(system_number,var,0);
      for (ID k = 0; k<AngInt.n_slices ; k++ )
	(it->second)[k] = (*sol_dir[k])(n_dof);
      
    }
  }
                  	
  
  //  //Error
  ID iter = 0;
  double err1 = 0.0;
  double err2 = 0.0;
  double err3 = 0.0;
  double err = 0.0;
  equilibrium_energy->close();
  double old_norm = equilibrium_energy->l2_norm();

  
  if  (myopts.max_iter == 0)
  {
  }
  else 
  {
  
    do {

      
      equilibrium_energy_new->zero();
     
      for (ID i = 0; i<dim ; i++ )
	thermal_flux[i]->zero();
       

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
      //(my_system->solution)->scale(1.0/myopts.scale);

      //Equilibrium---------------------
      NumericVector<Number>* temp2 = (my_system->solution)->clone().release();
      temp2->scale(d_omega);
      //equilibrium_energy_new->add(*(my_system->solution));
      equilibrium_energy_new->add(*temp2);
      delete temp2;
      //------------------------------------

      sol_dir[k]->zero();
      sol_dir[k]->add(*(my_system->solution));
     
      for (ID i = 0; i<dim ; i++ )
      {
	NumericVector<Number>* temp = (my_system->solution)->clone().release();
        temp->scale(IntDir(i));
	thermal_flux[i]->add(*temp);
	delete temp;

      }
     
      
    }
  
    equilibrium_energy_new->scale(1.0/AngInt.total_angle);
    equilibrium_energy->zero();
    equilibrium_energy->add(*equilibrium_energy_new);
    //    eq_mean = (equilibrium_energy->sum()/equilibrium_energy->size());

 //CORE--------------------Loop over direction---------------
   //  //   //Fill the boundary value
    {
      BoundaryData::iterator it(bv.begin());
      const  BoundaryData::iterator end(bv.end());
      
      const unsigned int system_number = my_system->number();
      const unsigned int var = my_system->variable_number("T");
      
      for ( ; it != end; ++it)
      {
	
	const unsigned int  n_dof = (it->first)->dof_number(system_number,var,0);
	for (ID k = 0; k<AngInt.n_slices ; k++ ) 
	  (it->second)[k] = (*sol_dir[k])(n_dof);
	
      }
    }
 
    //-----------------------------------------------------------

    
    //Compute Error----------------------------
     double sol_norm = equilibrium_energy->l2_norm();
     double err1 = abs(sol_norm - old_norm)/sol_norm;
     old_norm = sol_norm;
     double power_dissipated = 0;
     double power_emitted = 0;
     double err_diss = 0;
     calculate_power_dissipated(power_dissipated,err_diss);
     power_emitted = calculate_power_emitted();
     

     if (power_emitted == 0.0)
     {
       err2 = err_diss;
       
     }
     else
       err2 = std::abs((power_dissipated - power_emitted)/(power_dissipated + power_emitted));

     err = err1 + err2;
  
     if  (SimulationOptions::verbose() > 1)
       std::cout<<"Iter:   "<< iter <<"  Error:"<<err<<std::endl;
     
     if  (SimulationOptions::verbose() > 1)
     {
       std::cout<<"Error: PowerBalance:  "    <<err2 <<std::endl;
       std::cout<<"Error: Energy Error:  "<<err1<<std::endl;
       std::cout<<"Error: PowerDissipated  "    <<power_dissipated <<std::endl;
       std::cout<<"Error: PowerEmitted:  "<<power_emitted<<std::endl;
       std::cout<<" "<<std::endl;

     }
     
     iter +=1;
    
  } while (err > myopts.max_error & iter < myopts.max_iter);
  }
    


  

}



//--------------------------------------------------------------------------------//
MicroHeatBalance::~MicroHeatBalance()
{
  //   //Release pointers
  for (ID i = 0; i<dim; i++ )
    delete thermal_flux[i];
  
  for (ID k = 0; k<AngInt.n_slices ; k++ )
    delete sol_dir[k];


  delete  equilibrium_energy;

 delete  equilibrium_energy_new;

}
//---------------------------------------------------------------------------------//
MicroHeatBalance::MicroHeatBalance(const ModelOptions& options)
  : SimulationInterface(options)
{
 
 
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

  QGauss qrule (dim, FIFTH); //may be could be decreased (CHECK!!!)

  // quadrature rule
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

  AutoPtr<FEBase>  fe_face(build_finite_element(dim,fe_type,true));

  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // o cout<<"Start loop over lattice thermal conductivity"<<endl;f the element.map
  QGauss qface(dim-1, SIXTH);

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
 
  std::cout<<"t_0: "<<t_0<<std::endl;
  std::cout<<"s_0: "<<s_0<<std::endl;

  std::vector< std::map< ID, double > > temp_sol;

  for ( ; el != end_el ; ++el)   //loop over elements
  {

    const Elem* elem = *el;

    dof_map.dof_indices (elem, dof_indices);

    const unsigned int n_dofs = dof_indices.size();

    double e_0 = 0.0;

    fe->reinit(elem);

    Ke.resize(n_dofs,n_dofs);

    Fe.resize(n_dofs);

    Fe.zero();

    Ke.zero();

    const unsigned int num_sides = elem->n_sides();

    ID subdomain = elem->subdomain_id();

    const Material* mat = _device->get_material(subdomain);

    HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
 
    //HeatModel* heat_model = (  dynamic_cast<HeatModel*> ( get_physical_model(subdomain)  ) );

    // heat_model->set_element(elem);

    //heat_model->set_side(-1);

    heat_model->re_init();

    std::vector<double> heat_source;
    heat_model->get_total_heat_source(elem,q_point,heat_source);

    vg =  heat_model->get_phonon_group_velocity();
    tg =  heat_model->get_phonon_scattering();
    cg =  heat_model->get_lattice_thermal_capacity(); //j/cm3K

   
    vg *=s_0;
    vg /=t_0;
    tg *=t_0;
    
    std::vector<Point> p(1);
    p[0] = elem->centroid();


    get_solution_secure(elem,q_point,e0_var,temp_sol);
    //   get_solution_secure(elem,p,e0_var,temp_sol);
    //e_0 = temp_sol[0].find(E0)->second;

    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {//Loop over quadrature points
      
      e_0 = temp_sol[qp].find(E0)->second;

      
      //  std::cout<<e_0<<std::endl;
      for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
      { // loop over test function
	
	for (unsigned int p2=0; p2<n_dofs; p2++)
	{//loop over basis functions
	  
	  double value1 = 0.0;
	  double value2 = 0.0;
	
	  value1 = vg * JxW[qp] * (dphi[p2][qp] * IntDir) * phi[p1][qp];
	 
	  value2 = 1.0/tg * JxW[qp] *  phi[p1][qp] * phi[p2][qp] * d_omega;

 	  Ke(p1,p2) += value1;
	  Ke(p1,p2) += value2;

	  //  std::cout<<"v1: "<<value1<<std::endl;
	  // std::cout<<"v2: "<<value2<<std::endl;


	} //loop over basis functions

	double value1 = 0.0;
	double value2 = 0.0;

	value1 =  1.0/tg * JxW[qp] * e_0 * phi[p1][qp]  * d_omega;
	//std::cout<<"v3: "<<value1<<std::endl;
 
	value2 =  JxW[qp]  *  heat_source[qp] * phi[p1][qp]/t_0 * d_omega;// /AngInt.total_angle;

	Fe(p1) += value1;
 	Fe(p1) += value2;

      }//end Loop over quadrature points
    } // end loop over test function


    

    
    //Boundary conditions and source
    //The loop over element is the only loop that is surviving at this point

    for (unsigned int side = 0; side<num_sides; side++)
    {

      const ElementSide elside(elem->top_parent(), side);

      if (se.is_on_boundary(elside))
      {
	heat_model->re_init();    
	fe_face->reinit(elem,side);

	{
	  //Default Boundary Properties-----
	  double in = dir * normal[0];
	  
	  if (in <0.0)
	  {
	    for(unsigned int n = 0; n< n_dofs; ++n)
	    {
	      if (elem->is_node_on_side(n,side))
	      {
		//for (unsigned int nc = 0; nc < n_dofs; nc++)
		//  Ke(n,nc) = 0.0;
		
		//Ke(n,n) = 1;
		
		Ke(n,n) += 1e36;

		const unsigned int  n_dof = (elem->get_node(n))->dof_number(system_number,var,0); 
		double value =  bv[elem->get_node(n)][vec_spec];
	
	        value = e_0;
		Fe(n) += value *1e36; //From temp to energy	
		
	      }
	    }
	    
	  }
	}
        //---------------------------------


	Boundary* bd = se.get_boundary(elside);
      
	if (bd != NULL)
	{
	  if (bd->get_boundary_properties( get_id() ) != NULL )

	  {
	    
	    //heat_model->set_side(side);
	    // heat_model->re_init();
	    // fe_face->reinit(elem,side);
	    
	    ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );
	    
	    switch (contact->get_type())
	    {
	    case  ThermalContact::Reservoir:
	      {
		
		double in = dir * normal[0];
		double temp = (dynamic_cast<Reservoir*> (contact) )->get_temperature();
		
		
		if (in < 0.0)
	       {
		 for(unsigned int n = 0; n< n_dofs; ++n)
		 {
		   if (elem->is_node_on_side(n,side))
		   {
		  
		       Ke(n,n) += 1e36;
		       Fe(n) += 1e36* temp;
		       
		   		  
		     //  for (unsigned int nc = 0; nc < n_dofs; nc++)
		     // Ke(n,nc) = 0.0;
		     
		     // Ke(n,n) = 1;
		     
		    
		       //Fe(n) = temp * cg /AngInt.total_angle * myopts.scale; //From temp to energy

		     //		     Fe(n) = temp;// * cg /AngInt.total_angle * myopts.scale;
		    

  
		   }
                 }
	       }
	      }
	      break;
	      
	    case  ThermalContact::FourierBTE:
	      {
		
		double in = dir * normal[0];
                     
		if (in < 0.0)
		{
		  
		  for(unsigned int n = 0; n< n_dofs; ++n)
		  {
		    if (elem->is_node_on_side(n,side))
		    {
		      
		      const Elem* neighbour = (elside.first)->neighbor(side);
		      double temp = (dynamic_cast<FourierBTE*> (contact) )->get_temperature(neighbour,elem->point(n));
		      
		      for (unsigned int nc = 0; nc < n_dofs; nc++)
			Ke(n,nc) = 0.0;
		      
		      Ke(n,n) = 1;
		      
		     
		      Fe(n) =  temp;
		      
		    }
		  }
		}
		
		break;
	      }
	      
	     case  ThermalContact::Specular:	      
	      {
		
	
		double in = dir * normal[0];
				//	std::cout<<in<<std::endl;
		if (in < -0.1)
		{
		  for(unsigned int n = 0; n< n_dofs; ++n)
		  {
		    if (elem->is_node_on_side(n,side))
		    {
		      for (unsigned int nc = 0; nc < n_dofs; nc++)
			Ke(n,nc) = 0.0;
		      
		      Ke(n,n) = 1;
		      
		      //   std::cout<< bv[elem->get_node(n)][vec_spec]<<std::endl;
		      Fe(n) = bv[elem->get_node(n)][vec_spec]; //From temp to energy
		       //Fe(n) = e_0; //From temp to energy
		      

		    }
		  }
		}
		
		break;
		
	      }
	    case  ThermalContact::Diffusive:	      
	      {
		double temp = (dynamic_cast<Diffusive*> (contact) )->get_temperature();
		double emittivity = (dynamic_cast<Diffusive*> (contact) )->get_emittivity();

               
		double in = dir * normal[0];
		if (in < -0.1)
		{

		  for(unsigned int n = 0; n< n_dofs; ++n)
		  {
		    if (elem->is_node_on_side(n,side))
		    {
		      for (unsigned int nc = 0; nc < n_dofs; nc++)
			Ke(n,nc) = 0.0;
		      
		      Ke(n,n) = 1;
		      
                      double en = temp;
		      double value = emittivity * en + (1.0 - emittivity) * bv[elem->get_node(n)][vec_spec];
                      
		      Fe(n) = value;
		     

		    }
		  }
		}
		
		break;
		
	      }
 
	    }//switch

	  }//  if (bd->get_boundary_properties( get_id() ) != NULL )
	    //}// if (is_boundary != NULL)
	  //	else
	  //	{


	  //  std::cout<<"ecco"<<std::endl;
// 	  //Default Boundary Properties
// 	  double in = dir * normal[0];
	  
// 	  if (in < -0.1)
// 	  {
// 	    for(unsigned int n = 0; n< n_dofs; ++n)
// 	    {
// 	      if (elem->is_node_on_side(n,side))
// 	      {
// 		for (unsigned int nc = 0; nc < n_dofs; nc++)
// 		  Ke(n,nc) = 0.0;
		
// 		Ke(n,n) = 1;
		
		
// 		Fe(n) = bv[elem->get_node(n)][vec_spec]; //From temp to energy
		
		
// 	      }
// 	    }
// 	  }
			

//	}// if (is_boundary != NULL)
	}
      } //if side is on boundary
 
    }// for (unsigned int side = 0; side<num_sides; side++)


    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices);


  } //End Loop over elements
  //   system.matrix->print_matlab("Matr.m");
  //   system.rhs->print();
  std::cout<<"DONE"<<std::endl;
} //do assembly




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
	 Jq(i) += vg * (*thermal_flux[i])(dof_indices[alpha]) * phi[alpha][n]/t_0 /AngInt.total_angle;
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







//----------------------------------------------------------------------------------//
void MicroHeatBalance::build_nodal_results (const std::set< std::string > &variables,
				     std::vector< double > &results,
				     std::vector< std::string > &legend)
{


  legend.resize(0);
  unsigned int n_vars = 0;

  int Temp = -1;
  if (variables.count("LatticeTemp") ||
      variables.count("thermal")   )
  {
    Temp = n_vars;
    legend.push_back("LatticeTemp");
    n_vars++;
  }

   int eq = -1;
   if (variables.count("EqEnergy") ||
       variables.count("thermal")   )
   {
     eq = n_vars;
     legend.push_back("EqEnergy");
     n_vars++;
   }

   int J_Q = -1;
   if (variables.count("ThermalFlux") ||
       variables.count("thermal")   )
     {
       J_Q  = n_vars;
       legend.resize(legend.size() + dim +1);

       switch (dim)
       {
       case 3:
         legend[J_Q + 2] = "J_z";
         n_vars++;
       case 2:
         legend[J_Q + 1] = "J_y";
         n_vars++;
         legend[J_Q + dim] = "modJ";
         n_vars++;
       default:
         legend[J_Q] = "J_x";
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

    const unsigned int nn  = mesh->n_nodes();

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

      ID subdomain = elem->subdomain_id();
      const Material* mat = _device->get_material(subdomain);
      HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
      heat_model->re_init();
      cg =  heat_model->get_lattice_thermal_capacity();


      //Get thermal flux
      std::vector< std::map< ID, double > > jq_solution;
      get_solution_secure(elem,JQ_var,jq_solution);

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {

	unsigned int id =  (elem->node(n) * n_vars) ;

	if (Temp != -1)
	{
          double e_0 = (*equilibrium_energy)(dof_indices[n]);
          results[id+Temp]  =  e_0;// / cg;
	  //results[id+Temp]  = e_0;
	}


	if (J_Q != -1)
	{
	  Tensor1 Jq(0);
	  Jq(1) = jq_solution[n].find(JQX)->second;
	  Jq(2) = jq_solution[n].find(JQY)->second;
	  Jq(3) = jq_solution[n].find(JQZ)->second;



	  switch (dim)
	  {
	    case 3:
	      results[id + J_Q + 2] = Jq(3);
	    case 2:
	      results[id + J_Q + 1] =  Jq(2);
	      results[id + J_Q + dim] = norm(Jq);
	    default:
	      results[id + J_Q ] = Jq(1);
	  }
	}


	if (eq != -1)
	{
	  
	  double val = (*equilibrium_energy)(dof_indices[n]);
	  
	  results[id+eq]  = val / AngInt.total_angle;

	}

         if (par != -1)
           for (ID k = 0;k<  AngInt.n_slices;k++)
	   {
	  
             results[id+par+k] = (*sol_dir[k])(dof_indices[n]);

	   }


      }

    }

}
void
MicroHeatBalance::build_integrated_quantities_description(
    const std::set<std::string>& names,
    std::vector<std::string>& legend,
    std::vector<std::string>& description)
{


  if (names.count("PowerDissipated"))
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
MicroHeatBalance::build_integrated_quantities(const set<string>& names,
    vector<double>& values)
{


  if (names.count("PowerDissipated"))
  {

    // double power = calculate_power_dissipated();

    //  values.resize(1,power);


  }
}



void
MicroHeatBalance::calculate_power_dissipated(double& power_dissipated, double& error)
{

  // we only do something if we are on processor 0
  //  if (libMesh::processor_id() != 0)
    //    return 0;

  TiberLinearSystem& system = *my_system;
  const unsigned int  var = system.variable_number("T");
  DofMap& dof_map =  system.get_dof_map();
  FEType fe_type = dof_map.variable_type(var);
  const SimulationEnvironment& env = get_environment();
  //Surface function
  AutoPtr<FEBase>  fe_face(build_finite_element(dim,fe_type,true));
  QGauss qface(dim-1, SIXTH);
  fe_face->attach_quadrature_rule(&qface);
  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
  const std::vector<Point>& qface_point = fe_face->get_xyz();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& face_normals = fe_face->get_normals();


 //std::vector<unsigned int> dof_indices;

  MeshBase::const_element_iterator el =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =     mesh->active_local_elements_end();

//  const SimulationEnvironment& env = get_environment();

  power_dissipated = 0.0;
  double  power_dissipated_abs = 0.0;

  for ( ; el != end_el ; ++el)
  {

    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      ElementSide side(top_parent, s);

      if (env.is_outer_boundary(side))
      {

	fe_face->reinit(elem, s);
	std::vector< std::map< ID, double > > jq_solution;

        get_solution_secure(elem,qface_point,JQ_var,jq_solution);

	RealGradient P(0);

	for (unsigned int qp = 0; qp <  qface.n_points(); qp++)
	{
	
	  P(0) = jq_solution[qp].find(JQX)->second;
	  P(1) = jq_solution[qp].find(JQY)->second;
	  P(2) = jq_solution[qp].find(JQZ)->second;

	 
          double power = JxW_face[qp] * (P * face_normals[qp]); 

          power_dissipated_abs  += std::abs(power); 
	  power_dissipated += power;


	}
      }
    } // end loop over elem sides
  } // end loop over elements

  if (power_dissipated_abs == 0.0)
    error = 1.0;
  else
    error = 2.0 * std::abs(power_dissipated/power_dissipated_abs);

  // return power_dissipated;
}




void
MicroHeatBalance::do_print_info(void)
{

  string space("  ");
  //cout << space << "linear solver is: petsc" <<std::endl;

}

void
MicroHeatBalance::AngularIntegrator::compute_directions(void)
{

  n_slices = theta_slices * phi_slices;
  directions.resize(n_slices);
  integrate_directions.resize(n_slices);
  d_omega.resize(n_slices);
  dir.resize(n_slices);
  theta_vec.resize(n_slices);
  phi_vec.resize(n_slices);


  double min_theta, max_theta, min_phi, max_phi;


  switch (dim)
  {

  case 1 :


    min_theta = 0.0;
    max_theta = M_PI * 0.5;

    min_phi = M_PI * 0.5;
    max_phi = M_PI * 1.5;

    weight = 4.0;

    break;

  case 2 :


    min_theta = 0.0;
    max_theta = M_PI * 0.5;

    min_phi = 0;
    max_phi = M_PI * 2.0;

    weight = 2.0;

    break;

  case 3 :


    min_theta = 0.0;
    max_theta = M_PI;

    min_phi = 0.0;
    max_phi = M_PI * 2.0;

    weight = 1.0;


    break;
  }


  d_theta =  (max_theta - min_theta) / theta_slices;
  d_phi =  (max_phi - min_phi) / phi_slices;
  double theta, phi;



  ID k = 0;
  for (ID n_phi = 0; n_phi < phi_slices; n_phi++)
  {
    phi = min_phi + d_phi * 0.5 + d_phi * n_phi;

    for (ID n_theta = 0; n_theta < theta_slices; n_theta++)
    {
      theta = min_theta + d_theta * 0.5 + d_theta * n_theta;


      d_omega[k] = weight * 2.0 * sin(theta) * sin (0.5 * d_theta) * d_phi;


      dir[k](0) = sin(theta) * sin(phi);
      dir[k](1) = sin(theta) * cos(phi);
      dir[k](2) = cos(theta);


      directions[k](0) = weight * sin(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
      directions[k](1) = weight * cos(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
      directions[k](2) = weight * 0.5 * d_phi * sin(2.0 * theta) * sin(d_phi);

      theta_vec[k] = theta;
      phi_vec[k] = phi;

      k++;
    }
  }


}


void
MicroHeatBalance::AngularIntegrator::compute_alternative_directions(void)
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

    weight = 1.0;



    break;
  case 2 :

    theta_slices = 1;
    //phi_slices = 2;

    min_theta = 0.0;
    max_theta = M_PI;

    min_phi = M_PI * 0.5;
    max_phi = M_PI * 2.0 + M_PI * 0.5;

    weight = 1.0;


    break;
  case 3 :


    min_theta = 0.0;
    max_theta = M_PI;

    min_phi = 0.0;
    max_phi = M_PI * 2.0;

    weight = 1.0;


    break;
  }


  total_angle = 4.0 * M_PI;
  n_slices = theta_slices * phi_slices;
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


      d_omega[k] = weight * 2.0 * sin(theta) * sin (0.5 * d_theta) * d_phi;


      dir[k](0) = sin(theta) * sin(phi);
      dir[k](1) = sin(theta) * cos(phi);
      dir[k](2) = cos(theta);

      //directions[k] = dir[k];
      directions[k](0) = weight * sin(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
      directions[k](1) = weight * cos(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
      directions[k](2) = weight * 0.5 * d_phi * sin(2.0 * theta) * sin(d_theta);

      theta_vec[k] = theta;
      phi_vec[k] = phi;

      k++;
    }
  }

  spec.resize(n_slices);
  spec[0] = 4;
  spec[1] = 5;
  spec[2] = 6;
  spec[3] = 7;
  spec[4] = 0;
  spec[5] = 1;
  spec[6] = 2;
  spec[7] = 3;
}


void
MicroHeatBalance::AngularIntegrator::compute_very_alternative_directions(void)
{
  std::cout<<"Dimension:   "<<dim<<std::endl;
   switch (dim)
  {

  case 1 :
    
    n_slices = 2;
    directions.resize(n_slices);
    integrate_directions.resize(n_slices);
    d_omega.resize(n_slices);
    dir.resize(n_slices);
    theta_vec.resize(n_slices);
    phi_vec.resize(n_slices);
    spec.resize(n_slices);    

    total_angle = 2;
    
    d_omega[0] = 1;
    d_omega[1] = 1;
    
    dir[0](0) = 0;
    dir[0](1) = 0;
    dir[0](2) = 1;
    
    dir[1](0) = 0;
    dir[1](1) = 0;
    dir[1](2) = -1;

    directions = dir;

    spec[0] = 1;
    spec[1] = 0;
    
    
    break;
  case 2:
   
    n_slices = 4;
  
    total_angle = 4;
   
    directions.resize(n_slices);
    integrate_directions.resize(n_slices);
    d_omega.resize(n_slices);
    dir.resize(n_slices);
    theta_vec.resize(n_slices);
    phi_vec.resize(n_slices);
    spec.resize(n_slices);  


    
    d_omega[0] = 1;
    d_omega[1] = 1;
    d_omega[2] = 1;
    d_omega[3] = 1;

    dir[0](0) = 1;
    dir[0](1) = 0;
    dir[0](2) = 0;

    dir[1](0) = -1;
    dir[1](1) = 0;
    dir[1](2) = 0;

    dir[2](0) = 0;
    dir[2](1) = 1;
    dir[2](2) = 0;

    dir[3](0) = 0;
    dir[3](1) = -1;
    dir[3](2) = 0;


    directions[0](0) = 1;
    directions[0](1) = 0;
    directions[0](2) = 0;
    
    directions[1](0) = -1;
    directions[1](1) = 0;
    directions[1](2) = 0;

    directions[2](0) = 0;
    directions[2](1) = 1;
    directions[2](2) = 0;
    
    directions[3](0) = 0;
    directions[3](1) = -1;
    directions[3](2) = 0;

    spec[0] = 1;
    spec[1] = 0;
    spec[2] = 3;
    spec[3] = 2;

    break;
    
  case 3:
   
    n_slices = 6;
    //total_angle = 2.0 * M_PI ;
    total_angle = 6.0;
   
    directions.resize(n_slices);
    integrate_directions.resize(n_slices);
    d_omega.resize(n_slices);
    dir.resize(n_slices);
    theta_vec.resize(n_slices);
    phi_vec.resize(n_slices);
    spec.resize(n_slices);  

    // d_omega.resize(n_slices,1.0/n_slices);
   
    double value = 1.0;
    d_omega[0] = value;
    d_omega[1] = value;
    d_omega[2] = value;
    d_omega[3] = value;
    d_omega[4] = value;
    d_omega[5] = value;
   

   
    dir[0](0) = 1;
    dir[0](1) = 0;
    dir[0](2) = 0;

    dir[1](0) = -1;
    dir[1](1) = 0;
    dir[1](2) = 0;

    dir[2](0) = 0;
    dir[2](1) = 1;
    dir[2](2) = 0;

    dir[3](0) = 0;
    dir[3](1) = -1;
    dir[3](2) = 0;

    dir[4](0) = 0;
    dir[4](1) = 0;
    dir[4](2) = 1;

    dir[5](0) = 0;
    dir[5](1) = 0;
    dir[5](2) = -1;

    directions = dir;
    
 
    spec[0] = 1;
    spec[1] = 0;
    spec[2] = 3;
    spec[3] = 2;
    spec[4] = 5;
    spec[5] = 4;

    break;
  }


  


}



void
MicroHeatBalance::AngularIntegrator::compute_very_alternative_directions2(void)
{
  std::cout<<"Dimension:   "<<dim<<std::endl;
   switch (dim)
  {

  case 1 :
    
    n_slices = 2;
    directions.resize(n_slices);
    integrate_directions.resize(n_slices);
    d_omega.resize(n_slices);
    dir.resize(n_slices);
    theta_vec.resize(n_slices);
    phi_vec.resize(n_slices);
    spec.resize(n_slices);    

    total_angle = 2;
    
    d_omega[0] = 0.5;
    d_omega[1] = 0.5;
    
    dir[0](0) = 0;
    dir[0](1) = 0;
    dir[0](2) = 1;
    
    dir[1](0) = 0;
    dir[1](1) = 0;
    dir[1](2) = -1;

    directions = dir;

    spec[0] = 1;
    spec[1] = 0;
    
    
    break;

  case 2:
   
    n_slices = 4;
    total_angle = 2.0 * M_PI ;

    //total_angle = 4;
   
    directions.resize(n_slices);
    integrate_directions.resize(n_slices);
    d_omega.resize(n_slices);
    dir.resize(n_slices);
    theta_vec.resize(n_slices);
    phi_vec.resize(n_slices);
    spec.resize(n_slices);  

    d_omega.resize(n_slices,0.0);
    
    d_omega[0] = total_angle/n_slices;
    d_omega[1] = total_angle/n_slices;
    d_omega[2] = total_angle/n_slices;
    d_omega[3] = total_angle/n_slices;

    dir[0](0) = 1;
    dir[0](1) = 0;
    dir[0](2) = 0;

    dir[1](0) = -1;
    dir[1](1) = 0;
    dir[1](2) = 0;

    dir[2](0) = 0;
    dir[2](1) = 1;
    dir[2](2) = 0;

    dir[3](0) = 0;
    dir[3](1) = -1;
    dir[3](2) = 0;

    //  directions = dir * 4.0  * M_PI;

    //Real ange = 1.0 * M_PI;

    directions[0](0) = M_PI;
    directions[0](1) = 0;
    directions[0](2) = 0;
    
    directions[1](0) = -M_PI;
    directions[1](1) = 0;
    directions[1](2) = 0;

    directions[2](0) = 0;
    directions[2](1) = M_PI;
    directions[2](2) = 0;
    
    directions[3](0) = 0;
    directions[3](1) = -M_PI;
    directions[3](2) = 0;

    spec[0] = 1;
    spec[1] = 0;
    spec[2] = 3;
    spec[3] = 2;

    break;
    
  case 3:
   
    n_slices = 6;
    //total_angle = 2.0 * M_PI ;
    total_angle = 4.0 * M_PI;
   
    directions.resize(n_slices);
    integrate_directions.resize(n_slices);
    d_omega.resize(n_slices);
    dir.resize(n_slices);
    theta_vec.resize(n_slices);
    phi_vec.resize(n_slices);
    spec.resize(n_slices);  

    //d_omega.resize(n_slices,1.0/n_slices);


    d_omega[0] = n_slices/total_angle;
    d_omega[1] = n_slices/total_angle;
    d_omega[2] = n_slices/total_angle;
    d_omega[3] = n_slices/total_angle;
    d_omega[4] = n_slices/total_angle;
    d_omega[5] = n_slices/total_angle;

    dir[0](0) = 1;
    dir[0](1) = 0;
    dir[0](2) = 0;

    dir[1](0) = -1;
    dir[1](1) = 0;
    dir[1](2) = 0;

    dir[2](0) = 0;
    dir[2](1) = 1;
    dir[2](2) = 0;

    dir[3](0) = 0;
    dir[3](1) = -1;
    dir[3](2) = 0;

    dir[4](0) = 0;
    dir[4](1) = 0;
    dir[4](2) = 1;

    dir[5](0) = 0;
    dir[5](1) = 0;
    dir[5](2) = -1;

    directions = dir;
    
 
    spec[0] = 1;
    spec[1] = 0;
    spec[2] = 3;
    spec[3] = 2;
    spec[4] = 5;
    spec[5] = 4;

    break;
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

    std::cout<<"Direction: "<<k<<std::endl;
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

  fourier_system = TiberLinearSystem::create(get_equation_systems(),
					     "fourier", get_solver_options());

 

  fourier_system->add_variable("T", FIRST);

 
  
  fourier_system->attach_assemble_function(assemble_macro_heat_matrix);

  
  
  fourier_system->init();


  fourier_system->solve();



  //Convert from Temperature to energy
  
  



  equilibrium_energy->add(*(fourier_system->solution));

   


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
    double kappa = vg * vg * tg /3.0;

    //s_0 = 1.0;
    //
kappa /=s_0;
  					
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
	Fe(p1) +=JxW[qp] * heat_source[qp] * phi[p1][qp] /s_0/s_0/s_0;
    

     
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
		  Fe(n) = (dynamic_cast<Reservoir*> (contact) )->get_temperature();
		 
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
		  const Elem* neighbour = (elside.first)->neighbor(side);
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


double
MicroHeatBalance::calculate_power_emitted(void)
{

  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return 0;

  //Get dof map
  TiberLinearSystem& system = *my_system;
  const unsigned int  var = system.variable_number("T");
  DofMap& dof_map =  system.get_dof_map();
  FEType fe_type = dof_map.variable_type(var);

 // const unsigned int dim = mesh->mesh_dimension();


  //Fe Build
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true));
  QGauss qrule (dim, FIFTH);
  fe -> attach_quadrature_rule (&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  std::vector<unsigned int> dof_indices;
  MeshBase::const_element_iterator el =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =     mesh->active_local_elements_end();

  const SimulationEnvironment& env = get_environment();
  double PowerEmitted = 0.0;

  for ( ; el != end_el ; ++el)
  {

    const Elem* elem = *el;

    fe->reinit(elem);

    ID subdomain = elem->subdomain_id();
    dof_map.dof_indices(elem, dof_indices,var);
    const Material* mat = _device->get_material(subdomain);
    HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
    //    heat_model->set_element(elem);
    // heat_model->set_side(-1);
    heat_model->re_init();
    std::vector<double> heat_source;
    heat_model->get_total_heat_source(elem,q_point,heat_source);

    for (unsigned int qp = 0; qp <  qrule.n_points(); qp++)
      PowerEmitted  += JxW[qp] * heat_source[qp] /s_0/s_0/s_0;


  } // end loop over elements


  return  PowerEmitted;
}


  //Equilibrium Energy

//   SimulationInterface* macro_sim;
//   macro_sim = SimulationInterface::find_simulation(myopts.macro_sim);

//   std::set<ID> T_var;
//   ID TT = 0;
//   if (macro_sim != NULL)
//   {
//     std::vector< std::map< ID, double > > temp_sol; TT =  macro_sim->get_variable_id("temperature");
//     T_var.insert(TT);
//   }

//   DofMap& dof_map = my_system->get_dof_map();
//   //equilibrium_energy->add(0.0);
//   std::vector<unsigned int> dof_indices;
//     MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
//     const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();


//        for ( ; it != end; ++it)
//        {
//          const Elem* elem = *it;

//          dof_map.dof_indices (elem, dof_indices);

//          std::vector<Point> pp(1);
//          pp[0]=elem->centroid();


//          if (macro_sim != NULL)
//          {
//            //Per elementi
//            std::vector< std::map< ID, double > > temp_sol;
//            macro_sim->get_solution(elem,pp,T_var,temp_sol);
//            //eq_energy[elem] = ((temp_sol[0].find(TT)-> second)-myopts.ref_temp)/(4 * M_PI);


//            //Per nodi
//          macro_sim->get_solution(elem,T_var,temp_sol);

//          for (unsigned n = 0 ; n< elem->n_nodes(); ++n)
//          {
//            //double energy = ((temp_sol[n].find(TT)-> second)-myopts.ref_temp)/(4 * M_PI);
// 	   double energy = (temp_sol[n].find(TT)-> second);
//            equilibrium_energy->set(dof_indices[n],energy);
//          }

//         }
//            else
//            {
// 	     // eq_energy[elem] = (myopts.eq_temp - myopts.ref_temp)/(4 * M_PI);
//              for (unsigned n = 0 ; n< elem->n_nodes(); ++n)
//              {
//                //double energy = (myopts.eq_temp - myopts.ref_temp)/(4 * M_PI);
//                double energy = myopts.eq_temp;
	      
//                equilibrium_energy->set(dof_indices[n],energy);
//              }
//            }
//        }


//   //FIll the guess direction energy---------------
//   {
//     sol_dir.resize(AngInt.n_slices,(my_system->solution)->clone().release());
//     DofMap& dof_map = my_system->get_dof_map();
//     //equilibrium_energy->add(0.0);
//     std::vector<unsigned int> dof_indices;
//     MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
//     const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
    
//     for ( ; it != end; ++it)
//     {
//       const Elem* elem = *it;
//       dof_map.dof_indices (elem, dof_indices);
//       for (ID n = 0 ; n< elem->n_nodes(); ++n)
//       {
// 	double eq_energy = (*equilibrium_energy)(dof_indices[n]);
	
// 	for (unsigned int k = 0; k<AngInt.n_slices; k++ )
// 	{
// 	  //double dir_energy = eq_energy * AngInt.d_omega[k] / (4.0 * M_PI);
// 	  double dir_energy = eq_energy;
// 	  //	  sol_dir[k]->set(dof_indices[n],dir_energy);
// 	  sol_dir[k] = (equilibrium_energy)->clone().release();
// 	}
//       }
//     }
//   }
  //file<<std::setw(10)<<iter<<"     "<<std::setw(10)<<err<<"\n";

  //if  (SimulationOptions::verbose() > 1)
    // {
    // std::cout<<"Iter:   "<< iter <<"  Error Jq  :"<<err<<std::endl;
    //std::cout<<"Iter:   "<< iter <<"  Error Eq  :"<<err2<<std::endl;
    // }
// void MicroHeatBalance::compute_total_flux(void)
// {


//    thermal_flux_x->zero();
//    thermal_flux_y->zero();
//    thermal_flux_z->zero();



//     ID var = my_system->variable_number("T");
//     const unsigned int system_number = my_system->number();

//     MeshBase::const_node_iterator  nd  = mesh->active_nodes_begin();
//     const MeshBase::const_node_iterator nd_end = mesh->active_nodes_end();


//     for ( ;  nd != nd_end ; ++nd)
//     {
//       Node* node = *nd;
//       const unsigned int  n_dof = node->dof_number(system_number,var,0);

//       for (unsigned int k = 0; k<AngInt.n_slices; k++ )
//       {

//       d_omega = AngInt.d_omega[k];
//       IntDir = AngInt.directions[k];

//      double sol = (*sol_dir[k])(n_dof);

//      //      double Jx =  IntDir(0) * sol;
//      // thermal_flux_x->add(n_dof,Jx);

//      double Jx =  vg * dir(0) * sol;
//      thermal_flux_x->add(n_dof,Jx);

//       double Jy =  IntDir(1) * sol;
//       thermal_flux_y->add(n_dof,Jy);

//       double Jz =  IntDir(2) * sol;
//       thermal_flux_z->add(n_dof,Jz);

//       }




//     }


// }


// void MicroHeatBalance::compute_total_equilibrium_energy(void)
// {

//   equilibrium_energy->zero();

//   ID var = my_system->variable_number("T");
//    const unsigned int system_number = my_system->number();

//    MeshBase::const_node_iterator  nd  = mesh->active_nodes_begin();
//    const MeshBase::const_node_iterator nd_end = mesh->active_nodes_end();


//    for ( ;  nd != nd_end ; ++nd)
//    {
//      Node* node = *nd;
//      const unsigned int  n_dof = node->dof_number(system_number,var,0);

//      for (unsigned int k = 0; k<AngInt.n_slices; k++ )
//      {

//        d_omega = AngInt.d_omega[k];
//        double sol = (*sol_dir[k])(n_dof);
//        double eq = d_omega * sol/(4.0 * M_PI);
//        //equilibrium_energy->add(n_dof,eq);

//        double discr_sol = (*sol_dir[k])(n_dof)* 0.5;

//        equilibrium_energy->add(n_dof,discr_sol);

//      }
//    }

//    equilibrium_energy->close();

// }



// void MicroHeatBalance::compute_flux(void)
// {


//       ID var = my_system->variable_number("T");
//       const unsigned int system_number = my_system->number();

//       MeshBase::const_node_iterator  nd  = mesh->active_nodes_begin();
//       const MeshBase::const_node_iterator nd_end = mesh->active_nodes_end();


//         for ( ;  nd != nd_end ; ++nd)
//         {
//           Node* node = *nd;
//           const unsigned int  n_dof = node->dof_number(system_number,var,0);
//           double sol = (*(my_system->solution))(n_dof);

//           //double Jx = d_omega * IntDir(0) * sol;
//           //thermal_flux_x->add(n_dof,Jx);

// 	  double Jx = vg * dir(0) * sol;
//           thermal_flux_x->add(n_dof,Jx);

//           double Jy = d_omega * IntDir(1) * sol;
//           thermal_flux_y->add(n_dof,Jy);

//           double Jz = d_omega * IntDir(2) * sol;
//           thermal_flux_z->add(n_dof,Jz);

//         }




// }






// //--------------------------------------------------------------------------------//
// void MicroHeatBalance::compute_equilibrium_energy(void)
// {

//   DofMap& dof_map = my_system->get_dof_map();
//    //equilibrium_energy->add(0.0);
//    std::vector<unsigned int> dof_indices;
//     MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
//      const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();

//         std::vector< std::map< ID, double > > temp_sol;
//         for ( ; it != end; ++it)
//         {
//           const Elem* elem = *it;

//           dof_map.dof_indices (elem, dof_indices);

//           std::vector<Point> pp(1);
//           pp[0]=elem->centroid();

//           get_solution_secure(elem,pp,e0_var,temp_sol);

//           eq_energy[elem] = temp_sol[0].find(E0)-> second;

//          // eq_energy[elem] = 0.0;

//         }


// }
