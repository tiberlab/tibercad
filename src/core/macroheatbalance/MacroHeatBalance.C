#include "MacroHeatBalance.h"
#include "DriftDiffusion.h"
#include "BoundaryProperties.h"
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

#include "LatticeThermalConductivity.h"
#include "ZbLatticeThermalConductivity.h"
#include "ThermoelectricPower.h"
#include "HeatModel.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Boundary.h"
#include "Reservoir.h"
#include "SimulationOptions.h"

using namespace std;
MacroHeatBalance* MacroHeatBalance::static_this;
Device* MacroHeatBalance::_device;
//-----------------------------------------------------------------//


void MacroHeatBalance::parse_options( )
{ 

  const ModelOptions& sim_opt = get_options();

  opt.kappa_solve = sim_opt.get_option("kappa_solve", "no_self");
 
  if (opt.kappa_solve.compare("self_consistent") == 0 )
  {  
    opt.max_error = sim_opt.get_option("max_error",1e-2);
  }

   
}

void MacroHeatBalance::do_init( ) 
{
  

  const ModelOptions& sim_opt = get_options();

  SimulationEnvironment& si = get_environment();   

  _device = &( si.get_device() );

  mesh = & (_device->get_mesh());

  dim = mesh->mesh_dimension();


  double mesh_units = _device->get_mesh_units();

  opt.work_units = sim_opt.get_option("Work_length_units", 1e-2);

  opt.length_scale = mesh_units/opt.work_units;



  equation_systems = & (get_equation_systems());

  system_name = get_equation_system_name();

  equation_systems->add_system<LinearImplicitSystem> (system_name);

  my_system = &( equation_systems->get_system<LinearImplicitSystem>(system_name)  );

  my_system->add_variable("T", FIRST);


   // Insert the pointer to function that LibMesh library has to use    
  my_system->attach_assemble_function (assemble_heat_matrix);

   // Initialize the data structures for the equation system.
  my_system->init();	




  //Inizialize the solution to temperature of simulation options

    MeshBase::const_node_iterator       nd     = mesh->active_nodes_begin();
    const MeshBase::const_node_iterator nd_el  = mesh->active_nodes_end();

    unsigned int number_of_points = 0;
    for ( ; nd != nd_el ; ++nd)  number_of_points++;

    for ( unsigned int n = 0 ; n != number_of_points; ++n)
    {
      (*(my_system->solution)).set(n,SimulationOptions::temperature); 
    }
   //-----------------------------------------------------------------

   
  //------init is done---------------------------------------------------------------------//

}
//-------------------------------------------------------------------------------//
void  MacroHeatBalance::do_solve()
{
  
  //  std::cout<<"Solving macro heat balance: start"<<std::endl;
  
  parse_options();
  
  static_this = this;
  
  my_system->solve();
  
  
  if (opt.kappa_solve.compare("self_consistent") == 0)
  {
    cout<<endl;
    cout<<"Start loop over lattice thermal conductivity"<<endl; 
    cout<<endl;
    
    double norm_error;
    
    norm_error = opt.max_error + 1;
    
    //Inizialize the old_soluction----------------------------------------
    vector<double> old_solution((*(my_system->solution)).size());
    

    for ( unsigned int n = 0 ; n != (*(my_system->solution)).size(); ++n) 
    { old_solution[n] = (*(my_system->solution))(n);}
    
     //-------------------------------------------------------------------  
    
    for (; norm_error> opt.max_error;)
    {
      
      
      my_system->solve();
      
       
      //Compute error and old_solution
      norm_error = 0.0;
      for ( unsigned int n = 0 ; n !=  (*(my_system->solution)).size(); ++n)
      {
	norm_error += (old_solution[n]-(*(my_system->solution))(n)) * (old_solution[n]-(*(my_system->solution))(n));
	old_solution[n] = (*(my_system->solution))(n);
        
      }
      
      
      norm_error = sqrt(norm_error);
      
      //------------------------ 
      
       cout<<"Error_norm = " <<norm_error<<endl;
       
     } //end for 
    cout<<endl;
    cout<<"End loop over lattice thermal conductivity"<<endl;	
    cout<<endl;   
    
  }//end if
  
  
 
  if (get_options().get_option("update_temperature",true))
  {


    
    std::vector<unsigned int> dof_indices;
    
    DofMap& dof_map = my_system->get_dof_map();
    
    MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh->active_elements_end();
    
    double  T=0.0;
    
    for ( ; el != end_el ; ++el)   //loop over elements
    {  
      const Elem* elem = *el;
      
      dof_map.dof_indices (elem, dof_indices); 
      
      T=0.0;
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {   
	T +=  (*(my_system->solution))(dof_indices[n]); 
      }
      
      T /= elem->n_nodes();
      _device->set_temperature(elem,T);
      
    }
    ////////////////////////////////////
    
  }
  
  
  
  //  std::cout<<"Solving macro heat balance: finish"<<std::endl;

}


//--------------------------------------------------------------------------------//
MacroHeatBalance::~MacroHeatBalance()
{
  equation_systems->delete_system(system_name);

}
//---------------------------------------------------------------------------------//
MacroHeatBalance::MacroHeatBalance()
{
  

}
//----------------------------------------------------------------------------------//
PhysicalModel*   MacroHeatBalance::create_physical_model (const ModelOptions &options) const 
                    throw (ModelErrorException)
{
  
  HeatModel* model = dynamic_cast<HeatModel*> ( PhysicalModelInterface::create("thermal",options) );

  if (model == NULL) 
    throw ModelErrorException("MacroHeatBalance: Thermal physical model is not created" );

  return model;      

}
//----------------------------------------------------------------------------------//

BoundaryProperties* MacroHeatBalance::create_boundary_model (const ModelOptions &options) const 
                    throw (ModelErrorException)

{

   const string& modelname = options.get_option("type", "Heat_reservoir");
 
   ThermalContact* model = ThermalContact::create(modelname, options);

   if (model == NULL)  
     throw ModelErrorException("MacroHeatBalance: No such boundary model: " + modelname);

  return model;

}
//----------------------------------------------------------------------------------//
MacroHeatBalance*  MacroHeatBalance::create (void)
{
  return new MacroHeatBalance();
}

//----------------------------------------------------------------------------------//
void MacroHeatBalance::build_nodal_results (const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				     std::vector< std::string > &legend)
{

 

  if (variables.find("T") != variables.end())
  {
    legend.resize(1);
    legend[0] = "Temperature[K]";



    MeshBase::const_node_iterator       nd     = mesh->active_nodes_begin();
    const MeshBase::const_node_iterator nd_el  = mesh->active_nodes_end();

    unsigned int number_of_points = 0;
    for ( ; nd != nd_el ; ++nd)  number_of_points++;

    results.resize(number_of_points, 0.0);

    
    MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();

    std::vector<unsigned int> dof_indices;

    DofMap& dof_map = my_system->get_dof_map();

    for ( ; it != end; ++it)
    { 
      const Elem* elem = *it;

      dof_map.dof_indices (elem, dof_indices); 

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
	unsigned int id =  elem->node(n);
        
	results[id]  =  (*(my_system->solution))(dof_indices[n]);       
	//   std::cout<<results[id]<<std::endl;

      }
    }
  }
}





double MacroHeatBalance::get_temperature_element(const Elem* elem) const
{
        

       std::vector<unsigned int> dof_indices;

       DofMap& dof_map = my_system->get_dof_map();

       dof_map.dof_indices (elem, dof_indices); 

       double T = 0.0;
       for (unsigned int n = 0; n < elem->n_nodes(); n++)
       {   
	  T +=  (*(my_system->solution))(dof_indices[n]); 
       }
      
       T /= elem->n_nodes();

       return T;

}



std::vector<double>  MacroHeatBalance::get_temperature_node(const Elem* elem)
{

      std::vector<double> Tloc(elem->n_nodes());

       std::vector<unsigned int> dof_indices;

       DofMap& dof_map = my_system->get_dof_map();

       dof_map.dof_indices (elem, dof_indices); 

       for (unsigned int n = 0; n < elem->n_nodes(); n++)
       {   
	  Tloc[n] =  (*(my_system->solution))(dof_indices[n]); 
       }
    
       return Tloc;

}



//----------------------------------------------------------------------------------//
void MacroHeatBalance::assemble_heat_matrix(EquationSystems& es,
				     const std::string& system_name)
{



   static_this->do_assemble( es, system_name);

}

//----------------------------------------------------------------------------------//
void MacroHeatBalance::do_assemble(EquationSystems& es, const std::string& system_name)
{

  SimulationEnvironment& se = get_environment(); 

  LinearImplicitSystem& system = *my_system;
  
  const unsigned int uvar = system.variable_number("T");
  
  DofMap& dof_map = system.get_dof_map();
  
  FEType fe_type = dof_map.variable_type(uvar);
  
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  //AutoPtr<FEBase> fe (build_finite_element(dim, fe_type));
  
  QGauss qrule (dim, FIFTH); //may be could be decreased (CHECK!!!)
  
  // quadrature rule  
  fe -> attach_quadrature_rule (&qrule);
  
  
  // Declare a special finite element object for
  // boundary integration.

  AutoPtr<FEBase>  fe_face(FEBase::build(dim, fe_type));






  // AutoPtr<FEBase>  fe_face(build_finite_element(dim, fe_type));
  
  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // o cout<<"Start loop over lattice thermal conductivity"<<endl;f the element.map
  QGauss qface(dim-1, THIRD);
  
  // Tell the finite element object to use our
  // quadrature rule.
  
  fe_face -> attach_quadrature_rule (&qface);
  


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

 
  std::vector<unsigned int> dof_indices;
  
  DenseMatrix<Number>  Ke;

  DenseVector<Number>  Fe;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  
  //Model Variables

  Tensor2Sym kappa; 

  double eTEpower;

  double hTEpower;

  //------------------------------
  

  //Define the vectors which store the Drift Diffusion solutions

   

  std::vector<DriftDiffusion::Solution>  potentials;

  std::vector<DriftDiffusion::Currents>   currents;


  //----------------------------------------------------------LatticeThermalConductivity-------//
  

  //-----------------------------------------------------------------//
  //My Jacobian. It is to pass to our work units
  
  double my_Jacobian = 1.0;
  for (short i = 1; i <= dim; i++)  my_Jacobian *= opt.length_scale;
  //----------------------------------------------------------------//


 
   

  for ( ; el != end_el ; ++el)   //loop over elements
  {  

  
    //Inizialize the element environment
    const Elem* elem = *el;

    dof_map.dof_indices (elem, dof_indices);  

    const unsigned int n_dofs   = dof_indices.size();
    
    fe->reinit (elem);
    
    Ke.resize(n_dofs, n_dofs);

    Fe.resize(n_dofs);
    
    Fe.zero();  
     
    Ke.zero();

    //Init e read from heat model

    init_heat_model(elem);
  
    heat_model->get_thermal_conductivity(kappa);

    eTEpower = heat_model->get_electrons_thermoelectric_power();	
   
    hTEpower = heat_model->get_holes_thermoelectric_power();

    if (heat_model->get_joule_opt())
    {
      heat_model->get_dd_solution(q_point,potentials,currents); 
      
    } 
    
  
  
 
    
      
     for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
      
    { // loop over test function

      //!let us check if it belongs to a reservoir boundary
      const Node* nd = elem->get_node(p1);
	
      Boundary* bd = se.get_boundary(nd); 
	
      if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) )

      { //if belongs to boundary

	ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );
	if (contact->get_type() == ThermalContact::Reservoir)
	{//heat reservoir---{//loop over basis functions
	    
	  Ke(p1,p1) = 1.0;
	  
	  Fe(p1) = ( dynamic_cast<Reservoir*> (contact) )->get_temperature();
        
	  
	}//end reservoir
      } 
      else 
      {
	for (unsigned int qp=0; qp<qrule.n_points(); qp++)  
	{//Loop over quadrature points 

	  //------------------------
	  //volume integration for matrix
	  // - \int_V \frac{\patial \phi_{\alpha}{\partial x_i}\kappa_{i,j} \frac{\par tial \phi_{\beta}{\partial x_j} dx
	  	  
	   for (unsigned int p2=0; p2<n_dofs; p2++) 

	  {//loop over basis functions
            

	    double value = 0.0;
	    
	    
	    for (short i = 0; i < dim; i++) 
	    {//loop over direction (1); test function derivative

	      for (short j = 0; j < dim; j++)
	      {//loop over direction (2); basis function derivative

		double kappa_value;
		if (i < j) 
		  kappa_value = kappa(j+1, i+1);
		else
		  kappa_value = kappa(i+1, j+1);

	   
	          value += -JxW[qp] * kappa_value * dphi[p1][qp](i) * dphi[p2][qp](j) /(opt.length_scale * opt.length_scale);
                 
		
	      }//end loop over direction (2)

	      //This includes the Peltier-Thomson effects   
             
	      if (heat_model->get_peltier_thomson_opt())
	      {
		//cout<<"Peltier"<<endl;
		
		value += -JxW[qp]* phi[p1][qp] * dphi[p2][qp](i)*
		  (eTEpower*currents[qp].jn(i) + hTEpower*currents[qp].jp(i) ) /(opt.length_scale);

		//cout<<(opt.length_scale)<<endl;
	      }
	      
	       
	    }//end loop over direction (1)												
	    
	   value *= my_Jacobian;
	    

           Ke(p1,p2) += value;
	    
	  } //loop over basis functions
	  
	 
	   if (heat_model->get_joule_opt())
	   {
	     //Inclusion of jolue effect  
	     
	     for (short i = 0; i < dim; i++) 
	     { Fe(p1) -= dphi[p1][qp](i) * JxW[qp] *
		 ( currents[qp].jn(i)*potentials[qp].fermi_e + potentials[qp].fermi_h * currents[qp].jp(i) )
		 / opt.length_scale * my_Jacobian;
	     }
	     
	   }   
	   
	   
	   
	}//end Loop over quadrature points  
	

      } //end if it belongs to boundary
 

    } // end loop over test functions
      
     //The loop over element is the only loop that is surviving at this point
  
     
     if (heat_model->get_joule_opt())
     {
       for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
       {
	 //!let us check if it belongs to a reservoir boundary
	 const Node* nd = elem->get_node(p1);
	 
	 Boundary* bd = se.get_boundary(nd); 
	 
	 bool belongs_to_reservoir = false;
	 
	 const unsigned int num_sides = elem->n_sides();
	 
	 if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) )
	 {
	   ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );    
	   
	   if (contact->get_type() == ThermalContact::Reservoir)   belongs_to_reservoir = true;
	 }
	 
	 
	 
	 if ( !belongs_to_reservoir ) 
	 {//not fixed temperature2
	   
	   const unsigned int num_sides = elem->n_sides();
	   
	   for (unsigned int side = 0; side<num_sides; side++)
	   {
	     const ElementSide elside(elem->top_parent(), side);
	     
	     if ( (heat_model->get_dd_environment()).is_on_boundary(   elside   ) ) //if belongs to a boundary of current(!) simulation
	     {
	      
	       std::vector<DriftDiffusion::Solution>  face_potentials;   
	       
	       std::vector<DriftDiffusion::Currents>   face_currents;
		   

	       if (dim > 1)
	       {
		 
		 fe_face->reinit(elem, side);
		 
		 const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
		 
		 const std::vector<Real>& JxW_face = fe_face->get_JxW();
		 
		 const std::vector<Point>& qface_point = fe_face->get_xyz();
		 
		 const std::vector<Point> & normal = fe_face->get_normals();
		 
		    
		 heat_model->get_dd_solution(qface_point,face_potentials,face_currents);  
		     
		 
		 for (short i = 0; i < 3; i++)
		 {
		   for (unsigned int qp=0; qp < qface.n_points(); qp++)
		   {
		     Fe(p1) += (JxW_face[qp] * phi_face[p1][qp]) * normal[qp](i) * 
		       ( face_currents[qp].jn(i)*face_potentials[qp].fermi_e + face_potentials[qp].fermi_h * face_currents[qp].jp(i) ) * 
		       (my_Jacobian/opt.length_scale);
		   }
                 }


	       }
	       else //dim = 1
	       {
		 
		 if (p1== side)
		 {
		   std::vector<double> normal(3);
		   Point p = elem->point(side);
		   Point pc = elem->centroid();
		   
		   
		   
		   const double temp = sqrt((p(0) - pc(0)) * (p(0) - pc(0)) 
					    +(p(1) - pc(1)) * (p(1) - pc(1)) + 
					    (p(2) - pc(2)) * (p(2) - pc(2)));
		   
		   
		   for (short i = 0; i < 3; i++)
		     normal[i] = (p(i) - pc(i))/temp;
		   
		   
		   std::vector<Point> qface_point(1);
		   
		   qface_point[0] = elem->point(p1);
		   
		   
	
		   heat_model->get_dd_solution(qface_point,face_potentials,face_currents);  		   
		   
		   for (short i = 0; i < 3; i++)
		   {  
		     Fe(p1) +=   normal[i] * 
		       ( face_currents[0].jn(i) * face_potentials[0].fermi_e + face_potentials[0].fermi_h * face_currents[0].jp(i) )   ;
		     
		     
		   }
		   
		   
		   
		   
		 } //if (dim > 1)
		 
	       }
	     } // if ( (dd_simul->get_environment()).is_on_boundary(   elside   ) ) //if belongs to a boundary of current(!) simulation
	   } //  for (unsigned int side = 0; side<num_sides; side++)
	   
	 } // if ( !belongs_to_reservoir ) 
	 
       }// for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
       
     } // (dd_simul != NULL)


     dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
     system.matrix->add_matrix (Ke, dof_indices);
     system.rhs->add_vector    (Fe, dof_indices); 
      
      
  } //End Loop over elements
    

 
  //   system.matrix->print_matlab("Matr.m");
  //  system.rhs->print();
    
} //do assembly

void  MacroHeatBalance:: init_heat_model(const Elem* elem)
{
       DofMap& dof_map = my_system->get_dof_map(); 
       
       std::vector<unsigned int> dof_indices;

       dof_map.dof_indices (elem, dof_indices);  

       ID subdomain = elem->subdomain_id();
    
       const Material* mat = _device->get_material(subdomain);
    
       heat_model =  (  dynamic_cast<HeatModel*> (  mat ->get_model(get_id()) )  );
    

       double Tloc = 0.0;
       for (unsigned int n = 0; n < elem->n_nodes(); n++)
       {   
	 Tloc +=  (*(my_system->solution))(dof_indices[n]); 
       }
      
       Tloc /= elem->n_nodes();


       //Configure model phase
   
       heat_model->set_element(elem);     

       heat_model->set_temperature(Tloc);

       //Update model for a given element
       heat_model->re_init();
     
}
