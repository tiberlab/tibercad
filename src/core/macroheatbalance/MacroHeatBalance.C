#include "MacroHeatBalance.h"
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
#include "HeatModel.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Boundary.h"
#include "Reservoir.h"
#include "FluxContact.h"
#include "SimulationOptions.h"

using namespace std;
MacroHeatBalance* MacroHeatBalance::static_this;
Device* MacroHeatBalance::_device;
//-----------------------------------------------------------------//


void MacroHeatBalance::parse_options( )
{ 

  const ModelOptions& sim_opt = get_options();

  opt.kappa_solve = sim_opt.get_option("kappa_temperature_sc", "false");
 
  if (opt.kappa_solve.compare("true") == 0 )
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

  // we calculate in cm!
  double mesh_units = get_scaling().get_calc_mesh_units() / sim_opt.get_option("Work_length_units", 1e-2);
  get_scaling().set_calc_mesh_units(mesh_units);
  
  equation_systems = & (get_equation_systems());

  system_name = get_equation_system_name();

  equation_systems->add_system<LinearImplicitSystem> (system_name);

  my_system = &( equation_systems->get_system<LinearImplicitSystem>(system_name)  );

  my_system->add_variable("T", FIRST);

  // my_system->add_variable("J", FIRST);

   // Insert the pointer to function that LibMesh library has to use    
  my_system->attach_assemble_function (assemble_heat_matrix);

   // Initialize the data structures for the equation system.
  my_system->init();	


  //Inizialize the solution to temperature of simulation options


  my_system->solution->zero();
   
  my_system->solution->add(SimulationOptions::temperature);
     
    
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
  //equation_systems->delete_system(system_name);

}
//---------------------------------------------------------------------------------//
MacroHeatBalance::MacroHeatBalance()
{
  

}
//----------------------------------------------------------------------------------//
PhysicalModel*
MacroHeatBalance::create_physical_model(const ModelOptions &options,
    const Material* mat) const throw (ModelErrorException)
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
 
   const string& modelname = options.get_option("type", "heat_reservoir");
 
  
   ThermalContact* model = ThermalContact::create(modelname, options);

   if (model == NULL)  
     throw ModelErrorException("MacroHeatBalance: No such boundary model: " + modelname);

  return model;

}
//----------------------------------------------------------------------------------//
MacroHeatBalance*  MacroHeatBalance::create (void)
{
  return new MacroHeatBalance;
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

  //AutoPtr<FEBase>  fe_face(FEBase::build(dim, fe_type));


  AutoPtr<FEBase>  fe_face(build_finite_element(dim,fe_type,true));
  
  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // o cout<<"Start loop over lattice thermal conductivity"<<endl;f the element.map
  QGauss qface(dim-1, THIRD);
  
  // Tell the finite element object to use our
  // quadrature rule.
  
  fe_face->attach_quadrature_rule(&qface);
  

  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
  
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  
  const std::vector<Point>& qface_point = fe_face->get_xyz();
  
  const std::vector<Point>& normal = fe_face->get_normals();
  
 
  std::vector<unsigned int> dof_indices;
  
  DenseMatrix<Number>  Ke;

  DenseVector<Number>  Fe;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  
  //Model Variables

  Tensor2Sym kappa; 

  double eTEpower;

  double hTEpower;

  std::vector<double> heat_source;

  std::vector<RealGradient> flux_heat_source;
 
 
 
  ThermalContact* contact; 
  //----------------------------------------------------------LatticeThermalConductivity-------//

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
 
    init_heat_model(elem);
 
    heat_model->get_total_heat_source(q_point,heat_source);

    heat_model->get_total_flux_heat_source(q_point,flux_heat_source);

    heat_model->get_thermal_conductivity(kappa);
    
    for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
    { // loop over test function
      
      bool reservoir_contact = false;
      
      //!let us check if it belongs to a reservoir boundary
      const Node* nd = elem->get_node(p1);
      
      Boundary* bd = se.get_boundary(nd); 
      
      if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) )
      { //if belongs to boundary
	
	contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );
	
	if (contact->get_type() == ThermalContact::Reservoir)
	{//heat reservoir---{//loop over basis functions
	  reservoir_contact = true;
	}
      }
       
      if (reservoir_contact)
      {//heat reservoir---
	
	Ke(p1,p1) = 1.0;
	
	Fe(p1) = ( dynamic_cast<Reservoir*> (contact) )->get_temperature();
	
      } 
      else 
      {
	for (unsigned int qp=0; qp<qrule.n_points(); qp++)  
	{//Loop over quadrature points 
	  
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
		
		
		value += JxW[qp] * ( dphi[p1][qp](i) * kappa_value - flux_heat_source[qp](j) ) * dphi[p2][qp](j);
		
	      }//end loop over direction (2)
	      
	    }//end loop over direction (1)											         
	    
	    Ke(p1,p2) += value;
	    
	  } //loop over basis functions
	  
	  Fe(p1) +=JxW[qp] * heat_source[qp] * phi[p1][qp];
	  
	}//end Loop over quadrature points  
	
	
      } //end if it belongs to boundary
      
      
    } // end loop over test functions
    



  //    //The loop over element is the only loop that is surviving at this point

  

   //   //Contact Resistance
            
//      for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
//       {
//        //!let us check if it belongs to a reservoir boundary
//        const Node* nd = elem->get_node(p1);
       
//        Boundary* bd = se.get_boundary(nd); 
       
//        if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) )
//        {
// 	 ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );    
	 
//  	 if (contact->get_type() == ThermalContact::Neumann)  
// 	 {

// 	   double rho_e = ( dynamic_cast<FluxContact*> (contact) )->get_electrons_resistivity();
	
//            double rho_h = ( dynamic_cast<FluxContact*> (contact) )->get_holes_resistivity();
            
	    

//  	   const unsigned int num_sides = elem->n_sides();
	 
//  	   for (unsigned int side = 0; side<num_sides; side++)
//  	   {
//  	     const ElementSide elside(elem->top_parent(), side);
	     
//  	     if ( (heat_model->get_dd_environment()).is_on_boundary(   elside   ) ) //if belongs to a boundary of current(!) simulation
//  	     {
	     
//  	       fe_face->reinit(elem, side);
	       
// 	       heat_model->get_dd_solution_secure(qface_point,QfermiE,QfermiH,JE,JH);
	       
// 	       for (short i = 0; i < 3; i++)
// 	       {
// 		 for (unsigned int qp=0; qp < qface.n_points(); qp++)
// 		 {
		   
// 		   Fe(p1) -= Constants::e * JxW_face[qp] * phi_face[p1][qp] * (JE[qp](i)*JE[qp](i) * rho_e + JH[qp](i) * JH[qp](i) * rho_h); 
		     
// 		 }
// 	       }
	       
// 	     } // if on boudary of dd
	     
// 	   } // side
      
//        }//   if (contact->get_type() == ThermalContact::Neumann) 
//       }// if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) )
     
//    }// for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
  
  
  dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
  system.matrix->add_matrix (Ke, dof_indices);
  system.rhs->add_vector    (Fe, dof_indices); 
  
  // system.matrix->print_matlab("matr.m"); 
  
  
  } //End Loop over elements
  
  
  //   system.matrix->print_matlab("Matr.m");
  //   system.rhs->print();
    
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




ID
MacroHeatBalance::convert_variable_name_to_id(const string& variable_name) const
{

  ID id = INVALID_ID;


    if (variable_name == "temperature" )
       id  = TEMPERATURE;
    

  return id;
}



void
MacroHeatBalance::get_solution_secure(const Elem* elem,
    const set<ID>& ids, vector<map<ID, double> >& values)
{

  LinearImplicitSystem& system = *my_system;

  DofMap& dof_map =  system.get_dof_map();

  vector<unsigned int> dof_indices;

  dof_map.dof_indices(elem, dof_indices);  

  for (unsigned int n = 0; n < elem->n_nodes(); n++)
  {
    
     if (ids.count(TEMPERATURE))
      values[n][TEMPERATURE] = (*(system.solution))(dof_indices[n]);

  }
   

}



void
MacroHeatBalance::get_solution_secure(const Elem* elem, const vector<Point>& p,
    const set<ID>& ids, vector<map<ID, double> >& values)
{

  unsigned int np = p.size();
  values.resize(np);
  if ((np == 0) || (ids.size() == 0)) return;

  LinearImplicitSystem& system = *my_system;

  DofMap& dof_map =  system.get_dof_map();

  const unsigned int uvar = system.variable_number("T");

  FEType fe_type = dof_map.variable_type(uvar);
  
  //  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true)); 

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  vector<Point> points(np);

  //std::cout<<p.size()<<std::endl;
  FEInterface::inverse_map(dim, fe_type, elem, p, points);
 
  fe->reinit(elem, &points);

  vector<unsigned int> dof_indices;

  dof_map.dof_indices (elem, dof_indices);  

  const unsigned int n_dofs   = dof_indices.size();

  for (unsigned int n = 0; n < np; n++)
  {
     double T = 0;

    //do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
   
      T  += phi[i][n] * (*(system.solution))(dof_indices[i]);
    }
    
     if (ids.count(TEMPERATURE))
      values[n][TEMPERATURE] = T;
  }
   
}










void
MacroHeatBalance::build_elemental_results(const std::set<std::string>& variables,
					  std::vector<double>& results, std::vector<std::string>& legend)
{

  
  // we only do something if we are on processor 0
  // TODO parallelize
  if (libMesh::processor_id() != 0)
    return;

  const set<string>::const_iterator varend(variables.end());
  
  vector<ID> ids;
  unsigned int nm; 
  


  legend.resize(0); 
  unsigned int n_vars = 0;  
  
  int HS = -1;
  if (variables.find("HeatSource") != varend)
  {
    const Device& device = *(_device);

    HS = n_vars;

    std::map<ID, std::map<ID,std::string> > heat_source_ids;
    
    MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();

    MeshBase::const_element_iterator it_end =    mesh->active_local_elements_end(); 

    // assert(it_end != mesh->active_local_elements_end());
    
    HeatModel* heat_model = NULL;
    
    const Elem* elem = *it;
    
    ID subdomain = elem->subdomain_id();
    
    heat_model= dynamic_cast<HeatModel*>(
					 device.get_material(subdomain)->get_model(get_id()));
    
    nm = heat_model->get_heat_source_IDs(ids);
    
    for (int i = 0; i < nm; i++)
    {
         
      std::vector<std::string> local_legend  =
 	heat_model->get_heat_source_model(ids[i])->get_legend(); 
      
      unsigned int ns_loc = local_legend.size();

      unsigned int ns = 0;
      for (unsigned int n_loc = 0;  n_loc < ns_loc; ++n_loc )
      {
          
	legend.push_back(local_legend[n_loc]);
        
	n_vars++;
	ns++;
	 
      }
      if (ns>1)
      {
	n_vars++;
	legend.push_back( heat_model->get_heat_source_model(ids[i])->get_default_name() ); 
      }
    }
    if (nm>1)
    {
      n_vars++;
      legend.push_back("TotalHeatSource"); 
    }
    
  }

  const unsigned int nn  = mesh->n_nodes();
    
  //results.resize(nn * n_vars,0.0);
  results.resize(0,0.0);

  legend.resize(n_vars);
 
  std::vector<unsigned int> dof_indices;
  DofMap& dof_map = my_system->get_dof_map();
  
  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
  
 
  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  { 
    

    const Elem* elem = *it;

    unsigned int id = n_vars * elem_number;
    
    dof_map.dof_indices (elem, dof_indices); 

    std::vector<Point> _node(1);

    _node[0]=(elem->centroid());
      
    if (HS != -1)
    {
           
      ID subdomain = elem->subdomain_id();
	
      const Material* mat = _device->get_material(subdomain);
	
      HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
	
      heat_model->set_element(elem);     
	
      heat_model->re_init();
	
      unsigned int k = 0;
      
      std::vector<double> total_heat_source(nm);

      total_heat_source.clear();

      total_heat_source.resize(nm);
	
      for (int i = 0; i < nm; i++)
      {
	  
	std::vector<std::vector<double> > heat_source;
	 
	heat_model->get_heat_source_model(ids[i])->get_heat_source_output(_node,elem,heat_source);
	 
	unsigned int ns = 0;
      
	for (unsigned int n_loc = 0 ; n_loc<heat_source[0].size() ; ++n_loc)
	{ 
	   
	  //results[id + HS + k] = heat_source[0][n_loc];
          results.push_back(heat_source[0][n_loc]);
	 
	  total_heat_source[i] =  total_heat_source[i] + heat_source[0][n_loc];

	  //std::cout<<id + HS + k<<std::endl;
          //std::cout<< results.size() <<std::endl;
       
	  ++k;
	  ++ns;
	  
	}
	
	//Write the total heat of the i-heat source
	if (ns>1)
	{
	  results.push_back( total_heat_source[i]);
	  //  results[id + HS + k] = total_heat_source[i];
          // std::cout<<id + HS + k<<std::endl;
	}
	
      }
      
      //Write the total heat
      if (nm>1)
      {
	double total_heat = 0.0;
	
	for (int i = 0; i < nm; i++)
	{
	  total_heat = total_heat + total_heat_source[i];
	  
	}
	results[id + HS + k] = total_heat;
	
	
      }//	if (nm>1)
      
    } //if (HS != -1)

      
        elem_number++;
  } //over elements

  results.resize(elem_number * n_vars);
}





//----------------------------------------------------------------------------------//
void MacroHeatBalance::build_nodal_results (const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				     std::vector< std::string > &legend)
{

  

  
 
  const set<string>::const_iterator varend(variables.end());
  

  legend.resize(0); 
  unsigned int n_vars = 0;

  int Temp = -1;
  if (variables.find("T") != varend)
  {
    Temp = n_vars; 
    legend.push_back("Temperature[K]");
    n_vars++;
  }

  const unsigned int nn  = mesh->n_nodes();
  
  results.resize(nn * n_vars,0.0);
  
  std::vector<unsigned int> dof_indices;
  DofMap& dof_map = my_system->get_dof_map();
  
  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
  
 
  for ( ; it != end; ++it)
  { 
    const Elem* elem = *it;
    
    dof_map.dof_indices (elem, dof_indices); 
     
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    {
      
      unsigned int id =  (elem->node(n) * n_vars) ;
      
      if (Temp != -1)
      {
	
	results[id+Temp]  =  (*(my_system->solution))(dof_indices[n]);       
	
      }
      
    }  
     
  }

}
