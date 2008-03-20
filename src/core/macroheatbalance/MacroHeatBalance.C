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
#include "SimulationEnvironment.h"
#include "ThermalResistance.h"

using namespace std;
MacroHeatBalance* MacroHeatBalance::static_this;
Device* MacroHeatBalance::_device;
//-----------------------------------------------------------------//


void MacroHeatBalance::parse_options( )
{ 
 
  //const ModelOptions& sim_opt = SimulationInterface::get_options();

  const ModelOptions& sim_opt = get_options();

  myopts.integration_order = static_cast<libMeshEnums::Order>(sim_opt.get_option("integration_order", 5));

  myopts.work_units = sim_opt.get_option("Work_length_units", 1e-2);

  myopts.kappa_solve = sim_opt.get_option("kappa_temperature_sc", "false");

  if (myopts.kappa_solve.compare("true") == 0 )
  {  
    myopts.max_error = sim_opt.get_option("max_error",1e-2);
  }

   
}

void MacroHeatBalance::do_init( ) 
{
  const ModelOptions& sim_opt = get_options();

  SimulationEnvironment& si = get_environment();   

  _device = &( si.get_device() );

  mesh = & (_device->get_mesh());

  dim = mesh->mesh_dimension();

  double mesh_units = get_scaling().get_calc_mesh_units() / sim_opt.get_option("Work_length_units", 1e-2);

  get_scaling().set_calc_mesh_units(mesh_units);
  
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


  my_system->solution->zero();
   
  my_system->solution->add(SimulationOptions::temperature);

  heat_legend = "Wq";

  JQ_var.insert(JQX);

  JQ_var.insert(JQY);

  JQ_var.insert(JQZ);


   //-----------------------------------------------------------------

   
  //------init is done---------------------------------------------------------------------//

}
//-------------------------------------------------------------------------------//
void  MacroHeatBalance::do_solve()
{
  
   
  parse_options();
  
  static_this = this;
  
  my_system->solve();
  

  
 //   if (opt.kappa_solve.compare("self_consistent") == 0)
//   {
//     cout<<endl;
//     cout<<"Start loop over lattice thermal conductivity"<<endl; 
//     cout<<endl;
    
//     double norm_error;
    
//     norm_error = opt.max_error + 1;
    
//     //Inizialize the old_soluction----------------------------------------
//     vector<double> old_solution((*(my_system->solution)).size());
    

//     for ( unsigned int n = 0 ; n != (*(my_system->solution)).size(); ++n) 
//     { old_solution[n] = (*(my_system->solution))(n);}
    
//      //-------------------------------------------------------------------  
    
//     for (; norm_error> opt.max_error;)
//     {
      
      
//       my_system->solve();
      
       
//       //Compute error and old_solution
//       norm_error = 0.0;
//       for ( unsigned int n = 0 ; n !=  (*(my_system->solution)).size(); ++n)
//       {
// 	norm_error += (old_solution[n]-(*(my_system->solution))(n)) * (old_solution[n]-(*(my_system->solution))(n));
// 	old_solution[n] = (*(my_system->solution))(n);
        
//       }
      
      
//       norm_error = sqrt(norm_error);
      
//       //------------------------ 
      
//        cout<<"Error_norm = " <<norm_error<<endl;
       
//      } //end for 
//     cout<<endl;
//     cout<<"End loop over lattice thermal conductivity"<<endl;	
//     cout<<endl;   
    
//   }//end if
  
  
 

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


  std::vector<double> heat_source;

  std::vector<RealGradient> flux_power;


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

    heat_model->get_total_power_flux(q_point,flux_power,false);

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

	Fe(p1) = (dynamic_cast<Reservoir*> (contact) )->get_temperature();
	
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
		
		value += JxW[qp] * dphi[p1][qp](i) * kappa_value * dphi[p2][qp](j); 
		
	      }//end loop over direction (2)
	      
	    }//end loop over direction (1)											         

            
	    Ke(p1,p2) += value;
	    
	  } //loop over basis functions
	  
	  Fe(p1) +=JxW[qp] * heat_source[qp] * phi[p1][qp];
             
	  //for (short i = 0; i < dim; i++)
	  //Fe(p1) += JxW[qp] * flux_power[qp](i) * dphi[p1][qp](i);
	  
	}//end Loop over quadrature points  
	
	
      } //end if it belongs to boundary
      
      
    } // end loop over test function


    //Thermal Resistance

 
    //The loop over element is the only loop that is surviving at this point
    
    const unsigned int num_sides = elem->n_sides();
    
    for (unsigned int side = 0; side<num_sides; side++)
    {
      const ElementSide elside(elem->top_parent(), side);
      
      Boundary* bd = se.get_boundary(elside);
      
      if  (bd != NULL)
      { 
	fe_face->reinit(elem,side);
	
        if (bd->get_boundary_properties( get_id() ) != NULL ) 
        {
	  
	  ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );    
	  
	  if (contact->get_type() == ThermalContact::ThermalResistance)   
	  { 
	   
	    
	    double rth =  ( dynamic_cast<ThermalResistance*> (contact) )->get_thermal_resistance();
	    
	    double Text =  ( dynamic_cast<ThermalResistance*> (contact) )->get_external_temperature();
	    
	    
	    for (unsigned int qp=0; qp < qface.n_points(); qp++)
	    {
	      
	      for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
	      {
		double Fe_plus = 0.0; 
		
		Fe_plus =  JxW_face[qp] * 1/rth * phi_face[p1][qp] * Text;
		
		Fe(p1) += Fe_plus;
		
		for (unsigned int p2=0; p2<n_dofs; p2++) //bases test
		{
		  double val_plus = 0.0; 
		  val_plus  =  JxW_face[qp] * 1/rth * phi_face[p1][qp] * phi_face[p2][qp];
		  
		  
		  Ke(p1,p2) += val_plus;
		  
		  
		}// (unsigned int p2=0; p2<n_dofs; p2++)
		
		
	      }//for (unsigned int p1=0; p1<n_dofs; p1++) 
	      
	      
	    }// for (unsigned int qp=0; qp < qface.n_points(); qp++)
	    
	  }//if (contact->get_type() == ThermalContact::ThermalResistance)  

        }// if (bd->get_boundary_properties( get_id() ) != NULL )


        //Other simulations
	//heat_model->get_total_power_flux(qface_point,flux_power,true);

	//for (unsigned int qp=0; qp < qface.n_points(); qp++)
	//{
	  
	// for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
	  // {

	   
	// double Fe_surf = JxW_face[qp] * phi_face[p1][qp] * flux_power[qp] * normal[qp];
	//   Fe(p1) -= Fe_surf;
	      

	//}
	//	}
	
      }// if (bd != NULL)
      
    }// for (unsigned int side = 0; side<num_sides; side++)
    



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
    if (variable_name == "JQx" )
       id  = JQX;
    if (variable_name == "JQy" )
       id  = JQY;
    if (variable_name == "JQz" )
       id  = JQZ;
    

  return id;
}



void
MacroHeatBalance::get_solution_secure(const Elem* elem,
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
MacroHeatBalance::get_solution_secure(const Elem* elem, const vector<Point>& p,
    const set<ID>& ids, vector<map<ID, double> >& values)
{

 
  unsigned int np = p.size();
  values.resize(np);
  if ((np == 0) || (ids.size() == 0)) return;

  LinearImplicitSystem& system = *my_system;

  DofMap& dof_map =  system.get_dof_map();

  const NumericVector<double>&  solution = *(system.solution);

  const unsigned int var = system.variable_number("T");

  FEType fe_type = dof_map.variable_type(var);

  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true)); 

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  
  const std::vector<std::vector<RealGradient> >&  dphi = fe->get_dphi();

  vector<Point> points(np);

  FEInterface::inverse_map(dim, fe_type, elem, p, points);
 
  fe->reinit(elem, &points);

  vector<unsigned int> dof_indices;

  dof_map.dof_indices (elem, dof_indices);  

  const unsigned int n_dofs   = dof_indices.size();

        

  Tensor2Sym kappa; 
  init_heat_model(elem);
  heat_model->get_thermal_conductivity(kappa);


  for (unsigned int n = 0; n < np; n++)
  {


    double T = 0.0;
    std::vector<double> Jq(3);
    Jq.clear();

    for (unsigned int p = 0; p < n_dofs; p++)
    {

      T  += phi[p][n] * solution(dof_indices[p]);
      
      for (unsigned i = 0; i<dim; ++i)
      {
        for (unsigned j = 0; j<dim; ++j)
	{

	  double kappa_value = 0.0;
	  if (i < j) 
	    kappa_value = kappa(j+1, i+1);
	  else
	    kappa_value = kappa(i+1, j+1);
	  
	  Jq[i] -= kappa_value *  dphi[p][n](j) *  solution(dof_indices[p]);
        

	}
      }
      
    }


     if (ids.count(TEMPERATURE))
     {
      values[n][TEMPERATURE] = T;
      // std::cout<<T<<std::endl;

     }

     if (ids.count(JQX))
      values[n][JQX] = Jq[0];

     if (ids.count(JQY))
       values[n][JQY] = Jq[1];

     if (ids.count(JQZ))
       values[n][JQZ] = Jq[2];

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
   
  unsigned int n_vars = 0;  
  
  std::vector<unsigned int> W;

  const unsigned int nn  = mesh->n_active_elem();
  const unsigned int dim = mesh->mesh_dimension();
  legend.resize(variables.size());


  int PF = -1;
  if (variables.find("PowerFlux") != varend)
  {



     unsigned int k = 0;

     PF = n_vars;
   
     W.push_back(n_vars);

     legend.resize(legend.size() + dim);

     switch (dim)
     {
     case 3:
       legend[W[k] + 2] = heat_legend + "_z";
       n_vars++;
     case 2:
       legend[W[k] + 1] =  heat_legend + "_y";
       n_vars++;
       legend[W[k] + dim] = "mod" + heat_legend ;
       n_vars++;
     default:
       legend[W[k] ] = heat_legend + "_x";
       n_vars++;
     }
     ++k;
     //Other power fluxes

    const Device& device = *(_device);
    
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
      std::vector<std::string> flux_legend  =
	heat_model->get_heat_source_model(ids[i])->get_flux_legend(); 
     

 
      for (unsigned int n = 0;  n < flux_legend.size(); ++n )
      {

	W.push_back(n_vars);
        legend.resize(legend.size() + dim);

	switch (dim)
	{
	case 3:
	  legend[W[k] + 2] = flux_legend[n] + "_z";
	  n_vars++;
	case 2:
	  legend[W[k] + 1] =  flux_legend[n] + "_y";
	  n_vars++;
	  legend[W[k] + dim] = "mod" + flux_legend[n] ;
	  n_vars++;
	default:
	  legend[W[k] ] = flux_legend[n] + "_x";
	  n_vars++;
	}
	++k;
	
      }
      
    }
   
    if ( k > 1)
    {
      W.push_back(n_vars);
      legend.resize(legend.size() + dim);
      
      switch (dim)
      {
      case 3:
	legend[W[k] + 2] =  "W_z";
	n_vars++;
      case 2:
	legend[W[k] + 1] =  "W_y";
	n_vars++;
	legend[W[k] + dim] = "modW";
	n_vars++;
      default:
	legend[W[k] ] = "W_x";
	n_vars++;
      }
      ++k;
    }
  }

  int HS = -1;
  if (variables.find("HeatSource") != varend)
  {
    const Device& device = *(_device);

    HS = n_vars;

    std::map<ID, std::map<ID,std::string> > heat_source_ids;
    
    MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();

    MeshBase::const_element_iterator it_end =    mesh->active_local_elements_end(); 

    //assert(it_end != mesh->active_local_elements_end());
    
    HeatModel* heat_model = NULL;
    
    const Elem* elem = *it;
    
    ID subdomain = elem->subdomain_id();
    
    heat_model= dynamic_cast<HeatModel*>(
					 device.get_material(subdomain)->get_model(get_id()));
    
    nm = heat_model->get_heat_source_IDs(ids);
    
  
    for (int i = 0; i < nm; i++)
    {
         
      std::vector<std::string> source_legend  =
 	heat_model->get_heat_source_model(ids[i])->get_source_legend(); 
      
      unsigned int legsize = source_legend.size();

      legend.resize(legend.size() + legsize);

      for (unsigned int n = 0;  n < legsize; ++n )
      {
        
	legend[n_vars]=source_legend[n];
        
	n_vars++;
	 
      }
        
    }
  }

  legend.resize(n_vars);

  //for (unsigned int n=0; n<legend.size(); ++n)
  //std::cout<<legend[n]<<std::endl;

 
 
    
  results.resize(nn * n_vars,0.0);

  LinearImplicitSystem& system = *my_system;

  const unsigned int  var = system.variable_number("T");

  DofMap& dof_map =  system.get_dof_map();

  FEType fe_type = dof_map.variable_type(var);
 
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true)); 

  QGauss qrule(dim, libMeshEnums::CONSTANT);

  fe->attach_quadrature_rule(&qrule);

  const vector<vector<RealGradient> >& dphi = fe->get_dphi();
 
  std::vector<unsigned int> dof_indices;
 
  
  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
  
 
  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  { 

    const Elem* elem = *it;

    fe->reinit(elem);

    dof_map.dof_indices (elem, dof_indices);

    unsigned int n_dofs = dof_indices.size();

    unsigned int id = n_vars * elem_number;
   
    std::vector<Point> _node(1);

    _node[0]=(elem->centroid());

    ID subdomain = elem->subdomain_id();
    
    const Material* mat = _device->get_material(subdomain);
    
    HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
    
    heat_model->set_element(elem);     
    
    heat_model->re_init();

    if (PF != -1)
    {
  
      double Px_tot = 0.0;
      double Py_tot = 0.0;
      double Pz_tot = 0.0;
     
      std::vector< std::map< ID, double > > jq_solution;

      get_solution(elem,_node,JQ_var,jq_solution); 


      //Thermal flux
      unsigned int k = 0;  

      double Px = jq_solution[0].find(JQX)->second;
      double Py = jq_solution[0].find(JQY)->second;
      double Pz = jq_solution[0].find(JQZ)->second;
      
      Px_tot += Px;
      Py_tot += Py; 
      Pz_tot += Pz;

      // std::cout<<"Pq"<<Px<<std::endl;

      switch (dim)
      {
      case 3:
	results[id + W[k] + 2] = Pz;
      case 2:
	results[id + W[k] + 1] = Py;
	results[id + W[k] + dim] = sqrt(Px * Px + Py * Py + Pz * Pz);
      default:
	results[id + W[k] ] = Px;
      }

      ++k;
      
      
      //Other power flux


      std::vector<std::vector<RealGradient> > power_flux;
      
      power_flux.clear();
      
         
      
      for (int i = 0; i < nm; i++)
      {
	
	heat_model->get_heat_source_model(ids[i])->get_power_fluxes(_node,elem,power_flux,false);

	
	for (unsigned int n_p = 0 ; n_p<power_flux[0].size() ; ++n_p)
	{ 
	  
	  double Px = power_flux[0][n_p](0);
	  double Py = power_flux[0][n_p](1);
	  double Pz = power_flux[0][n_p](2);


	  Px_tot += Px;
	  Py_tot += Py;
	  Pz_tot += Pz;

	  switch (dim)
	  {
	  case 3:
	    results[id + W[k] + 2] = Pz;
	  case 2:
	    results[id + W[k] + 1] = Py;
	    results[id + W[k] + dim] = sqrt(Px * Px + Py * Py + Pz * Pz);
	  default:
	    results[id + W[k] ] = Px;
	  }
          ++k;

	}
      }
      
      if (k>1)
      {
       

	switch (dim)
	{
	case 3:
	  results[id + W[k] + 2] = Pz_tot;
	case 2:
	  results[id + W[k] + 1] = Py_tot;
	  results[id + W[k] + dim] = sqrt(Px_tot * Px_tot + Py_tot * Py_tot + Pz_tot * Pz_tot);
	default:
	  results[id + W[k] ] = Px_tot;
	}
	
      }


    }
    
    if (HS != -1)
    {
            
  
  
      std::vector<double> total_heat_source(nm);

      total_heat_source.clear();

      total_heat_source.resize(nm);
	
      unsigned int k = 0;
     
      for (int i = 0; i < nm; i++)
      {
	  
	std::vector<std::vector<double> > heat_source;
       
  
	heat_model->get_heat_source_model(ids[i])->get_heat_sources(_node,elem,heat_source);
               
	for (unsigned int n_s = 0 ; n_s<heat_source[0].size() ; ++n_s)
	{ 
	  
          results[id + HS + k] = heat_source[0][n_s];
	  
          ++k;
	  
	}
      }
	
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

;

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

void
MacroHeatBalance::build_integrated_quantities_description(
    const std::set<std::string>& names,
    std::vector<std::string>& legend,
    std::vector<std::string>& description)
{
  const set<string>::const_iterator varend(names.end());

  if (names.find("Power") != varend)
  {
    legend.resize(1);

    description.resize(1);

    const unsigned int dim = mesh->mesh_dimension();
  
    ostringstream s;
    s << "Power Dissiapted. Units W";
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
MacroHeatBalance::build_integrated_quantities(const set<string>& names,
    vector<double>& values)
{
  const set<string>::const_iterator varend(names.end());

  if (names.find("Power") != varend)
  {
    
    calculate_power_surfint();

     values.resize(1,_power);


  }
}




void
MacroHeatBalance::calculate_power_surfint(void)
{

  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return;

  LinearImplicitSystem& system = *my_system;
  
  const unsigned int  var = system.variable_number("T");

  const unsigned int dim = mesh->mesh_dimension();

  DofMap& dof_map =  system.get_dof_map();

  FEType fe_type = dof_map.variable_type(var);
 
  AutoPtr<FEBase> fe_face (build_finite_element(dim, fe_type)); 

  // libMeshEnums::Order integration_order;
  //if (dim == 1)
  //  integration_order = libMeshEnums::CONSTANT;
  //else
  //  integration_order = _options.integration_order;
  
  //QGauss qface(dim - 1,integration_order);

  QGauss qface(dim - 1,FIFTH);

  fe_face->attach_quadrature_rule(&qface);

  const vector<Real>& JxW = fe_face->get_JxW();
  
  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe_face->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe_face->get_phi();

  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe_face->get_dphi();

  // the face normals
  const vector<Point>& face_normals = fe_face->get_normals();

  std::vector<unsigned int> dof_indices;


  MeshBase::const_element_iterator el =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =     mesh->active_local_elements_end();

  const SimulationEnvironment& env = get_environment();

  _power = 0.0;     

  for ( ; el != end_el ; ++el) 
  {
    
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    ID subdomain = elem->subdomain_id();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices,var);
  
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      ElementSide side(top_parent, s);
      
      if (env.is_boundary(side))
      {
	
	init_heat_model(elem);
	
	//	assert(_heat_model != NULL);

	// only for dim > 1 we need to integrate
        if (dim > 1)
        {
	  fe_face->reinit(elem, s);

	  std::vector< std::map< ID, double > > jq_solution;

	  get_solution(elem,q_point,JQ_var,jq_solution); 
     
          int phi_size = phi.size();
    
      	  RealGradient P;

          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
     
	    P(0) = jq_solution[qp].find(JQX)->second;
	    P(1) = jq_solution[qp].find(JQY)->second;
	    P(2) = jq_solution[qp].find(JQZ)->second;
           
	    _power += JxW[qp] * P * face_normals[qp]; 
              
          } 

        }
        else
        {

          vector<Point> v(1, elem->point(s));
          fe_face->reinit(elem, &v);

    	  std::vector< std::map< ID, double > > jq_solution;

	  get_solution(elem,v,JQ_var,jq_solution); 
          
	  double P = jq_solution[0].find(JQX)->second;

          // what is the outer normal in this point??
          // Idea: if x(s) > x(centroid), normal is +1
          //       else it is -1
          double x_c = elem->centroid()(0);
          double x_s = elem->point(s)(0);
          if (x_s < x_c)
          {
            P = - P;
          }
	  // std::cout<<P<<std::endl;
          _power += P;

        }
      }
    } // end loop over elem sides
  } // end loop over elements

 

}



