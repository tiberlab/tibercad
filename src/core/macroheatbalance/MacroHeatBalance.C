// $Id$

#include "MacroHeatBalance.h"
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
#include "ThermalSurfaceResistance.h"
#include "ThermalSurfaceConductance.h"
#include "FluxContact.h"

//To Canc
#include "PhononModel.h"
//

//To Cut
extern "C" void zheev_(char *jobz, char *uplo, int *n, complex<double> *a,int *lda,
 double *w, complex<double> *work, int *lwork, double *rwork, int *info);
//
extern "C" void dsyev_(char* jobz,  char* uplo,int* n, double* a,int* lda, double* w,double* work,int* lwork, int* info );


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

  //To Canc
   myopts.mode = sim_opt.get_option("SimulationMode","none");


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
  

  my_system = TiberLinearSystem::create(get_equation_systems(),
      get_equation_system_name(), get_solver_options());


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
  //to be cut

  //
   
  //------init is done---------------------------------------------------------------------//

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
  
 
  //To Canc
  std::string mode = options.get_option("SimulationMode","none");

  if (mode.compare("PhononDispersion") == 0 )
  {

    std::cout<<"PhononModel...";
    
    PhononModel* model = dynamic_cast<PhononModel*> ( PhysicalModelInterface::create("phonon",options) );
    
   if (model == NULL) 
     throw ModelErrorException("MacroHeatBalance: PhononDispersion physical model is not created" );
     return model;      
  }
  else
  {
  HeatModel* model = dynamic_cast<HeatModel*> ( PhysicalModelInterface::create("thermal",options) );

  if (model == NULL) 
  throw ModelErrorException("MacroHeatBalance: Thermal physical model is not created" );

  return model;      
  }
  // 

  
  //To UnComm
  // HeatModel* model = dynamic_cast<HeatModel*> ( PhysicalModelInterface::create("thermal",options) );
  //if (model == NULL) 
  // throw ModelErrorException("MacroHeatBalance: Thermal physical model is not created" );
  //



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
 
 
  std::vector<unsigned int> dof_indices;
  
  DenseMatrix<Number>  Ke;

  DenseVector<Number>  Fe;

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  
  //Model Variables

  Tensor2Sym kappa; 



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


    ID subdomain = elem->subdomain_id();
    
    const Material* mat = _device->get_material(subdomain);
    
    HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
    
    heat_model->set_element(elem);     

    heat_model->set_side(-1);

    heat_model->re_init();   


    std::vector<double> heat_source;
    heat_model->get_total_heat_source(q_point,heat_source);

    std::vector<RealGradient> flux_power;
    heat_model->get_total_power_flux(q_point,flux_power);

    heat_model->get_thermal_conductivity(kappa);

    

 
    for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
    { // loop over test function
      
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
	
	//Fe(p1) +=JxW[qp] * heat_source[qp] * phi[p1][qp];

	Fe(p1) += JxW[qp] * flux_power[qp] * dphi[p1][qp];
	
      }//end Loop over quadrature points  
    } // end loop over test function
    
    //Boundary conditions and source
    
    //The loop over element is the only loop that is surviving at this point
    
    const unsigned int num_sides = elem->n_sides();
    
    for (unsigned int side = 0; side<num_sides; side++)
    {
      
      const ElementSide elside(elem->top_parent(), side);
      
      heat_model->set_side(side);
	
      heat_model->re_init();
	
      fe_face->reinit(elem,side);
	
      //Source. Must be before boundary condition  
      heat_model->get_total_power_flux(qface_point,flux_power);
	
        for (unsigned int qp=0; qp < qface.n_points(); qp++)
       {
      	for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
      	{
      	  double Fe_surf = JxW_face[qp] * phi_face[p1][qp] * flux_power[qp] * normal[qp];
	 
      	   Fe(p1) -= Fe_surf;
	  
      	}
       }
     
      Boundary* bd = se.get_boundary(elside);

      if (bd != NULL)
      { 
	if (bd->get_boundary_properties( get_id() ) != NULL ) 
	{
	 

	  ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );    
	  
	  switch (contact->get_type())
	  { 
	  case  ThermalContact::Reservoir:
	    
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
	    break;
	    
	  case  ThermalContact::ThermalSurfaceResistance :
	    {
	      double r_surf =  ( dynamic_cast<ThermalSurfaceResistance*> (contact) )->get_thermal_surface_resistance();
	      double temp =  ( dynamic_cast<ThermalSurfaceResistance*> (contact) )->get_temperature();
	      
	      for (unsigned int qp=0; qp < qface.n_points(); qp++)
	      {
		for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
		{
		  double Fe_plus = 0.0; 
		  
		  Fe_plus = JxW_face[qp] * 1/r_surf * phi_face[p1][qp] * temp;
		  Fe(p1) += Fe_plus;
		  
		  for (unsigned int p2=0; p2<n_dofs; p2++) //bases test
		  {
		    double val_plus = 0.0;
		    val_plus  =  JxW_face[qp] * 1/r_surf * phi_face[p1][qp] * phi_face[p2][qp];
		    Ke(p1,p2) += val_plus;
		    
		  }// (unsigned int p2=0; p2<n_dofs; p2++)
		}//for (unsigned int p1=0; p1<n_dofs; p1++) 
	      }// for (unsigned int qp=0; qp < qface.n_points(); qp++)
	    }
	    break;
	    
	    
	  case  ThermalContact::ThermalSurfaceConductance :
	    {
	      double g_surf =  ( dynamic_cast<ThermalSurfaceConductance*> (contact) )->get_thermal_surface_conductance();
	      double temp =  ( dynamic_cast<ThermalSurfaceConductance*> (contact) )->get_temperature();
	      
	      for (unsigned int qp=0; qp < qface.n_points(); qp++)
	      {
		for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
		{
		  double Fe_plus = 0.0; 
		  
		  Fe_plus = JxW_face[qp] * g_surf * phi_face[p1][qp] * temp;
		  Fe(p1) += Fe_plus;
		  
		  for (unsigned int p2=0; p2<n_dofs; p2++) //bases test
		  {
		    double val_plus = 0.0;
		    val_plus  =  JxW_face[qp] * g_surf * phi_face[p1][qp] * phi_face[p2][qp];
		    Ke(p1,p2) += val_plus;
		    
		  }// (unsigned int p2=0; p2<n_dofs; p2++)
		}//for (unsigned int p1=0; p1<n_dofs; p1++) 
	      }// for (unsigned int qp=0; qp < qface.n_points(); qp++)
	    }
	    break;
	    
	    
	  case  ThermalContact::FluxContact :
	    {        
	      double heat_flux =  ( dynamic_cast<FluxContact*> (contact) )->get_heat_flux();
	      
	      for (unsigned int qp=0; qp < qface.n_points(); qp++)
	      {
		
		for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T
		{
		  double fe_plus = 0.0;
		  fe_plus =  JxW_face[qp] * heat_flux * phi_face[p1][qp];
		  
		  Fe(p1) += fe_plus;
		  
		}//for (unsigned int p1=0; p1<n_dofs; p1++) 
		
	      }// for (unsigned int qp=0; qp < qface.n_points(); qp++)
	    }
	    break;
	    
	  }//switch 
	  
	}//  if (bd->get_boundary_properties( get_id() ) != NULL )
	
	
      }// if (is_boundary != NULL)
      
    }// for (unsigned int side = 0; side<num_sides; side++)
    

    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices); 
  
    
  } //End Loop over elements
  //   system.matrix->print_matlab("Matr.m");
  //   system.rhs->print();
    
} //do assembly




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

  TiberLinearSystem& system = *my_system;

  DofMap& dof_map = system.get_dof_map();

  const NumericVector<double>& solution = *(system.solution);

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

        

  Tensor2Sym kappa; 
  
  
  ID subdomain = elem->subdomain_id();
  const Material* mat = _device->get_material(subdomain);
  HeatModel* heat_model = (dynamic_cast<HeatModel*>(mat->get_model(get_id())));
  heat_model->set_element(elem);     
  heat_model->set_side(-1);
  heat_model->re_init(); 
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
       values[n][TEMPERATURE] = T;
     
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

  if (myopts.mode.compare("PhononDispersion") == 0 )
  {
    build_elemental_results_phonon_disp(variables,results,legend);
  } 
  else
  {
    build_elemental_results_original(variables,results,legend);
  }
}


void
MacroHeatBalance::build_elemental_results_original(const std::set<std::string>& variables,
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

 
  
    const Device& device = *(_device);
    
    //  HS = n_vars;
    
    std::map<ID, std::map<ID,std::string> > heat_source_ids;
    
    MeshBase::const_element_iterator it_temp =    mesh->active_local_elements_begin();
    
    MeshBase::const_element_iterator it_end =    mesh->active_local_elements_end(); 
    
    //assert(it_end != mesh->active_local_elements_end());
    
    HeatModel* heat_model = NULL;
    
    const Elem* elem = *it_temp;
    
    ID subdomain = elem->subdomain_id();
    
    heat_model= dynamic_cast<HeatModel*>(
					 device.get_material(subdomain)->get_model(get_id()));
    
    nm = heat_model->get_heat_source_IDs(ids);
    
   
    std::vector<std::set<ID> > source_index(nm); 
    
    for (int i = 0; i < nm; i++)
    {	  
      
      std::map<ID,std::string> source_legend =  
	heat_model->get_heat_source_model(ids[i])->get_source_legend(variables);
      
      std::map<ID,std::string>::iterator leg(source_legend.begin());
      std::map<ID,std::string>::iterator leg_end(source_legend.end());
      
      for (;leg != leg_end; leg++)
      {
	
	legend.resize(legend.size() + 1);
	legend[n_vars]=leg->second;
	source_index[i].insert(leg->first);
        n_vars++;
      }    
    }
  
  if (variables.count("TotalHeat")  ||
      variables.count("HeatSource") ||
      variables.count("thermal"))
  {
    legend.resize(legend.size() + 1);
    legend[n_vars]="TotalHeat";
    n_vars++;
  }

  int HS = -1;
  if (n_vars>0)
    HS=0;

    ID PF_temp = n_vars;

     unsigned int k = 0;

     if (variables.count("thermal") ||
         variables.count("ThermalFlux")      ||
         variables.count("PowerFlux") )
     {
     
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
    }

   
    //Other fluxes
     std::vector<std::set<ID> > flux_index(nm);
     for (int i = 0; i < nm; i++)
     {
       std::map<ID,std::string> flux_legend  =
	 heat_model->get_heat_source_model(ids[i])->get_flux_legend(variables); 

       std::map<ID,std::string>::iterator leg(flux_legend.begin());
       std::map<ID,std::string>::iterator leg_end(flux_legend.end());
       
       for (; leg != leg_end; leg++)
       {
	 W.push_back(n_vars);
	 legend.resize(legend.size() + dim);
	 
         flux_index[i].insert(leg->first);
	 std::string label = leg->second;
	 
	 switch (dim)
	 {
	 case 3:
	   legend[W[k] + 2] = label + "_z";
	   n_vars++;
	 case 2:
	   legend[W[k] + 1] =  label + "_y";
	   n_vars++;
	   legend[W[k] + dim] = "mod" + label;
	   n_vars++;
	 default:
	   legend[W[k] ] = label + "_x";
	   n_vars++;
	 }
	 ++k;
       } 
     }
     
     if (variables.count("thermal") ||
	 variables.count("PowerFlux")      ||
	 variables.count("TotalFlux") )
     {
       
       W.push_back(n_vars);
       legend.resize(legend.size() + dim);
       
       switch (dim)
       {
       case 3:
	 legend[W[k] + 2] = "W_z";
	 n_vars++;
       case 2:
	 legend[W[k] + 1] = "W_y";
	 n_vars++;
	 legend[W[k] + dim] ="modW";
	 n_vars++;
       default:
	 legend[W[k] ] ="W_x";
	 n_vars++;
       }
       ++k;
       
     }
    
   
     int PF = -1;  
     if (n_vars>PF_temp)    
       PF = PF_temp;
      

    int Kappa = -1;
    int Kappa_xx = -1;
    int Kappa_zz = -1;
    if (variables.count("thermal") ||
        variables.count("LatticeThermalCond") )
    {
      Kappa = 0;
      Kappa_xx = n_vars;
      legend.resize(legend.size() + 1);
      legend[n_vars]="kappa_xx";
      n_vars++;

      Kappa_zz = n_vars;
      legend.resize(legend.size() + 1);
      legend[n_vars]="kappa_zz";
      n_vars++; 
    }



  legend.resize(n_vars);
  
  results.resize(nn * n_vars,0.0);

  TiberLinearSystem& system = *my_system;

  const unsigned int  var = system.variable_number("T");

  DofMap& dof_map =  system.get_dof_map();

  FEType fe_type = dof_map.variable_type(var);
 
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true)); 

  QGauss qrule(dim, libMeshEnums::CONSTANT);

  fe->attach_quadrature_rule(&qrule);

  const vector<vector<RealGradient> >& dphi = fe->get_dphi();
 
  std::vector<unsigned int> dof_indices;

  Tensor2Sym kappa; 
  
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

    heat_model->set_side(-1);

    heat_model->re_init();

    heat_model->get_thermal_conductivity(kappa);

    if (HS != -1)
    {
      
      unsigned int k = 0;

      for (int i = 0; i < nm; i++)
      {
        std::vector<std::map<ID, double> > heat_sources;  
	
        heat_model->get_heat_source_model(ids[i])->get_heat_sources(_node,source_index[i],heat_sources);

	std::map<ID,double>::iterator  it_s(heat_sources[0].begin());
	std::map<ID,double>::iterator  it_end(heat_sources[0].end());


	for (;it_s != it_end; it_s++)         
        {
          results[id + HS + k] = it_s->second;
          ++k; 
	}

      }

      //  //Include Total Heat source 
      if (variables.count("TotalHeat")  ||
	  variables.count("HeatSource") ||
	  variables.count("thermal"))
      {
        std::vector< double > total_heat_source; 

	heat_model->get_total_heat_source(_node,total_heat_source);
	results[id + HS + k] =  total_heat_source[0];

      }
	
    } //if (HS != -1)

    
    if (PF != -1)
    {
      unsigned int k = 0; 

      std::vector< std::map< ID, double > > jq_solution;
      get_solution(elem,_node,JQ_var,jq_solution); 
 	
      double Pqx = jq_solution[0].find(JQX)->second;
      double Pqy = jq_solution[0].find(JQY)->second;
      double Pqz = jq_solution[0].find(JQZ)->second;

      if (variables.count("thermal") ||
	  variables.count("ThermalFlux")      ||
	  variables.count("PowerFlux") )
      {
      
	switch (dim)
	{
	case 3:
	  results[id + W[k] + 2] = Pqz;
	case 2:
	  results[id + W[k] + 1] = Pqy;
	  results[id + W[k] + dim] = sqrt(Pqx * Pqx + Pqy * Pqy + Pqz * Pqz);
	default:
	  results[id + W[k] ] = Pqx;
	}
	
	++k;
      }
      
      //Other power flux
      std::vector<std::map<ID,RealGradient> > power_flux;
      
      for (int i = 0; i < nm; i++)
      {
	heat_model->get_heat_source_model(ids[i])->get_power_fluxes(_node,flux_index[i],power_flux);


       	std::map<ID,RealGradient>::iterator  it_s(power_flux[0].begin());
	std::map<ID,RealGradient>::iterator  it_end(power_flux[0].end());
 
	for (;it_s != it_end; it_s++)         
        {

          double Px = (it_s->second) (0);
	  double Py = (it_s->second) (1);
	  double Pz = (it_s->second) (2);

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
       
      }//loop over models

      if (variables.count("TotalFlux")  ||
	  variables.count("thermal")    ||
	  variables.count("PowerFlux") )
      {
	
	std::vector<RealGradient > total_power_flux; 

	heat_model->get_total_power_flux(_node,total_power_flux);

	double Px_tot = total_power_flux[0](0) + Pqx; 
	double Py_tot = total_power_flux[0](1) + Pqy; 
	double Pz_tot = total_power_flux[0](2) + Pqz; 

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

    if (Kappa != -1)
    {
      
      results[id + Kappa_xx] = kappa(1,1);
      results[id + Kappa_zz] = kappa(3,3); 

    }



    elem_number++;
  } //over element

  results.resize(elem_number * n_vars);
}


void
MacroHeatBalance::build_elemental_results_phonon_disp(const std::set<std::string>& variables,
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
  const Device& device = *(_device);

  
  int PD = -1;
  if (variables.count("RamanShift") ||
      variables.count("PhononVariables")  )
  {
    PD = n_vars;
    legend.resize(legend.size() + 1);
    legend[n_vars]="dw_1";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="dw_2";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="dw_3";
    n_vars++;
  }

  
  //To cut
  int PFS = -1;
  if (variables.count("PhononEnergy") ||
      variables.count("PhononVariables") )
  {
    PFS = n_vars;
    legend.resize(legend.size() + 1);
    legend[n_vars]="w_1";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="w_2";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="w_3";
    n_vars++;
  }
  ///

  int EV = -1;
  int EV1 = -1;
  int EV2 = -1;
  int EV3 = -1;
  if (variables.count("PhononPolarization")  ||
      variables.count("PhononVariables")  )
  {
    EV = 1;
    legend.resize(legend.size() + dim);
    EV1 = n_vars;
    switch (dim)
    {
      case 3:
        legend[EV1 + 2] = "E1_z";
        n_vars++;
      case 2:
        legend[EV1 + 1] = "E1_y";
        n_vars++;
        legend[EV1 + dim] = "modE1";
        n_vars++;
      default:
        legend[EV1] = "E1_x";
        n_vars++;
    }

    legend.resize(legend.size() + dim);
    EV2 = n_vars;
    switch (dim)
    {
      case 3:
        legend[EV2 + 2] = "E2_z";
        n_vars++;
      case 2:
        legend[EV2 + 1] = "E2_y";
        n_vars++;
        legend[EV2 + dim] = "modE2";
        n_vars++;
      default:
        legend[EV2] = "E2_x";
        n_vars++;
    }

    legend.resize(legend.size() + dim);
    EV3 = n_vars;
    switch (dim)
    {
      case 3:
        legend[EV3 + 2] = "E3_z";
        n_vars++;
      case 2:
        legend[EV3 + 1] = "E3_y";
        n_vars++;
        legend[EV3 + dim] = "modE3";
        n_vars++;
      default:
        legend[EV3] = "E3_x";
        n_vars++;
    }

  }

 //To cut
  int OA = -1;
  if (variables.count("OverAll") ||
      variables.count("PhononVariables") )
  {
    OA = n_vars;
    legend.resize(legend.size() + 1);
    legend[n_vars]="ParallelOverAll";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="CrossOverAll";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="NoPolOverAll";
    n_vars++;
    
  }
  ///

  //To cut
  int RI = -1;
  if (variables.count("RamanIntensity") ||
      variables.count("PhononVariables")  )
  {
    RI = n_vars;

    legend.resize(legend.size() + 1);
    legend[n_vars]="ParallelRamanIntensity_1";
    n_vars++;

    legend.resize(legend.size() + 1);
    legend[n_vars]="ParallelRamanIntensity_2";
    n_vars++;
    
    legend.resize(legend.size() + 1);
    legend[n_vars]="ParallelRamanIntensity_3";
    n_vars++;
  
    
    legend.resize(legend.size() + 1);
    legend[n_vars]="CrossRamanIntensity_1";
    n_vars++;

    legend.resize(legend.size() + 1);
    legend[n_vars]="CrossRamanIntensity_2";
    n_vars++;
    
    legend.resize(legend.size() + 1);
    legend[n_vars]="CrossRamanIntensity_3";
    n_vars++;
  
    legend.resize(legend.size() + 1);
    legend[n_vars]="NoPolRamanIntensity_1";
    n_vars++;

    legend.resize(legend.size() + 1);
    legend[n_vars]="NoPolRamanIntensity_2";
    n_vars++;
    
    legend.resize(legend.size() + 1);
    legend[n_vars]="NoPolRamanIntensity_3";
    n_vars++;

  }
  ///


    legend.resize(n_vars);
    
    results.resize(nn * n_vars,0.0);
 
  
  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
  
  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  { 
   unsigned int id = n_vars * elem_number;
    const Elem* elem = *it;

    if (PFS != -1)
    {
      std::vector<double> sol = PD_full_sol[elem];
      results[id + PFS]   = sol[0];
      results[id + PFS+1] = sol[1]; 
      results[id + PFS+2] = sol[2]; 
    }

    if (PD != -1)
    {
      std::vector<double> sol = PD_sol[elem];
      results[id + PD]   = sol[0];
      results[id + PD+1] = sol[1]; 
      results[id + PD+2] = sol[2]; 
    }

    if (EV != -1)
    {
      std::vector< std::vector <double> >  sol = sol_eignvectors[elem];

      //std::cout<< sol[0][0]<<std::endl;
 
      RealGradient eign(0);




      eign(0) = sol[0][0];
      eign(1) = sol[1][0];
      eign(2) = sol[2][0];

      switch (dim)
      {
        case 3:
          results[id + EV1 + 2] = eign(2);
        case 2:
          results[id + EV1 + 1] = eign(0);
          results[id + EV1 + dim] = eign.size();
        default:
          results[id + EV1] = eign(1); 
      
      }

 
      eign(0) = sol[0][1];
      eign(1) = sol[1][1];
      eign(2) = sol[2][1];

      switch (dim)
      {
        case 3:
          results[id + EV2 + 2] = eign(2);
        case 2:
          results[id + EV2 + 1] = eign(0);
          results[id + EV2 + dim] = eign.size();
        default:
          results[id + EV2] = eign(1); 
      
      }

 
      eign(0) = sol[0][2];
      eign(1) = sol[1][2];
      eign(2) = sol[2][2];

      switch (dim)
      {
        case 3:
          results[id + EV3 + 2] = eign(2);
        case 2:
          results[id + EV3 + 1] = eign(0);
          results[id + EV3 + dim] = eign.size();
        default:
          results[id + EV3] = eign(1); 
      
      }




    }

    if (OA != -1)
    {
      std::vector<double> sol = OverAll[elem];
      results[id + OA]    = sol[0];
      results[id + OA+1]  = sol[1]; 
      results[id + OA+2]  = sol[2];
          }
  

    if (RI != -1)
    {
      std::vector<std::vector <double> > sol = Intensity[elem];
      results[id + RI]   = sol[0][0];
      results[id + RI+1] = sol[0][1];
      results[id + RI+2] = sol[0][2];

      results[id + RI +3] = sol[1][0];
      results[id + RI +4] = sol[1][1];
      results[id + RI +5] = sol[1][2];
 
      results[id + RI +6] = sol[2][0];
      results[id + RI +7] = sol[2][1];
      results[id + RI +8] = sol[2][2];
     
 
    }
  

    elem_number++;
  } //over element

  results.resize(elem_number * n_vars);
}


//----------------------------------------------------------------------------------//
void MacroHeatBalance::build_nodal_results (const std::set< std::string > &variables, 
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
}
void
MacroHeatBalance::build_integrated_quantities_description(
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
    s << "Power Dissiapated. Units W";
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


  if (names.count("PowerDissipated"))
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

  TiberLinearSystem& system = *my_system;
  
  const unsigned int  var = system.variable_number("T");

  const unsigned int dim = mesh->mesh_dimension();

  DofMap& dof_map =  system.get_dof_map();

  FEType fe_type = dof_map.variable_type(var);
 
  AutoPtr<FEBase> fe_face (build_finite_element(dim, fe_type)); 
 
  libMeshEnums::Order integration_order;

  //if (dim == 1)
  //integration_order = libMeshEnums::CONSTANT;
  //else
  //  integration_order = _options.integration_order;
  
  //QGauss qface(dim - 1,integration_order);

  QGauss qface(dim - 1, FIFTH);

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
      
      if (env.is_outer_boundary(side))
      {

	ID subdomain = elem->subdomain_id();
	
	const Material* mat = _device->get_material(subdomain);
	
	HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
	
	heat_model->set_element(elem);     
	
	heat_model->set_side(-1);
	
	heat_model->re_init();

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
    } // end loop over elem sides
  } // end loop over elements
}



void
MacroHeatBalance::do_print_info(void)
{

  string space("  ");
  //cout << space << "linear solver is: petsc" <<std::endl;
   
}

//-------------------------------------------------------------------------------//
void  MacroHeatBalance::do_solve()
{
   
  parse_options();
 
  //To Canc;
  std::cout<<"SimulationMode->"<<"  "<<myopts.mode<<std::endl;

  if (myopts.mode.compare("PhononDispersion") == 0 )
  {
   
    solve_phonon_dispersion();		 
  }
 else
 { 
    static_this = this;  
    my_system->set_options(get_solver_options());
    my_system->solve();

  }
 
  
}

void 
MacroHeatBalance::diagonalize_complex(void)
 {


   std::cout<<"Diagonalize...start"<<std::endl;


   int  n = 3;
   
   char jobz = 'V';
   char uplo = 'U';
   int info;
   int lwork = 3*n-1;
   
   double *b;
   complex<double> *work = new complex<double>[lwork]; //The work array to be used by zheev an its size
   double *rwork = new double[3*n-2];
   complex<double> *a = new complex<double>[n*n];

   for (int j=0; j < n; j++){
     for (int i=0; i < n; i++){
        double value;
       if (i < j) 
	 value = dynamical_matrix(j+1, i+1);
       else
	 value = dynamical_matrix(i+1, j+1);

   
       complex<double> s(value,0.0);
       a[n*j+i] =s;//passed as rods
     }
   }
    
   zheev_(&jobz, &uplo, &n, a, &n, b, work, &lwork, rwork, &info);//perform the diagonalization
     
   delete [] work;
   delete [] rwork;
   
   //sorting
   {  
     complex<double> ctemp;
     
     for(int j = 0; j < n; j++){
       int it = j;
       double temp = b[j];
      
       for(int i = j; i < n; i++){
	 if(b[i] > temp){
	   temp = b[i];
	   it = i;
	 }
       }

       b[it] = b[j];
       b[j] = temp;
      
       for(int k = 0; k < n; k++){
	 ctemp = a[n*j+k];
	 a[n*j+k] = a[n*it+k];
	 a[n*it+k] = ctemp;
       }
     }
     
   }
      //empty the autovector matrix
   _eignvectors.resize(n);
   E.resize(n);
   for (unsigned int k = 0; k<n; k++)
      _eignvectors[k].resize(n,0.0);
   
   //Put the value in the new style matrix
   for(int i=0; i < n; i++)
   {
     for(int j = 0; j < n; j++)
     {
        _eignvectors[i][j] = a[n*i+j].real();//passed as lines
     }
     E[i] = b[i];
   }
   
   delete [] a;
   delete [] b;

}          
//

void 
MacroHeatBalance::diagonalize_double(void)
 {


   // std::cout<<"Diagonalize...start"<<std::endl;

   //int n = 2;   
    //  D.resize(n);
//   for (unsigned int i = 0; i<n; i++)
//   {
//    D[i].resize(n,0.0);
//   }
 
  
//   D[0][0] = 1.0;
//    D[0][1] = 0.0;
//    D[1][0] = 0.0;
//    D[1][1] = 2.0; 

   //int  n = D[0].size();


   int n = 3;
   
   char jobz = 'V';
   char uplo = 'U';
   int info;
   int lwork = 3*n-1;

   double *b = new double[n];;
   double *work = new double[lwork]; //The work array to be used by zheev an its size
   double *a = new double[n*n];

   // std::cout<<"DM"<<std::endl;
   for (int j=0; j < n; j++){
    for (int i=0; i < n; i++){
       
     double value;
      if (i < j) 
        value = dynamical_matrix(j+1, i+1);
      else
        value = dynamical_matrix(i+1, j+1);

       a[n*j+i] = (double) value;
       //std::cout<< a[n*j+i] <<std::endl;
       }
    }

    dsyev_(&jobz, &uplo, &n, a, &n, b,work,&lwork,&info);

    //  std::cout<<"eignvalue:"<<std::endl;
    //std::cout<<std::sqrt(b[0])*8065.6<<std::endl;
    //std::cout<<std::sqrt(b[1])*8065.6<<std::endl;
    //std::cout<<std::sqrt(b[2])*8065.6<<std::endl;

   delete [] work;
    
   _eignvectors.resize(n);
   E.resize(n);
   for (unsigned int k = 0; k<n; k++)
      { _eignvectors[k].resize(n,0.0);}
   
   //Put the value in the new style matrix
   for(int i=0; i < n; i++)
    {
    for(int j = 0; j < n; j++)
     {
          _eignvectors[i][j] = a[n*i+j];//passed as lines
    }
     E[i] = std::sqrt(b[i]) *8065.6;
   }
   
   delete [] a;
   delete [] b;


}          
//

void MacroHeatBalance::solve_phonon_dispersion(void)
{

     std::vector<Tensor2Sym> raman_tensor;
     std::vector<Tensor1> light_polarization;
  

     mesh = & (_device->get_mesh());

     MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
     const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

      for ( ; el != end_el ; ++el)
     {   //loop over elements
       const Elem* elem = *el;

       ID subdomain = elem->subdomain_id();

       const Material* mat = _device->get_material(subdomain);

       PhononModel* phonon_model =  (  dynamic_cast<PhononModel*> (  mat -> get_model(get_id()) )  );

       phonon_model->set_element(elem);     

       phonon_model->re_init(); 

       //Full
       phonon_model->get_full_dynamical_matrix(dynamical_matrix);
       diagonalize_tensor();
       std::vector<double> full_sol  =  E;
       sol_eignvectors[elem] = _eignvectors;
          
        

       //Free
       phonon_model->get_free_dynamical_matrix(dynamical_matrix);
       diagonalize_tensor();
       std::vector<double> free_sol  =  E;
       
       std::vector<double> diff(3);
       for (unsigned int k = 0; k<3; k++)       
       {diff[k] = full_sol[k]-free_sol[k];}

       PD_full_sol[elem] = full_sol;
       PD_sol[elem] = diff; 

       phonon_model->get_raman_tensor(raman_tensor);

     phonon_model->get_light_polarization(light_polarization);
     Tensor1 e0       = light_polarization[0];
     Tensor1 es_paral = light_polarization[1];
     Tensor1 es_cross = light_polarization[2];
     Tensor1 es_nopol = light_polarization[3];

       std::vector<double> I_paral(3);
       std::vector<double> I_cross(3);
       std::vector<double> I_nopol(3);

       for(unsigned int nm =0; nm<3;nm ++)
       {
         
          Tensor2Gen D =  raman_tensor[0]*_eignvectors[0][nm] +
                          raman_tensor[1]*_eignvectors[1][nm] +
                          raman_tensor[2]*_eignvectors[2][nm];
       
         

          Tensor1 tens = D*e0;


          //Parallel
          double int_temp = es_paral*tens;
          I_paral[nm] = std::abs(int_temp) * std::abs(int_temp);

          //Cross
          int_temp = es_cross*tens;
          I_cross[nm] = std::abs(int_temp) * std::abs(int_temp);

          //NoPol
          int_temp = es_nopol*tens;
          I_nopol[nm] = std::abs(int_temp) * std::abs(int_temp);

       }
  
       //Paral
       double temp_paral = 0;
       for (unsigned int k = 0; k<3; k++)
         {temp_paral += diff[k] * I_paral[k];}
        
        double sum_paral =  I_paral[0] +  I_paral[1] + I_paral[2];
        if (sum_paral < 1e-13)
       {
            temp_paral = 0.0;
         }
       else
       {
        temp_paral /=sum_paral; 
        }


       //Cross
       double temp_cross = 0;
       for (unsigned int k = 0; k<3; k++)
         {temp_cross += diff[k] * I_cross[k];}

            
       double sum_cross =  I_cross[0] +  I_cross[1] + I_cross[2];
       if (sum_cross < 1e-13)
       {
            temp_cross = 0.0;
         }
       else
       {
        temp_cross /=sum_cross; 
        }

       //NoPol
       double temp_nopol = 0;
       for (unsigned int k = 0; k<3; k++)
         {temp_nopol += diff[k] * I_nopol[k];}

            
       double sum_nopol =  I_nopol[0] +  I_nopol[1] + I_nopol[2];
       if (sum_nopol < 1e-13)
       {
            temp_nopol = 0.0;
         }
       else
       {
        temp_nopol /=sum_nopol; 
       }


        //Write Output
        OverAll[elem].resize(3,0.0);
        OverAll[elem][0] = temp_paral;
        OverAll[elem][1] = temp_cross;
        OverAll[elem][2] = temp_nopol;

        //Intensity
        Intensity[elem].resize(3);
        Intensity[elem][0].resize(3,0.0);
        Intensity[elem][0] = I_paral;
        Intensity[elem][1].resize(3,0.0);
        Intensity[elem][1] = I_cross;
        Intensity[elem][2].resize(3,0.0);
        Intensity[elem][2] = I_nopol;

        
     }
     
}

void 
MacroHeatBalance::diagonalize_tensor(void)
 {

   int n = 3;
   
   Tensor2Gen  V;
   double landa1;
   double landa2;
   double landa3;
   dynamical_matrix.eigen(&landa1,&landa2,&landa3,&V);

      _eignvectors.resize(n);
      E.resize(n);
      for (unsigned int k = 0; k<n; k++)
      { _eignvectors[k].resize(n,0.0);}

    

    for(unsigned int i=0; i < n; i++)
    {
    for(unsigned int j = 0; j < n; j++)
     {
          
       _eignvectors[i][j] = V(i+1,j+1);

    }
    
  }

      double conv = 8065.6;

      E[0]=sqrt(landa1)*conv;
      E[1]=sqrt(landa2)*conv;
      E[2]=sqrt(landa3)*conv;
   

}          
