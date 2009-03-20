// $Id$


#include "Poisson.h"
#include "TiberLinearSystem.h"
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


#include "SimulationEnvironment.h"
#include "Material.h"
#include "Boundary.h"
#include "Reservoir.h"
#include "SimulationOptions.h"
#include "PoissonContact.h"
#include "Dirichlet.h"
#include "PoissonModel.h"
#include "Neumann.h"


// Module interface
//TIBER_MODULE(Poisson,poisson)


using namespace std;


Poisson* Poisson::static_this;
Device* Poisson::_device;
//-----------------------------------------------------------------//


void Poisson::parse_options( )
{ 

 
}


void  Poisson::do_init( ) 
{
    const ModelOptions& sim_opt = get_options();

  SimulationEnvironment& si = get_environment();   

  _device = &( si.get_device() );

  mesh = & (_device->get_mesh());

  dim = mesh->mesh_dimension();
  
  double mesh_units = 100 * get_scaling().get_calc_mesh_units();
  get_scaling().set_calc_mesh_units(mesh_units);

  my_system = TiberLinearSystem::create(get_equation_systems(),
      get_equation_system_name(), get_solver_options());


  my_system->add_variable("V", FIRST);


   // Insert the pointer to function that LibMesh library has to use    
  my_system->attach_assemble_function (assemble_poisson_matrix);

   // Initialize the data structures for the equation system.
  my_system->init();	


}
//-------------------------------------------------------------------------------//
void  Poisson::do_solve()
{
  
  parse_options();
  
  static_this = this;
  
  my_system->set_options(get_solver_options());
  my_system->solve();
 
  check_gauss();

}


//--------------------------------------------------------------------------------//
Poisson::~Poisson()
{

  //equation_systems->delete_system(system_name);

}
//---------------------------------------------------------------------------------//
Poisson::Poisson()
{
  

}
//----------------------------------------------------------------------------------//

PhysicalModel*   Poisson::create_physical_model (const ModelOptions &options,
    const Material* mat) const throw (ModelErrorException)
{
  
  PoissonModel* model = dynamic_cast<PoissonModel*> ( PhysicalModelInterface::create("poisson",options) );

  if (model == NULL) 
    throw ModelErrorException("Poisson: Physical model is not created" );

  return model;      

}
//----------------------------------------------------------------------------------//

BoundaryProperties* Poisson::create_boundary_model (const ModelOptions &options) const 
                    throw (ModelErrorException)

{

  
   const string& modelname = options.get_option("type", "Dirichlet");
  
   
   PoissonContact* model =  PoissonContact::create(modelname, options);

   if (model == NULL)  
     throw ModelErrorException("Poisson: No such boundary model: " + modelname);

  return model;

}




Poisson*  Poisson::create(void)
{
  return new Poisson;
}






//----------------------------------------------------------------------------------//
void Poisson::build_nodal_results (const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				     std::vector< std::string > &legend)
{

  if (variables.find("ElPotential") != variables.end())
  {
    legend.resize(1);
    legend[0] = "Potential[V]";

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
      

      }
    }
  }
}




void
Poisson::build_elemental_results(const std::set<std::string>& variables,
					  std::vector<double>& results, std::vector<std::string>& legend)
{

  // we only do something if we are on processor 0
  // TODO parallelize
  if (libMesh::processor_id() != 0)
    return;

  const set<string>::const_iterator varend(variables.end());
  
  vector<ID> ids;
 
  unsigned int n_vars = 0;  

  legend.resize(variables.size());
  const unsigned int nn  = mesh->n_active_elem();
  const unsigned int dim = mesh->mesh_dimension();
 
  int PC = -1; 
  if (variables.find("PoissonCharge") != varend)
  {
    PC = n_vars;
    legend.resize(legend.size()+1);
    legend[n_vars] = "PoissonCharge[Q/cm3]";
    n_vars +=1;
  }

  int TP = -1; 
                     
  if (variables.find("Polarization") != varend)
  {
    TP = n_vars;
  
    legend.resize(legend.size() + dim);
    switch (dim)
    {
    case 3:
      legend[TP + 2] = "P_z";
      n_vars++;
    case 2:
      legend[TP + 1] =  "P_y";
      n_vars++;
      legend[TP + dim] = "modP";
      n_vars++;
    default:
      legend[TP] = "P_x";
      n_vars++;
    }
  }
  

  //  int DC = -1; 
  //if (variables.find("DielectricConstant") != varend)
  // {
  //  DC = n_vars;
  //  legend.resize(legend.size()+2);
  //  legend[n_vars] = "DielectricConstant[F/cm]";
  //  n_vars +=1;
  // }
  legend.resize(n_vars);
  results.resize(nn * n_vars,0.0);
  

  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();

   std::vector<double> charge_density;
   Tensor2Sym epsilon(0);
   RealVectorValue _tot_pol(3,0.0);

  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  { 
    const Elem* elem = *it;

    unsigned int id = n_vars * elem_number;
   
    std::vector<Point> _node(1);

    _node[0]=(elem->centroid());

    init_poisson_model(elem);
    

    if (PC != -1)
    {
         poisson_model->get_charge_density(_node,charge_density); 
         results[id + PC] = charge_density[0];
    }

    
    if (TP != -1)
    {
      _tot_pol = poisson_model->get_total_polarization();
      
      switch (dim)
      {
      case 3:
	results[id + TP + 2] = _tot_pol(2);
      case 2:
	results[id + TP + 1] =  _tot_pol(1);
	results[id + TP + dim] = sqrt(_tot_pol(0)*_tot_pol(0)+_tot_pol(1)*_tot_pol(1)+_tot_pol(2)*_tot_pol(2));
      default:
	results[id + TP ] = _tot_pol(0);
      }
     
    }


    // if (DC != -1)
    // {
    // epsilon = poisson_model->get_dielectric_constant();   
    //  results[id + DC] = epsilon(1,1);
  
    // }

    elem_number++;
  } //over element

  results.resize(elem_number * n_vars);
}




double Poisson::get_potential_element(const Elem* elem) const
{
        

       std::vector<unsigned int> dof_indices;

       DofMap& dof_map = my_system->get_dof_map();

       dof_map.dof_indices (elem, dof_indices); 

       double V = 0.0;
       for (unsigned int n = 0; n < elem->n_nodes(); n++)
       {   
	  V +=  (*(my_system->solution))(dof_indices[n]); 
       }
      
       V /= elem->n_nodes();

       return V;

}



std::vector<double>  Poisson::get_potential_node(const Elem* elem)
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
void Poisson::assemble_poisson_matrix(EquationSystems& es,
				     const std::string& system_name)
{



   static_this->do_assemble( es, system_name);

}

//----------------------------------------------------------------------------------//
void Poisson::do_assemble(EquationSystems& es, const std::string& system_name)
{


   //Fem Initialization-------------------------------
  SimulationEnvironment& se = get_environment();
  //Commons
  LinearImplicitSystem& system = *my_system;
  DofMap& dof_map =  system.get_dof_map();
  const unsigned int uvar = system.variable_number("V");
  FEType fe_type = dof_map.variable_type(uvar);
  std::vector<unsigned int> dof_indices;
  //Volume function 
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true));
  QGauss qrule (dim, FIFTH); 
  fe -> attach_quadrature_rule (&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();


  //Surface function
  AutoPtr<FEBase>  fe_face(build_finite_element(dim,fe_type,true));
  QGauss qface(dim-1, SIXTH);
  fe_face->attach_quadrature_rule(&qface);
  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& qface_point = fe_face->get_xyz();
  const std::vector<Point>& normal = fe_face->get_normals();
  const std::vector<std::vector<RealGradient> >& dphi_face = fe->get_dphi();
  //--------------------------  

  
  DenseMatrix<Number>  Ke;
  DenseVector<Number>  Fe;

  bool node_on_boundary = false;
  bool side_on_boundary = false;
  bool is_dirichlet = false;
  PoissonContact* contact;

  //Model Variables
  std::vector<double> charge_density;

  double total_charge = 0.0;
  Tensor2Sym epsilon(0);

  RealVectorValue _tot_pol(3,0.0);

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  for ( ; el != end_el ; ++el)   //loop over elements
  {  
  
    //Inizialize the element environment
    const Elem* elem = *el;

    dof_map.dof_indices (elem, dof_indices);  

    const unsigned int n_dofs   = dof_indices.size();
    
    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);

    Fe.resize(n_dofs);
    
    Fe.zero();  
     
    Ke.zero();

    //Init e read from heat model---------------------

    init_poisson_model(elem);

    poisson_model->get_charge_density(q_point,charge_density); 

    epsilon = poisson_model->get_dielectric_constant();

    _tot_pol = poisson_model->get_total_polarization();

    for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
    { // loop over test function
      //-------------------------------------------------
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)  
      {//Loop over quadrature points 
	
	for (unsigned int p2=0; p2<n_dofs; p2++) 
	{//loop over basis functions
	  
	  double value = 0.0;    
	  
	  for (short i = 0; i < dim; i++) 
	  {//loop over direction (1); test function derivative
	    
	    for (short j = 0; j < dim; j++)
	    {//loop over direction (2); basis function derivative
	      
             // double value =  JxW[qp] * dphi[p1][qp] * (epsilon * dphi[p2][qp]);
	      double epsilon_value;
	      if (i < j) 
		epsilon_value = epsilon(j+1, i+1);
	      else
		epsilon_value = epsilon(i+1, j+1);
	      
	      value += JxW[qp] * epsilon_value * dphi[p1][qp](i) * dphi[p2][qp](j);
	    }//end loop over direction (2)       
	  }//end loop over direction (1)												
	  
	  Ke(p1,p2) += value;
	  
	} //loop over basis functions
	
	//Fe construction----	
	
	Fe(p1) += JxW[qp] * charge_density[qp] * phi[p1][qp]; 
	
	total_charge +=  JxW[qp] * charge_density[qp];
	
	//Add the polarization 
	
	// for (short i = 0; i < dim; i++) 
	Fe(p1) += JxW[qp] * (dphi[p1][qp] * _tot_pol) ;
	//--------------------
      }//end Loop over quadrature points  
    }   //test function
    
    
    const unsigned int num_sides = elem->n_sides();  
    
    //Surface Heat Source
    
    for (unsigned int side = 0; side<num_sides; side++) 
    {   
      const ElementSide elside(elem->top_parent(), side); 
      
      if (se.is_on_boundary(elside))
      {
	fe_face->reinit(elem,side);  
	for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
	  for (unsigned int qp = 0; qp < qface.n_points(); qp++) 
	    Fe(p1) -=   JxW_face[qp] * (_tot_pol * normal[qp]) * phi_face[p1][qp];
	
      }
      
    }
    
    //Boundary condition
    
    for (unsigned int side = 0; side<num_sides; side++)   
    {
      const ElementSide elside(elem->top_parent(), side); 
      
      Boundary* bd = se.get_boundary(elside);
      
      if (bd != NULL)
      { 
	if (bd->get_boundary_properties( get_id() ) != NULL ) 
	{
	  contact = dynamic_cast<PoissonContact*>( bd->get_boundary_properties (get_id()) );    
	  switch (contact->get_type())
	  { 
	  case  PoissonContact::Dirichlet:
	    
	    for(unsigned int n = 0; n< n_dofs; ++n)             
	    {
	      
	      if (elem->is_node_on_side(n,side))
	      {
		for (unsigned int nc = 0; nc < n_dofs; nc++)    
		  Ke(n,nc) = 0.0;
		
		Ke(n,n) = 1.0;
		Fe(n) = (dynamic_cast<Dirichlet*> (contact) )->get_potential();
		
	      }
	    }
	    break;
	  }//case
	}//if (bd->get_boundary_properties( get_id() ) != NULL ) 
      }//if (bd != NULL)
    }//side
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices); 
    
    
  } //End Loop over elements
    //std::cout<<total_charge<<std::endl;
  
} //do assembly




void  Poisson:: init_poisson_model(const Elem* elem)
{     

       ID subdomain = elem->subdomain_id();
    
       const Material* mat = _device->get_material(subdomain);
    
       poisson_model =  (  dynamic_cast<PoissonModel*> (  mat ->get_model(get_id()) )  );
   
       poisson_model->set_element(elem);     

       //Update model for a given element
       poisson_model->re_init();
     
}




void
Poisson::get_solution(const Elem* elem, const std::vector<Point>& p,
		      std::vector<double>& solution)
{
  unsigned int np = p.size();

  solution.resize(np);

  if (np == 0) return;


  const DofMap& dof_map = my_system->get_dof_map();

  const unsigned int u_var = my_system->variable_number("V");

  FEType fe_type = my_system->variable_type(u_var);

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  std::vector<unsigned int> dof_indices;


  // element shape functions

  const vector<vector<Real> >& phi = fe->get_phi();

  vector<Point> points(np);

  FEInterface::inverse_map(dim, fe_type, elem, p, points);
 
  fe->reinit(elem, &points);

  dof_map.dof_indices(elem, dof_indices, u_var);

  const unsigned int n_dofs = dof_indices.size();

  for (unsigned int n = 0; n < np; n++)
  {
      // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
      solution[n] += phi[i][n]  *  (*(my_system->solution))(dof_indices[i]); ;

   
  }
}





ID
Poisson::convert_variable_name_to_id(const string& variable_name) const
{

   ID id = INVALID_ID;

    if (variable_name == "potential" )
      id  = POTENTIAL;
    if (variable_name == "e_field_x" )
      id  = EX; 
    if (variable_name == "e_field_y" )
      id  = EY; 
    if (variable_name == "e_field_z" )
      id  = EZ; 

  return id;
}




void
Poisson::get_solution_secure(const Elem* elem,
    const set<ID>& ids, vector<map<ID, double> >& values)
{

  LinearImplicitSystem& system = *my_system;

  DofMap& dof_map =  system.get_dof_map();

  vector<unsigned int> dof_indices;

  dof_map.dof_indices (elem, dof_indices);  

  for (unsigned int n = 0; n < elem->n_nodes(); n++)
  {
    
    if (ids.count(POTENTIAL))
      values[n][POTENTIAL] = (*(system.solution))(dof_indices[n]);

  }
   

}

void
Poisson::get_solution_secure(const Elem* elem, const vector<Point>& p,
    const set<ID>& ids, vector<map<ID, double> >& values)
{

  unsigned int np = p.size();
  values.resize(np);
  if ((np == 0) || (ids.size() == 0)) return;

  LinearImplicitSystem& system = *my_system;

  DofMap& dof_map =  system.get_dof_map();

  const unsigned int uvar = system.variable_number("V");

  FEType fe_type = dof_map.variable_type(uvar);
   AutoPtr<FEBase>  fe(build_finite_element(dim,fe_type,true));


  // element shape functions
   const vector<vector<Real> >& phi = fe->get_phi();
   const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  

  vector<Point> points(np);

  FEInterface::inverse_map(dim, fe_type, elem, p, points);
 
  fe->reinit(elem, &points);

  vector<unsigned int> dof_indices;

  dof_map.dof_indices (elem, dof_indices);  

  const unsigned int n_dofs   = dof_indices.size();

  std::vector<double> E(3);
  E.clear();
  double V = 0.0;
  init_poisson_model(elem); 


  for (unsigned int n = 0; n < np; n++)
  {
    V = 0.0;
    E.clear();

    //do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      V  += phi[i][n] * (*(system.solution))(dof_indices[i]);
      
      for (unsigned j = 0; j<dim; ++j)
      {
	E[j] -=  dphi[i][n](j) * (*(system.solution))(dof_indices[i]);
      }
  
    }
    if (ids.count(POTENTIAL))
      values[n][POTENTIAL] = V;
    
    if (ids.count(EX))
      values[n][EX] = E[0];
    
    if (ids.count(EY))
      values[n][EY] = E[1];
    
    if (ids.count(EZ))
      values[n][EZ] = E[2];
    
  }
  
}

double
Poisson::check_gauss()
{
  //Fem Initialization-------------------------------
  SimulationEnvironment& se = get_environment();
  //Commons
  LinearImplicitSystem& system = *my_system;
  DofMap& dof_map =  system.get_dof_map();
  const unsigned int uvar = system.variable_number("V");
  FEType fe_type = dof_map.variable_type(uvar);
  std::vector<unsigned int> dof_indices;
  //Volume function 
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type,true));
  QGauss qrule (dim, FIFTH); 
  fe -> attach_quadrature_rule (&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();


  //Surface function
  AutoPtr<FEBase>  fe_face(build_finite_element(dim,fe_type,true));
  QGauss qface(dim-1, SIXTH);
  fe_face->attach_quadrature_rule(&qface);
  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& qface_point = fe_face->get_xyz();
  const std::vector<Point>& normal = fe_face->get_normals();
  const std::vector<std::vector<RealGradient> >& dphi_face = fe->get_dphi();

  //Start calculation

  double total_charge = 0.0;
  double total_flux = 0.0;
  std::vector<double> charge_density;
  Tensor2Sym epsilon(0);
  Tensor1 normals;

 
  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  for ( ; el != end_el ; ++el)   //loop over elements
  {  
    //Initialize for a given element
    const Elem* elem = *el;
    dof_map.dof_indices (elem, dof_indices);  
    const unsigned int n_dofs = dof_indices.size();
    fe->reinit(elem);
    init_poisson_model(elem); 
   

    //Charge Density
    poisson_model->get_charge_density(q_point,charge_density); 
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)  
      total_charge +=  JxW[qp] * charge_density[qp];
   
   

    //Flux
    epsilon = poisson_model->get_dielectric_constant();
    const unsigned int num_sides = elem->n_sides();
    for (unsigned int side = 0; side<num_sides; side++)
    {
      ElementSide elside(elem->top_parent(),side);
      if (se.is_on_boundary(elside))
      {
        fe_face->reinit(elem,side);

        
	std::vector< std::map< ID, double > > electric_field;

	std::set<ID> IDFIELD;
        IDFIELD.insert(EX);
        IDFIELD.insert(EY);
        IDFIELD.insert(EZ);
       

        get_solution_secure(elem,qface_point,IDFIELD,electric_field); 

       

 	for (unsigned int qp = 0; qp <  qface.n_points(); qp++)	
 	{

	   Tensor1 E(0);
           Tensor1 D(0);
	   double flux = 0.0;

 	   E(1) = electric_field[qp].find(EX)->second;
 	   E(2) = electric_field[qp].find(EY)->second;
 	   E(3) = electric_field[qp].find(EZ)->second;
	 
          
           D = epsilon * E;

	   for( unsigned int kt =0;kt<3;kt++)   
	     normals(kt+1)=normal[qp](kt);

	 
           flux = D*normals;
	  
	
	  total_flux =  total_flux + flux;
	  
	   if (dim> 1)
 	    total_flux *= JxW[qp]; 
 
	  
 	 }

      }

    }

 
  }


  std::cout<< "Total Charge:"<< total_charge<<std::endl;
  std::cout<<"Total flux:"<< total_flux<<std::endl;
}
