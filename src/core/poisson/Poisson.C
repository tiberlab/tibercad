#include "Poisson.h"
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


  double mesh_units = _device->get_mesh_units();

  opt.work_units = sim_opt.get_option("Work_length_units", 1e-2);

  opt.length_scale = mesh_units/opt.work_units;

  equation_systems = &(get_equation_systems());

  system_name = get_equation_system_name();

  equation_systems->add_system<LinearImplicitSystem> (system_name);

  my_system = &( equation_systems->get_system<LinearImplicitSystem>(system_name)  );

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
  
  my_system->solve();

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

 

  if (variables.find("V") != variables.end())
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
  
  SimulationEnvironment& se = get_environment(); 

  LinearImplicitSystem& system = *my_system;
  
  const unsigned int uvar = system.variable_number("V");
  
  DofMap& dof_map = system.get_dof_map();
  
  FEType fe_type = dof_map.variable_type(uvar);
  
  // Declare a special finite element object for
  // volume integration.

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  //AutoPtr<FEBase> fe (build_finite_element(dim, fe_type));
  QGauss qrule (dim, FIFTH); //may be could be decreased (CHECK!!!)
  
  // quadrature rule  
  fe->attach_quadrature_rule (&qrule);
  
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

 // The element Jacobian * quadrature weight at each integration point.   
   const std::vector<Real>& JxW_face = fe_face->get_JxW(); 


   //The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
   const std::vector<Point>& qface_point = fe_face->get_xyz();  


  // The element shape functions evaluated at the quadrature point
   const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi(); 
		 
    // The normal to side
   const std::vector<Point>& normal = fe_face->get_normals(); 

 
  std::vector<unsigned int> dof_indices;
  
  DenseMatrix<Number>  Ke;

  DenseVector<Number>  Fe;

  bool node_on_boundary = false;
  bool side_on_boundary = false;
  bool is_dirichlet = false;
  PoissonContact* contact;

  //Model Variables
  std::vector<double> charge_density;

  Tensor2Sym epsilon(0);

  Tensor1 bi_pol(0);

  //-----------------------------------------------------------------//
  //My Jacobian. It is to pass to our work units
  
  double my_Jacobian = 1.0;
  for (short i = 1; i <= dim; i++)  my_Jacobian *= opt.length_scale;
  //----------------------------------------------------------------//


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
  
    // charge_density = poisson_model->get_charge_density();

    poisson_model->get_charge_density(q_point,charge_density);

    epsilon = poisson_model->get_dielectric_constant();

    bi_pol = poisson_model->get_built_in_polarization();

    

    //-------------------------------------------------

      
    for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
    { // loop over test function

      //!let us check if it belongs to a reservoir boundary
      const Node* nd = elem->get_node(p1);
	
      Boundary* bd = se.get_boundary(nd); 
       
      node_on_boundary =  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  );	
      
      if ( node_on_boundary )
      { 
          contact = dynamic_cast<PoissonContact*>( bd->get_boundary_properties (get_id()) );
          is_dirichlet = (contact->get_type() == PoissonContact::Dirichlet);

	  //          cout<<(contact->get_type() == PoissonContact::Dirichlet)<<endl;
	  //  cout<<(contact->get_type() == PoissonContact::Dirichlet)<<endl;
 
      }
       
      if ( is_dirichlet )
      { 
	  
 	  Ke(p1,p1) = 1.0;
	  
	  Fe(p1) = ( dynamic_cast<Dirichlet*> (contact) )->get_potential();

          is_dirichlet = false;
      } 
      else 
      {
       
	for (unsigned int qp=0; qp<qrule.n_points(); qp++)  
	{//Loop over quadrature points 

	  //------------------------
	  //volume integration for matrix
	  // \int_V  \frac{\patial \phi_{\alpha}{\partial x_i} \epsilon_{ij} \frac{\partial \phi_{\beta}{\partial x_j} dx
	  	  
	   for (unsigned int p2=0; p2<n_dofs; p2++) 

	  {//loop over basis functions
            

	    double value = 0.0;
	    
	   
	    for (short i = 0; i < dim; i++) 
	    {//loop over direction (1); test function derivative

	      for (short j = 0; j < dim; j++)
	      {//loop over direction (2); basis function derivative

		double epsilon_value;

		if (i < j) 
		  epsilon_value = epsilon(j+1, i+1);
		else
		  epsilon_value = epsilon(i+1, j+1);
		
	   	  value += JxW[qp] * epsilon_value * dphi[p1][qp](i) * dphi[p2][ qp](j) /(opt.length_scale * opt.length_scale);
                 
		
	      }//end loop over direction (2)

	       
	    }//end loop over direction (1)												
	    
	   value *= my_Jacobian;
	    

           Ke(p1,p2) += value;

	    
	  } //loop over basis functions
	   
	   //Fe construction----	
	
	   Fe(p1) += JxW[qp] * charge_density[qp] * phi[p1][qp] *  my_Jacobian  * Constants::e / Constants::epsilon;

           //Add the polarization 

	      for (short i = 0; i < dim; i++) 
		Fe(p1) -= JxW[qp] * dphi[p1][qp](i) * bi_pol(i+1) * my_Jacobian * Constants::e / Constants::epsilon;
	  

	   //--------------------
	   
	}//end Loop over quadrature points  
	
 	const unsigned int num_sides = elem->n_sides();  
	
 	for (unsigned int side = 0; side<num_sides; side++) 
 	{
	  
 	  const ElementSide elside(elem->top_parent(), side); 
	  
 	  Boundary* bd =   se.get_boundary(elside);
	  
          side_on_boundary =  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL ) );  
	  
 	  if (side_on_boundary)  
	  { 
  	    contact =  dynamic_cast<PoissonContact*>( bd->get_boundary_properties (get_id()) ); 
	    
  	    if ( (contact->get_type() == PoissonContact::Neumann)) 
 	    { 
	      
 	      double  D_condition  = ( dynamic_cast<Neumann*> (contact) )->get_polarization();  

	       D_condition /=epsilon(1,1);
 
 	      fe_face->reinit(elem, side);  
	      
	      for (unsigned int qp = 0; qp < qface.n_points(); qp++) 
 	      { 
			      
		//Fe(p1) =  JxW_face[qp] *  D_condition  * phi_face[p1][qp] *  my_Jacobian; 
 
		
 	      }  // for (unsigned int qp = 0; qp < qface.n_points(); qp++)    
	      
	      
	      
 	    }//if Neumann 
	    
	  }//if it is boundary with associated model 
	  


	  
	   //Add piezopolarization surface 
          
           bd =   poisson_model->get_piezo_environment().get_boundary(elside); 
	  
           if (bd != NULL) 
           { 
	    fe_face->reinit(elem, side);
	    if (dim>1)
            {
	      
	      for (unsigned int qp = 0; qp < qface.n_points(); qp++)
	      {
		for (short i = 0; i < dim; i++) 
		{
		  Fe(p1) += JxW[qp] * phi[p1][qp] * bi_pol(i+1) * my_Jacobian * normal[qp](i) *  Constants::e/Constants::epsilon ; 
		}
	      }
	    }    
	    else
	    {
	      double x_c = elem->centroid()(0);
	      double x_s = elem->point(p1)(0);
	      double n_1d = 1.0;
	      
	      if (x_s < x_c)
	      { 
		n_1d = -1.0;
	      }  
	      for (unsigned int qp = 0; qp < qface.n_points(); qp++)
	      {
		Fe(p1) += JxW[qp] * phi[p1][qp] * bi_pol(0) * my_Jacobian * Constants::e * n_1d * Constants::e/Constants::epsilon; 
	      }
	      
	    }//i dim >1   
	   }// if (bd != NULL) 
 	}//side 
	
	//----------------------------------------
	
      } //end if it is not a dirichlet boundary
      
    } // end loop over test functions
    
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices); 

    
  } //End Loop over elements
  
  
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

  const unsigned int uvar = system.variable_number("T");

  FEType fe_type = dof_map.variable_type(uvar);
  
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  // element shape functions
   const vector<vector<Real> >& phi = fe->get_phi();

  vector<Point> points(np);

  FEInterface::inverse_map(dim, fe_type, elem, p, points);
 
  fe->reinit(elem, &points);

  vector<unsigned int> dof_indices;

  dof_map.dof_indices (elem, dof_indices);  

  const unsigned int n_dofs   = dof_indices.size();

  for (unsigned int n = 0; n < np; n++)
  {
     double V = 0;

    //do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
   
      V  += phi[i][n] * (*(system.solution))(dof_indices[i]);
    }
    
     if (ids.count(POTENTIAL))
      values[n][POTENTIAL] = V;
  }
   
}
