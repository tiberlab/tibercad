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

  equation_systems->delete_system(system_name);

}
//---------------------------------------------------------------------------------//
Poisson::Poisson()
{
  

}
//----------------------------------------------------------------------------------//

PhysicalModel*   Poisson::create_physical_model (const ModelOptions &options) const 
                    throw (ModelErrorException)
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
   
   const string& modelname = options.get_option("type", "Dirichelet");
 
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
  
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  //AutoPtr<FEBase> fe (build_finite_element(dim, fe_type));
  
  QGauss qrule (dim, FIFTH); //may be could be decreased (CHECK!!!)
  
  // quadrature rule  
  fe->attach_quadrature_rule (&qrule);
  

   // The element Jacobian * quadrature weight at each integration point.   
    const std::vector<Real>& JxW = fe->get_JxW();


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


  bool is_on_bounary = false;

  PoissonContact* contact;

  //Model Variables
   double charge_density;

   Tensor2Sym epsilon;

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
    
    fe->reinit(elem);

   

    Ke.resize(n_dofs, n_dofs);

    Fe.resize(n_dofs);
    
    Fe.zero();  
     
    Ke.zero();


    //Init e read from heat model---------------------

    init_poisson_model(elem);
  
    charge_density = poisson_model->get_charge_density();

    epsilon = poisson_model->get_dielectric_constant();// *  opt.work_units;

    

    //-------------------------------------------------

      
    for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
      
    { // loop over test function

      //!let us check if it belongs to a reservoir boundary
      const Node* nd = elem->get_node(p1);
	
      Boundary* bd = se.get_boundary(nd); 
       
      is_on_bounary =  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  );	
      
      if ( is_on_bounary )
      { //if belongs to boundary
	
	contact = dynamic_cast<PoissonContact*>( bd->get_boundary_properties (get_id()) );
	
	if (contact->get_type() == PoissonContact::Dirichlet)
	{//heat reservoir---
	    
 	  Ke(p1,p1) = 1.0;
	  
	  Fe(p1) = ( dynamic_cast<Dirichlet*> (contact) )->get_potential();
          
	}//end reservoir

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
           
	   Fe(p1) += JxW[qp] * charge_density * phi[p1][qp] *  my_Jacobian  * Constants::e / Constants::epsilon ;
	   
	   //--------------------
   
	}//end Loop over quadrature points  
	

      } //end if it belongs to boundary

 
    } // end loop over test functions

    //Add Neumann boundary condition
    
/*     if ( is_on_bounary ) */
/*     { */

/*       if (contact->get_type() == PoissonContact::Neumann) */
/*       { */
	
/* 	const unsigned int num_sides = elem->n_sides(); */
	
/* 	for (unsigned int side = 0; side<num_sides; side++) */
/* 	{ */
/* 	  const ElementSide elside(elem->top_parent(), side); */
	  
/* 	  if ( (se.is_on_boundary(   elside   ) ) //if belongs to a boundary of poisson simulation */
/* 	    { 		    */
	      
/* 	      if (dim > 1) */
		
/* 	       { */
		 
/* 		 fe_face->reinit(elem, side); */
		 
/* 		 const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi(); */
		 
/* 		 const std::vector<Real>& JxW_face = fe_face->get_JxW(); */
		 
/* 		 const std::vector<Point>& qface_point = fe_face->get_xyz(); */
		 
/* 		 const std::vector<Point> & normal = fe_face->get_normals(); */

		 
/* 		 for (short i = 0; i < 3; i++) */
/* 		 { */
/* 		   for (unsigned int qp=0; qp < qface.n_points(); qp++) */
/* 		   { */
/* 		     Fe(p1) += (JxW_face[qp] * phi_face[p1][qp]) * normal[qp](i) *  */
/* 		       ( face_currents[qp].jn(i)*face_potentials[qp].fermi_e + face_potentials[qp].fermi_h * face_currents[qp].jp(i) ) *  */
/* 		       (my_Jacobian/opt.length_scale); */
/* 		   } */
/*                  } */


/* 	       } */
/* 	       else //dim = 1 */
/* 	       { */
		 
/* 		 if (p1== side) */
/* 		 { */
/* 		   std::vector<double> normal(3); */
/* 		   Point p = elem->point(side); */
/* 		   Point pc = elem->centroid(); */
		   
		   
		   
/* 		   const double temp = sqrt((p(0) - pc(0)) * (p(0) - pc(0))  */
/* 					    +(p(1) - pc(1)) * (p(1) - pc(1)) +  */
/* 					    (p(2) - pc(2)) * (p(2) - pc(2))); */
		   
		   
/* 		   for (short i = 0; i < 3; i++) */
/* 		     normal[i] = (p(i) - pc(i))/temp; */
		   
		   
/* 		   std::vector<Point> qface_point(1); */
		   
/* 		   qface_point[0] = elem->point(p1); */
		   
		   
	
/* 		   heat_model->get_dd_solution(qface_point,face_potentials,face_currents);  		    */
		   
/* 		   for (short i = 0; i < 3; i++) */
/* 		   {   */
/* 		     Fe(p1) +=   normal[i] *  */
/* 		       ( face_currents[0].jn(i) * face_potentials[0].fermi_e + face_potentials[0].fermi_h * face_currents[0].jp(i) )   ; */
		     
		     
/* 		   } */
		   
		   
		   
		   
/* 		 } //if (dim > 1) */
		 
/* 	       } */
/* 	     } // if ( (dd_simul->get_environment()).is_on_boundary(   elside   ) ) //if belongs to a boundary of current(!) simulation */
/* 	   } //  for (unsigned int side = 0; side<num_sides; side++) */
	   
/* 	 } // if ( !belongs_to_reservoir )  */
	 
/*        }// for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T */
       
/*      } // (dd_simul != NULL) */























     //The loop over element is the only loop that is surviving at this point
  
   /*   if (poisson_model->get_strain_opt()) */

/*      { */
/*        for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable V */
/*        { */
/* 	 //!let us check if it belongs to a dirichlet boundary */

/* 	 const Node* nd = elem->get_node(p1); */
	 
/* 	 Boundary* bd = se.get_boundary(nd);  */
	 
/* 	 bool dirichlet_contact = false; */
	 
/* 	 if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) ) */
/* 	 { */
/* 	   PoissonContact* contact = dynamic_cast<PoissonContact*>( bd->get_boundary_properties (get_id()) );     */
	   
/* 	   if (contact->get_type() == PoissonContact::Dirichlet) dirichlet_contact = true; */
/* 	 } */
	 
/* 	 if ( !dirichlet_contact )  */
/* 	 {//not fixed potential */
	   
/* 	   const unsigned int num_sides = elem->n_sides(); */
	   
/* 	   for (unsigned int side = 0; side<num_sides; side++) */
/* 	   { */
/* 	     const ElementSide elside(elem->top_parent(), side); */
	     
/* 	     if ( (poisson_model->get_strain_environment()).is_on_boundary(elside) ) //if belongs to a boundary of current(!) simulation */
/* 	     { */
	      
/* 	       std::vector<DriftDiffusion::Solution>  face_potentials;    */
	       
/* 	       std::vector<DriftDiffusion::Currents>   face_currents; */
		   

/* 	       if (dim > 1) */
/* 	       { */
		 
/* 		 fe_face->reinit(elem, side); */
		 
/* 		 const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi(); */
		 
/* 		 const std::vector<Real>& JxW_face = fe_face->get_JxW(); */
		 

/* 		 const std::vector<Point>& qface_point = fe_face->get_xyz(); */
		 
/* 		 const std::vector<Point> & normal = fe_face->get_normals(); */
		 
		    
/* 		 heat_model->get_dd_solution(qface_point,face_potentials,face_currents);   */
		     
		 
/* 		 for (short i = 0; i < 3; i++) */
/* 		 { */
/* 		   for (unsigned int qp=0; qp < qface.n_points(); qp++) */
/* 		   { */
/* 		     Fe(p1) += (JxW_face[qp] * phi_face[p1][qp]) * normal[qp](i) *  */
/* 		       ( face_currents[qp].jn(i)*face_potentials[qp].fermi_e + face_potentials[qp].fermi_h * face_currents[qp].jp(i) ) *  */
/* 		       (my_Jacobian/opt.length_scale); */
/* 		   } */
/*                  } */


/* 	       } */
/* 	       else //dim = 1 */
/* 	       { */
		 
/* 		 if (p1== side) */
/* 		 { */
/* 		   std::vector<double> normal(3); */
/* 		   Point p = elem->point(side); */
/* 		   Point pc = elem->centroid(); */
		   
		   
		   
/* 		   const double temp = sqrt((p(0) - pc(0)) * (p(0) - pc(0))  */
/* 					    +(p(1) - pc(1)) * (p(1) - pc(1)) +  */
/* 					    (p(2) - pc(2)) * (p(2) - pc(2))); */
		   
		   
/* 		   for (short i = 0; i < 3; i++) */
/* 		     normal[i] = (p(i) - pc(i))/temp; */
		   
		   
/* 		   std::vector<Point> qface_point(1); */
		   
/* 		   qface_point[0] = elem->point(p1); */
		   
		   
	
/* 		   heat_model->get_dd_solution(qface_point,face_potentials,face_currents);  		    */
		   
/* 		   for (short i = 0; i < 3; i++) */
/* 		   {   */
/* 		     Fe(p1) +=   normal[i] *  */
/* 		       ( face_currents[0].jn(i) * face_potentials[0].fermi_e + face_potentials[0].fermi_h * face_currents[0].jp(i) )   ; */
		     
		     
/* 		   } */
		   
		   
		   
		   
/* 		 } //if (dim > 1) */
		 
/* 	       } */
/* 	     } // if ( (dd_simul->get_environment()).is_on_boundary(   elside   ) ) //if belongs to a boundary of current(!) simulation */
/* 	   } //  for (unsigned int side = 0; side<num_sides; side++) */
	   
/* 	 } // if ( !belongs_to_reservoir )  */
	 
/*        }// for (unsigned int p1=0; p1<n_dofs; p1++) //test functions of the variable T */
       
/*      } // (dd_simul != NULL) */



     //Restore the boundary information for the next element
     is_on_bounary = false;
     //---------------------------------------------------




     dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
     system.matrix->add_matrix (Ke, dof_indices);
     system.rhs->add_vector    (Fe, dof_indices); 
      
      
  } //End Loop over elements
    

 
  //  system.matrix->print_matlab("Matr.m");
  //  system.rhs->print();
    
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
