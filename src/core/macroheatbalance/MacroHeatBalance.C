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
#include "HeatModel.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "Boundary.h"
#include "Reservoir.h"


using namespace std;
MacroHeatBalance* MacroHeatBalance::static_this;
Device* MacroHeatBalance::_device;
//-----------------------------------------------------------------//


void MacroHeatBalance::parse_options( )
{ 

  const ModelOptions& sim_opt = get_options();
 
  opt.current_simulation = sim_opt.get_option("current_simulation", "no_current");

  if (opt.current_simulation != "no_current")
  {
    _dd_simul = dynamic_cast< DriftDiffusion* > ( find_simulation(opt.current_simulation ) );
 
    if (_dd_simul == NULL)
      throw InitFailedException("Unknown current model" + opt.current_simulation );
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
   
  //------init is done---------------------------------------------------------------------//

}
//-------------------------------------------------------------------------------//
void  MacroHeatBalance::do_solve()
{

   parse_options();

   static_this = this;

   my_system->solution->zero();

   my_system->solve();

}
//--------------------------------------------------------------------------------//
MacroHeatBalance::~MacroHeatBalance()
{
  equation_systems->delete_system(system_name);

}
//---------------------------------------------------------------------------------//
MacroHeatBalance::MacroHeatBalance()
{
  _dd_simul = NULL;

}
//----------------------------------------------------------------------------------//
PhysicalModel*   MacroHeatBalance::create_physical_model (const ModelOptions &options) const 
                    throw (ModelErrorException)
{
  
  HeatModel* model = HeatModel::create( );

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
    legend[0] = "Temperature [K]";



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

 QGauss qrule (dim, FIFTH); //may be could be decreased (CHECK!!!)

  // quadrature rule  
  fe -> attach_quadrature_rule (&qrule);
 

  // Declare a special finite element object for
  // boundary integration.
  AutoPtr<FEBase>  fe_face(FEBase::build(dim, fe_type));

 // Boundary integration requires one quadraure rule,
 // with dimensionality one less than the dimensionality
 // of the element.
 QGauss qface(dim-1, THIRD);
  
 // Tell the finite element object to use our
 // quadrature rule.
 
    fe_face -> attach_quadrature_rule (&qface);

 // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
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

  Tensor2Sym kappa;

  HeatModel* heat_model;

  const LatticeThermalConductivity* lattice_conductivity; 
  

  //-----------------------------------------------------------------//
  //My Jacobian. It is to pass to our work units
  
  double my_Jacobian = 1.0;
  for (short i = 1; i <= dim; i++)  my_Jacobian *= opt.length_scale;
  //----------------------------------------------------------------//

  
  for ( ; el != end_el ; ++el) 
    {  
    const Elem* elem = *el;

    ID subdomain = elem->subdomain_id();
    const Material* mat = _device->get_material(subdomain);

    heat_model =  (  dynamic_cast<HeatModel*> (  mat ->get_model(get_id()) )  ); 

    lattice_conductivity = heat_model->get_lattice_conductivity();

    lattice_conductivity->get_conductivity(kappa);

    

    dof_map.dof_indices (elem, dof_indices); 
    
    const unsigned int n_dofs   = dof_indices.size();
    
    fe->reinit (elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {//volumic integration loop
      for (unsigned int p1=0; p1<n_dofs; p1++) //basis functions of the variable T
      {
	//!let us check if it belongs to a reservoir bounday
	const Node* nd = elem->get_node(p1);
	
	Boundary* bd = se.get_boundary(nd); 
	
	if (  (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  ) )
	{
	  ThermalContact* contact = dynamic_cast<ThermalContact*>( bd->get_boundary_properties (get_id()) );
	  if (contact->get_type() == 0)
	  {
	   
	   
	    Ke(p1,p1) = 1.0;
	    
	    Fe(p1) = ( dynamic_cast<Reservoir*> (contact) )->get_temperature();
	  }
	}
	else
	{
	  for (unsigned int p2=0; p2<n_dofs; p2++) //test functions
	  {
	    double value = 0.0;
	    
	    for (short i = 0; i < dim; i++) //basis function derivative
	      for (short j = 0; j < dim; j++) //test function derivative
	      {
		double kappa_value;
		if (i > j) 
		  kappa_value = kappa(j+1, i+1);
		else
		  kappa_value = kappa(i+1, j+1);
		
		value += -JxW[qp] * kappa_value * dphi[p1][qp](i) * dphi[p2][qp](j) /(opt.length_scale * opt.length_scale);
		
	      }
	    
	    value *= my_Jacobian;
	   
	    Ke(p1,p2) += value;
	  }
	}
      } 
      
	
    }
    
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

    
    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices); 
    
    
    

  }

  // system.matrix->print_matlab("Matr.m");
 
}
