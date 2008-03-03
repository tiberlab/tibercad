#include "Macrostrain.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "BoundaryProperties.h"
#include "Boundary.h"
#include "MacrostrainBoundaryProperties.h"
#include "MacrostrainPressure.h"
#include "MacrostrainSubstrate.h"
#include "TiberPetscLinearSolver.h"




using namespace std;
//-----------------------------------------------------------------//


Device*  Macrostrain:: _device;
Macrostrain* Macrostrain::static_this;

ID Macrostrain::convert_variable_name_to_id(const std::string& variable_name) const
{
  
  ID id = INVALID_ID;
 
  if (variable_name == "") return id;

  if (variable_name == "eps_xx")
    id = EPS_XX;
  else if  (variable_name == "eps_yy")
    id = EPS_YY;
  else if  (variable_name == "eps_zz")
    id = EPS_ZZ;
  else if  (variable_name == "eps_xy" || variable_name == "eps_yx")
    id = EPS_XY;
  else if  (variable_name == "eps_yz" || variable_name == "eps_zy")
    id = EPS_YZ;
  else if  (variable_name == "eps_xz" || variable_name == "eps_zx")
    id = EPS_XZ;
  else if  (variable_name == "Px")
    id = P_X;
  else if  (variable_name == "Py")
    id = P_Y;
  else if  (variable_name == "Pz")
    id = P_Z;

  return id;
  
}
 
//-----------------------------------------------------------------//

    

 
//-----------------------------------------------------------------//

void Macrostrain::get_solution_secure(const Elem* elem,
				   const std::vector<Point>& p, const std::set<ID>& ids,
				   std::vector<std::map<ID, double> >& values)
{
  unsigned int np = p.size();
 
 

  Tensor2Sym strain_el = result_strain[elem];

  ID subdomain = elem->subdomain_id();
          
  const Material* mat = _device->get_material(subdomain);
  
  const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

  MacrostrainModel* macrostrain_model =  dynamic_cast<MacrostrainModel*>(   mat ->get_model(get_id())     );
      
  Tensor1 polariz = ( macrostrain_model->get_piezo() )-> get_polariz_cryst(strain_el); //crystal system

  polariz = (crystal_el->RotMatrix) * polariz; //calculation system
 
  for (unsigned int n = 0; n < np; n++)
  {
    if (ids.count(EPS_XX))
      values[n][EPS_XX] = strain_el(1,1);

    if (ids.count(EPS_YY))
      values[n][EPS_YY] = strain_el(2,2);

    if (ids.count(EPS_ZZ))
      values[n][EPS_ZZ] = strain_el(3,3);


    if (ids.count(EPS_XY))
      values[n][EPS_XY] = strain_el(2,1);

    if (ids.count(EPS_YZ))
      values[n][EPS_YZ] = strain_el(3,2);

    if (ids.count(EPS_XZ))
      values[n][EPS_XZ] = strain_el(3,1);
    

    if (ids.count(P_X))
      values[n][P_X] = polariz(1);


    if (ids.count(P_Y))
      values[n][P_Y] = polariz(2);

    if (ids.count(P_Z))
      values[n][P_Z] = polariz(3);


  }
  

}





//-----------------------------------------------------------------//
void Macrostrain::build_elemental_results(const std::set<std::string>& variables,
			     std::vector<double>& results, std::vector<std::string>& legend)
{

  std::vector<std::string> eps_names;
  std::vector<double> eps_data;  

  std::vector<std::string> pol_names;
  std::vector<double> pol_data;  

  prepare_strain_data_for_output( eps_names,  eps_data);
  prepare_polarization_data_for_output( pol_names,  pol_data);


  
  short num_var = 0;
  const set<string>::const_iterator varend = variables.end();
  
  const string strain_name("strain");
  const string pol_name("polarization");

  if (variables.find(strain_name) != varend) num_var += 6;  
  if (variables.find(pol_name) != varend) num_var +=3;

  unsigned int num_elem = eps_data.size()/6;



  results.resize(num_var * num_elem);
  legend.resize(num_var);

  

  if (variables.find(strain_name) != varend)
  {//we do strain
    for (short i = 0; i < 6; i++)
    {	
      legend[i] = eps_names[i];
      for (unsigned int j = 0; j < num_elem; j++)
	results[i + j * num_var ] = eps_data[i + j * 6];  
    }
    
  }
  

  
  if (variables.find(pol_name) != varend)
  {//now we do polarization
    for (short i1 = 0; i1 < 3; i1++)
    {  
      
      legend[i1 + num_var - 3] = pol_names[i1];
      for (unsigned int j = 0; j < num_elem; j++)
	results[i1 + num_var - 3 + j * num_var ] = pol_data[i1 + j * 3];  
      
      
    }
  }

}

//-------------------------------------------------------------------------//

PhysicalModel*
Macrostrain::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{
  const string& modelname = options.get_option("model", "macrostrain");

  MacrostrainModelInterface* model =
    MacrostrainModelInterface::create(modelname, options);

  if (model == NULL)
    throw ModelErrorException(
        "Macrostrain: No such physical model: " + modelname);

  return model;
}


//---------------------------------------------------------------------------//


BoundaryProperties*
Macrostrain::create_boundary_model(const ModelOptions& options) const
throw (ModelErrorException)
{
  //const string& modelname = options.get_option("BC_region_name", "");
  //if (modelname != "substrate")

 
  
  const string& modelname = options.get_option("type", "pressure");

  MacrostrainBoundaryProperties* model =
    MacrostrainBoundaryProperties::create(modelname, options);

  if (model == NULL)
    throw ModelErrorException(
        "Macrostrain: No such boundary model: " + modelname);

  return model;
  


}


//----------------------------------------------------------------------------//


//-----------------------------------------------------------------//


//-----------------------------------------------------------------//
void Macrostrain::parse_options( )
{ 

 const ModelOptions& opt = get_options();

 


 max_r_steps = opt.get_option("refinement_steps", 0);

 uniform_refinement = opt.get_option("uniform_refinement", false);


 refine_fraction = opt.get_option("refine_fraction", 0.25);

 coarsen_fraction = opt.get_option("coarsen_fraction", 0.0);

 max_ref_level = opt.get_option("max_refinement_level",10);
 tolerance  = opt.get_option("tolerance", 1e-10);  
 max_shape_steps = opt.get_option("number_shape_steps",0);
   

 calculate_atom_displacements = opt.get_option("calculate_atom_displacements", false);
 atom_structure_filename = opt.get_option("atom_structure_filename", "");
 atom_displacements_filename = opt.get_option("atom_displacements_filename","");
 atom_potential_filename = opt.get_option("atom_potential_filename","");
  

  
 unsigned int max_ksp_iterations = opt.get_option("max_iterations",1000);
 

 // assert(periodicity_x == false);

   
 equation_systems->parameters.set<Real>("linear solver tolerance") = tolerance; 
 
  




 SolverType solver_type;

 PreconditionerType pc_type;

 string ksptype;

 string pc; 

 if (dim == 1)
 {
   ksptype = "gmres";

   pc =  "ilu" ;
 }
 else
 { 
   ksptype = "bcgsl";

   pc =  "jacobi" ;
 }

 ksptype = opt.get_option("ksp_type", ksptype );

 if (ksptype == "") {}
 else if (ksptype == "bcgsl")
  solver_type  = BICGSTAB;
 else if (ksptype == "gmres")
  solver_type = GMRES;
 else if (ksptype == "bcgs")
  solver_type = BICG;
 else if (ksptype == "cg")
  solver_type = CG;
 else if (ksptype == "richardson")
  solver_type = RICHARDSON;



 pc = opt.get_option("pc_type", pc);

 if (pc == "") {}
 else if (pc == "ilu")
   pc_type = ILU_PRECOND;
 else if (pc == "composite")
   pc_type = USER_PRECOND;
 else if (pc == "jacobi")
   pc_type = JACOBI_PRECOND;
 else if (pc == "lu")
   pc_type = LU_PRECOND;
 else if (pc == "cholesky")
   pc_type = CHOLESKY_PRECOND;
 else if (pc == "eisenstat")
   pc_type = EISENSTAT_PRECOND;
  


 


 my_system->linear_solver->set_solver_type(solver_type);

 my_system->linear_solver->set_preconditioner_type(pc_type);


 my_solver->set_ksp_options ( tolerance, max_ksp_iterations);

 my_solver->init();



 
   
 bool monitor = opt.get_option("monitor", false);

 if (monitor)
 {
     
   int ierr; 

   KSP ksp_of_my_solver = (dynamic_cast< TiberPetscLinearSolver* > (  (my_system->linear_solver).get() )   )->get_ksp();

#if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) \
      && (PETSC_VERSION_SUBMINOR >= 2))
   ierr = KSPMonitorSet(ksp_of_my_solver,KSPMonitorDefault, PETSC_NULL,0);
#else
   ierr = KSPSetMonitor(ksp_of_my_solver,KSPDefaultMonitor, PETSC_NULL,0);
#endif

  
 }

 bool xmonitor = opt.get_option("xmonitor", false);

 
 if (xmonitor)
 {
     
#if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) \
      && (PETSC_VERSION_SUBMINOR >= 2))
   KSPMonitorLGCreate (NULL, get_name().c_str(),0,0,400,400, &_lg);
#else
   KSPLGMonitorCreate (NULL, get_name().c_str(),0,0,400,400, &_lg);
#endif

   int ierr; 

   KSP ksp_of_my_solver = (dynamic_cast< TiberPetscLinearSolver* > (  (my_system->linear_solver).get() )   )->get_ksp(); 

#if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) \
      && (PETSC_VERSION_SUBMINOR >= 2))
   ierr = KSPMonitorSet ( ksp_of_my_solver,KSPMonitorLG,_lg,0);
#else
   ierr = KSPSetMonitor ( ksp_of_my_solver,KSPLGMonitor,_lg,0);
#endif
 }


/*

 {
   int ierr; 
 
   KSP KSP_of_solver = (dynamic_cast< PetscLinearSolver<Real>* > (  (my_system->linear_solver).get() )  )->ksp();
   ierr = KSPSetType (KSP_of_solver, (char*) KSPBCGSL);    

   PC  PC_of_solver = (dynamic_cast< PetscLinearSolver<Real>* > (  (my_system->linear_solver).get() )  )->pc();
   ierr = PCSetType (PC_of_solver, (char*) PCJACOBI);   

 }

*/



 if (!grown_on_substrate)
 {
   vector<double> point;
   if ( ! opt.find_option("fixed_point1"))
   {
     throw InitFailedException( "Macrostrain: fixed_point1 is not defined");
   }
     

   opt.get_option("fixed_point1",point);
   for (short i = 0; i < 3; i++)  fixed_point1(i) = point[i];
   


   if (dim>1)
   {

     if ( ! opt.find_option("fixed_point2"))
     {
       throw InitFailedException( "Macrostrain: fixed_point2 is not defined");
     }

     opt.get_option("fixed_point2",point);
     for (short i = 0; i < 3; i++)  fixed_point2(i) = point[i];
     if (dim>2)
     {
       if ( ! opt.find_option("fixed_point3"))
       {
	 throw InitFailedException( "Macrostrain: fixed_point2 is not defined");
       }

       opt.get_option("fixed_point3",point);
       for (short i = 0; i < 3; i++)  fixed_point3(i) = point[i];
     }
   }
 }
  


  //-------------------------------------------------------------------//

  


    intermediate_output = opt.get_option("intermediate_output",false);
    
    output_strain_on_atoms = false;
    atom_output_type = "uptight";

 
    output_type = opt.get_option("output_type","GMV");
 

    {
      //-----------------------------------------------------------------//
      //potential on atoms
      //should be removed
      
      std::string poisson_model_name = opt.get_option("poisson_model_name","no_poisson");
      if ( poisson_model_name != "no_poisson" )
      {

    
	poisson_equation  = find_simulation ( poisson_model_name );
	
	if (poisson_equation == NULL)
	  throw InitFailedException( "Unknown poisson model " + poisson_model_name);

      }
    }
   

}


//-----------------------------------------------------------------//
void Macrostrain::do_init( ) 
{

  const ModelOptions& options = get_options();

  
  equation_systems = & (get_equation_systems());


  SimulationEnvironment& si = get_environment();   

  _device = &( si.get_device() );


  get_scaling().set_length_scaling(  _device->get_mesh_units() );


  uname_vec[0]="ux";
  uname_vec[1]="uy";
  uname_vec[2]="uz";

  mesh = &(equation_systems->get_mesh());
      
  dim = mesh->mesh_dimension();



  substrate_name = options.get_option("substrate","_no_substrate_");

  Boundary* substrate_boundary= si.get_boundary ( substrate_name) ;

  grown_on_substrate = false;

  substrate_crystal = NULL;

  if (substrate_boundary != NULL)
  {
    BoundaryProperties*  bp =  substrate_boundary->get_boundary_properties (get_id() );

    if (bp !=NULL) 
    {
      grown_on_substrate = true;

      MacrostrainSubstrate* subst = dynamic_cast<MacrostrainSubstrate*>(bp);
      if (subst == NULL)
        throw InitFailedException("Macrostrain: Boundary model of \'" +
            substrate_name + "\' is not substrate model.");

      Material* mat = subst->get_material();
      if (mat == NULL)
        throw InitFailedException("Macrostrain: Substrate \'" + substrate_name +
            "\' has no material.");
      
      substrate_crystal =  &(mat->get_rotated_crystal());
    }


  }
  

  if (!grown_on_substrate)
  {
    Point p;
    vector<double>  ref_point(3);
    ref_point[0] = 0; ref_point[1] = 0; ref_point[2] = 0;
    
    if (!options.find_option("reference_material_point"))
    {
      throw InitFailedException("Macrostrain: reference material point must be given");
    }

    options.get_option("reference_material_point", ref_point);


    for (short i = 0; i < 3; i++)  p(i) = ref_point[i];
    
    MeshBase::const_element_iterator el  = mesh->active_elements_begin();
    MeshBase::const_element_iterator end_el = mesh->active_elements_end();
    
    const Elem* elem1 = NULL;

    for ( ; ( el != end_el ) ; ++el)  
    {
      Elem* elem = *el;
      if (   may_belong_to_element(elem,  p) )
      {
	if (elem->contains_point(p))
	{
	  if (elem->contains_point(p))
	  {
	    elem1 = elem;
	    
	    break;
	  }
	}
      }
	  
    }
      
    ID subdomain = elem1->subdomain_id();
    const Material* mat = _device->get_material(subdomain);
    substrate_crystal  = &(mat->get_rotated_crystal());
   
    
  }


  periodicity[0] = options.get_option("periodicity_x",false);
  periodicity[1] = options.get_option("periodicity_y",false);
  periodicity[2] = options.get_option("periodicity_z",false);
 


  define_additional_variables();
 

  


  //--------------------------------------------------------------------------------------//
  //add new system
  system_name = get_equation_system_name ( );
  
 
  equation_systems->add_system<LinearImplicitSystem> (system_name);

  my_system = &( equation_systems->get_system<LinearImplicitSystem>(system_name)  );

  //--------------------------------------------------------------------------------------//

  //--------------------------------------------------------------------------------------//
  //add normal variables
	
  for (unsigned int i = 0; i <  3 ; i++)  
  {  
    my_system->add_variable(uname_vec[i], FIRST);
  }
    
 
  //---------------------------------------------------------------------------------------//
  //add aditional varables 
  
  
  if (number_of_add_var != 0)
  {
    FEType fe_type(CONSTANT,MONOMIAL);
    my_system->add_variable("fict", fe_type);
  }
  
  //---------------------------------------------------------------------------------------//
  
  my_system->attach_assemble_function (assemble_strain_matrix);



  //------------------------------------------------------
  NumericVector<Number>& old_solution = 
    my_system->add_vector("old solution");



  
  // Initialize the data structures for the equation system.
  my_system->init();	



   //---------------------------------------------------------------------//
 //define  max and min coordinates 


 

    
 unsigned int num_nodes = mesh->n_nodes();
 const Node& nd = mesh->node(0);
 for (unsigned i = 0; i < 3; i++)
 {
   min_coord[i] = nd(i);
   max_coord[i] = nd(i);
 }

 for (unsigned i = 1; i < num_nodes; i++)
 {
   const Node& nd = mesh->node(i);
   for (unsigned i = 0; i < 3; i++)
   {
     if (min_coord[i] < nd(i)) min_coord[i] = nd(i);
     if (max_coord[i] > nd(i)) max_coord[i] = nd(i);
     
   }
   
 }
   
  //-------------------------------------------------------------------//
 my_solver =  TiberLinearSolver::create ("petsc");
 my_system->linear_solver = AutoPtr< LinearSolver< Real > >(my_solver);

  //---------------------------------------------------------------------------------------------------------//




  //------init is done---------------------------------------------------------------------//
}


//--------------------------------------------------------------------//

void Macrostrain::define_fixed_nodes()
{
  
  fixed_node1 = find_nearest_node(fixed_point1);
 
  if (dim>1)
  {
    fixed_node2 = find_nearest_node(fixed_point2);
    if (dim>2)
    {
      fixed_node3 = find_nearest_node(fixed_point3);
    }
  }

}

//--------------------------------------------------------------------//

void Macrostrain::assemble_strain_matrix(EquationSystems& es,
				     const std::string& system_name)
{

  static_this->do_assemble( es, system_name);

}



//-----------------------------------------------------------------//

void Macrostrain::do_assemble(EquationSystems& es,
				     const std::string& system_name)

{ //

 
 

 

  int temp_i;
  temp_i = 0;

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Matrix Assembly",false);

 
  
  // Get a constant reference to the mesh object.
  const Mesh& mesh = es.get_mesh();

  unsigned int dim = mesh.mesh_dimension();

  // The dimension that we are running
   
  //dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem& system = *my_system;

 
  
  unsigned int uvar[3] ;
  unsigned int var_fict;

  for (unsigned int i = 0; i<= 3 - 1; i++) 
  {
    uvar[i] = system.variable_number(uname_vec[i]);
  }
  

  if (number_of_add_var !=0 ) var_fict = system.variable_number("fict");

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  // const DofMap& dof_map = system.get_dof_map();
  DofMap& dof_map = system.get_dof_map();
 
  

  FEType fe_type = dof_map.variable_type(uvar[0]);
 
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type, true)); //I need scaling here 
 


  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);
  
  // Tell the finite element object to use our quadrature rule.

  fe -> attach_quadrature_rule (&qrule);
 

  // Declare a special finite element object for
  // boundary integration.
  AutoPtr<FEBase>  fe_face(build_finite_element(dim, fe_type, true)); //I need scaling here 


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
  
  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".  These datatypes are templated on
  //  Number, which allows the same code to work for real
  // or complex numbers.
  

  //-------------------------------------------------------------
  //matrixes to built the system
 
  
  DenseMatrix<Number> Ke_total; 
  DenseVector<Number> Fe_total;

  //------------------------------------------------------------

  //-------------------------------------------------------------
  //submatrixes
  DenseSubMatrix<Number>   Ke_sub(Ke_total);//matrix for master equation that couples u and additional variables
            
        
  DenseSubVector<Number>   Fe_sub(Fe_total);//RHS for the master equation

  DenseSubVector<Number>   Fe_add_sub(Fe_total);//RHS for the additional equation
   
  DenseSubMatrix<Number>   Ke_u_add_sub(Ke_total);//matrix for master equation that couples u and additional variables
  
  DenseSubMatrix<Number>   Ke_add_u_sub(Ke_total);//matrix for superalttice equation that couples u and additional variables

 
  DenseSubMatrix<Number>   Ke_add_add_sub(Ke_total);//matrix for superalttice equation that couples additional variables only
  //--------------------------------------------------------------



  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  
  std::vector<unsigned int> dof_indices_component;
 
  std::vector<unsigned int> dof_indices_total;

  

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  //  std::vector<int> mapping(dof_map.n_dofs());
  
  const   unsigned int Number_of_act_el = mesh.n_active_elem(); //number of fict DOFs

  const   unsigned int Number_of_DOFs = dof_map.n_dofs(); //number of total DOFs
  

  
  //dof_map.print_dof_constraints();   	
 


  SimulationEnvironment& si = get_environment(); 
  

  

  Tensor1 vec1;
  Tensor1 vec2, vec3;

  Tensor2Sym eps_var;

  Tensor2Sym eps_const;

  Tensor2Sym stress_converse_piezo;


  Tensor2Sym C_kl;

  // double a_substrate[3];
 
  double lattice_factor;


  
  unsigned int el_number = 0; 
  
  system.matrix->zero();
  system.rhs->zero();

  Stiffness* C_tensor_el;
 
  
  MacrostrainModel* macrostrain_model;

 

  for ( ; el != end_el ; ++el) 
  {//el
    
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem* elem = *el;
    
    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices (elem, dof_indices);
   
    dof_map.dof_indices (elem, dof_indices_component, uvar[0]);
    const unsigned int n_dofs   = dof_indices_component.size() * 3; //in fact, could be  dof_indices.size() - 1, fict is not used 

    fe->reinit  (elem);
  

    Ke_total.resize (n_dofs + number_of_add_var, n_dofs + number_of_add_var);
    Fe_total.resize (n_dofs + number_of_add_var);

    dof_indices_total.resize(n_dofs + number_of_add_var);

      
    ID subdomain = elem->subdomain_id();
      

      
    const Material* mat = _device->get_material(subdomain);

    const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());
      
	
    macrostrain_model = dynamic_cast<MacrostrainModel*>(   mat ->get_model(get_id())     );

    C_tensor_el = macrostrain_model->get_stiffness();
    

    eps_const =  crystal_el->get_const_eps0(substrate_lat_const, eps0_var_log) 
      + eps0_of_elem[el_number] ;//+ substrate_shear;
    

    macrostrain_model->get_converse_piezo_stress(stress_converse_piezo, elem);

    


    double lat_const[3];
    crystal_el->get_lat_const(lat_const);
     
    //-------------------------------------------------------------//
    //master equation:                                             //
    // without converse piezo effect                               //
    //     d/dx_i  ( C_ijkl (du_k/dx_l + eps0_kl)) = 0  
    // 
    // with converse piezo effect  
    //     d/dx_i  ( C_ijkl (du_k/dx_l + eps0_kl) - d_k,ij Ek) =   0
    //                                                             //
    //                                                             //
    //-------------------------------------------------------------//
    
    for (unsigned int j = 0; j<=2; j++)
    {//loop over j: master equation discretization


      //Right Hand Side---------------------------------------------------------------------//

	 
	  
	      
      dof_map.dof_indices (elem, dof_indices_component, uvar[j]);
      const unsigned int n_u_dofs = dof_indices_component.size(); 
      Fe_sub.reposition (uvar[j]*n_u_dofs, n_u_dofs);
      
      const unsigned int num_sides =  elem->n_sides(); 
      //!first we calculate volume part
      //  /
      //  | C_ijkl eps0_kl df/dx_i dV
      //  /
      //

      //!first we calculate volume part
      //  /
      //  | (C_ijkl eps0_kl + d_k,ij Ek) df/dx_i dV
      //  /
      //
      


      for (unsigned int p1=0; p1<n_u_dofs; p1++)
      {//nodes loop
	if (!belongs_to_substrate(p1, elem))
	{//not a substrate
	  for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	  {//qp
	    vec1 = 0;
	    for (int i = 1; i<=dim; i++) vec1(i) = dphi[p1][qp](i-1);
	    
	    //-------------eps0 part------------------
	    for (unsigned int k = 0; k <= 2; k++)
	    {
	      vec2 = 0;
	      vec3 = 0;
	      for (int i = 1; i <=3; i++ ) 
	      {	     	
		if (k+1 > i)
		{
		  vec2(i) = eps_const(k+1,i);
		  vec3(i) = stress_converse_piezo(k+1,i);
		}
		else
		{
		  vec2(i) = eps_const(i,k+1);
		  vec3(i) = stress_converse_piezo(i,k+1); 
		}
	      }
	      
	     
	      Fe_sub(p1) -= JxW[qp]*(vec1 * ( C_tensor_el->get_subtensor(j+1,k+1)* vec2 - vec3 ) ) ;
	    } 
	  }
	       
	    
	  

	  //! may be there is external pressure or extended device, so we have to calculate a surface part
	  //
	  //
	  //
	  for (unsigned int side=0; side<num_sides; side++)
	  {//side loop
	    Boundary* bd = si.get_boundary(std::pair<const Elem*,  unsigned int> (elem,side));
	    
	    
	    if (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  )
	      if (  dynamic_cast<MacrostrainBoundaryProperties*>( bd->get_boundary_properties( get_id() ))->get_type() == "pressure" ) 
	      { 	 
		MacrostrainPressure* press =
		  dynamic_cast< MacrostrainPressure* > (bd->get_boundary_properties (get_id()) ) ;
		if (dim > 1)
		{
		  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
		  
		  const std::vector<Real>& JxW_face = fe_face->get_JxW();
		  
		  const std::vector<Point >& qface_point = fe_face->get_xyz();
		  
		  const std::vector<Point> & normal = fe_face->get_normals();
			        
		  fe_face->reinit(elem, side);
			         
		  for (unsigned int qp=0; qp<qface.n_points(); qp++)
		  {
		    //why minus? because pressure = -stress (points into the region)
		    Fe_sub(p1) -= ((JxW_face[qp] * phi_face[p1][qp])
				   * press->get_value() ) * normal[qp](j);
		  } 
		  
		}
		else
		{
		  double normal;
		  Point p = elem->point(side);
		  Point pc = elem->centroid();
		  if (p(0) > pc(0)) 
		    normal = 1.0;
		  else
		    normal = -1.0;
		  
		  
		  Fe_sub(p1) -=  press->get_value() * normal;
		}
		  
	      }
	      else if (dynamic_cast<MacrostrainBoundaryProperties*>( bd->get_boundary_properties( get_id() ))->get_type() == "extended")
	      {
		if (dim > 1)
		{
		  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
		  
		  const std::vector<Real>& JxW_face = fe_face->get_JxW();
		  
		  const std::vector<Point >& qface_point = fe_face->get_xyz();
		  
		  const std::vector<Point> & normal = fe_face->get_normals();
			        
		  fe_face->reinit(elem, side);
			         
		 
		  
		} 
		else
		{
		  
		}
	      }
	  }//end of side loop
	}
	else
	{
	  Fe_sub(p1) = 0.0;
	}	      
      }
   
	    	
			
      //---RHS of master equation is done---------------------------------------------------//
      //---now we do the matrix-------------------------------------------------------------//    
      
      //----- additional variables first, if any--------------------------------------------//
      for (unsigned int p1=0; p1<n_u_dofs; p1++)
      {
	if (!belongs_to_substrate(p1, elem))
	{
	  for (unsigned int i1 = 0; i1 < add_var.size()  ; i1++)
	  {
	    Ke_u_add_sub.reposition (uvar[j]*n_u_dofs, n_dofs, n_u_dofs, add_var.size());	  
	    
	    eps_var = crystal_el->get_var_eps0( add_var[i1].name );
	    
	    for (unsigned int qp=0; qp<qrule.n_points(); qp++) 
	    {
	      vec1 = 0;
	      for (int i = 1; i<=dim; i++) vec1(i) = dphi[p1][qp](i-1);
	      
	      for (unsigned int k = 0; k <= 2; k++)
	      {
		vec2 = 0;
		for (int i = 1; i <=3; i++ ) 
		{  
		  if (k+1 > i)
		    vec2(i) = eps_var(k+1,i);
		  else
		    vec2(i) = eps_var(i,k+1);
		}
		
		Ke_u_add_sub(p1,i1) += JxW[qp]*(vec1 * ( C_tensor_el->get_subtensor(j+1,k+1) * vec2  )) ;
	      }
	    }
	  } 
	}
      }
      //------ additional variables done --------------------
      //------ now normal variables -------------------------
      
      for (unsigned int k = 0; k<=2; k++)
      {//loop over k
	
	dof_map.dof_indices (elem, dof_indices_component, uvar[k]);
	const unsigned int n_u_dofs = dof_indices_component.size(); 
	Ke_sub.reposition (uvar[j]*n_u_dofs, uvar[k]*n_u_dofs, n_u_dofs, n_u_dofs);
	
	
	for (unsigned int p1=0; p1<n_u_dofs; p1++)
	{
	  for (unsigned int p2=0; p2<n_u_dofs; p2++)
	  {		      
	    double scal_prod;
	    
	    
	    for (unsigned int qp=0; qp<qrule.n_points(); qp++) 
	    {
	      vec1 = 0;
	      for (int i = 1; i<=dim; i++) vec1(i) = dphi[p1][qp](i-1) ;
	      
	      vec2 = 0;
	      for (int i = 1; i<=dim; i++) vec2(i) = dphi[p2][qp](i-1) ;
	      
	      scal_prod = vec1 * ( C_tensor_el->get_subtensor(j+1,k+1) *vec2) ;
	      
	      if (!belongs_to_substrate(p1, elem))
	      {
	
		 Ke_sub(p1,p2) += JxW[qp]*scal_prod ;
	
	      }
	      else
	      {
		Ke_sub(p1,p2) = delta(p1,p2)*delta(j,k);
	      }
	    }
	    
	  }
	  
	}
	
		
      }
		      
		      
    }//end of loop over j - end of master equation
    //---------------------------------------------------------------------------//
    //master equation is done                                                    //
    //---------------------------------------------------------------------------//
    //superlattice equations
    //---------------------------------------------------------------------------//
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)//qp
    {
      for (unsigned int eq_number =  0;   eq_number <  number_of_add_var; eq_number++)
      {
	
	dof_map.dof_indices (elem, dof_indices_component, uvar[0]);
	const unsigned int n_u_dofs = dof_indices_component.size();
	
	  
	if ( add_var[eq_number].lat_cons )
	{
	  unsigned int lat_index = add_var[eq_number].index1;
	  double lat_constants[3];
	  crystal_el-> get_lat_const(lat_constants);
	  
	      
	  double lat_const = lat_constants[lat_index - 1];
	      
	  lattice_factor = 1/lat_const;
	}
	else
	{
	  lattice_factor = 1.0;
	}
		  
	unsigned int lat_index1 = add_var[eq_number].index1;
	unsigned int lat_index2 = add_var[eq_number].index2;
	  
	  
	C_kl = C_tensor_el->get_another_subtensor(lat_index1,lat_index2);
	//----------------RHS------------------
	Fe_add_sub.reposition(n_dofs + eq_number,1);
	
	  
	Fe_add_sub(0) -=  JxW[qp] * doubleContraction(C_kl -  stress_converse_piezo, eps_const ) * lattice_factor  ;
	  
	//-------------------------------------
	  
	//-------------ux,uy,uz----------------
	for (unsigned int k = 0; k<=2; k++)
	{//loop over k
	  dof_map.dof_indices (elem, dof_indices_component, uvar[k]);
	  const unsigned int n_u_dofs = dof_indices_component.size(); 
	  Ke_add_u_sub.reposition (n_dofs, uvar[k]*n_u_dofs, add_var.size() , n_u_dofs);
	  
	  vec1 = 0;
	  for (unsigned int l = 1; l <= dim; l++)
	  {
	    if (l>= k+1 ) 
	      vec1(l) = C_kl(l,k+1);
	    else 
	      vec1(l) = C_kl(k+1,l);
	  }
	  
	  for (unsigned int p1=0; p1<n_u_dofs; p1++)
	  {
	    vec2 = 0; 
	    for (unsigned int l = 1; l <= dim; l++) vec2(l) = dphi[p1][qp](l-1) ;
	    
	    Ke_add_u_sub(eq_number, p1) += JxW[qp] * (vec1*vec2) * lattice_factor ;
	    
	  }
	}
	//------------------------------------------
	
	//-----------add_add---matrix---------------
		 
	for (unsigned int i1 = 0; i1 < add_var.size()  ; i1++)
	{
	  Ke_add_add_sub.reposition(n_dofs + eq_number,n_dofs + i1,1,1);
	  eps_var = crystal_el->get_var_eps0( add_var[i1].name );
	  
	  Ke_add_add_sub(0,0) += JxW[qp]  *  doubleContraction(eps_var,C_kl)   *  lattice_factor;
	  

	}
      }
      
    }  

    //------------------------------------------

		
	    
    


	  
    // end of superlattice equations
    //--------------------------------------------------------------------------

	  
   
     
       
    for (unsigned i =0 ; i < n_dofs; i++)
      dof_indices_total[i] = dof_indices[i];
  
       
    for (unsigned i =0 ; i <number_of_add_var ; i++)
      dof_indices_total[i+n_dofs] =  add_dofs_vector[i];
	
       
       

      
    dof_map.constrain_element_matrix_and_vector(Ke_total, Fe_total, dof_indices_total);
    
  
    system.matrix->add_matrix (Ke_total, dof_indices_total);
    system.rhs->add_vector    (Fe_total, dof_indices_total);
      
    // constraint is necessary for Ke_add!!
      

   

    el_number++;
  }


  //-------------------------------------------------------//
  //diagonal of the unused ficticious variables 
  if  (number_of_add_var != 0)
  {
    el     = mesh.active_elements_begin();
    el_number = 0;
    for ( ; el != end_el ; ++el) 
    { 
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
     
      unsigned int n = dof_indices.size();
      Ke_total.resize (n,n);


      if ( el_number >=  number_of_add_var )
      {
	Ke_sub.reposition(n - 1,n - 1, 1, 1);
	Ke_sub(0,0) += 1.0;
      }

      system.matrix->add_matrix (Ke_total, dof_indices);
      
      el_number++;
    }

  }

     
/*
{
  system.matrix->close();
  
  system.matrix->print_matlab("matr.m");

  system.rhs->close();

  system.rhs->print_matlab("rhs.m");
}
*/
  //-----------------------------------------------------------------------
  //Application of periodicity constraints



  //-----------------------------------------------------------------------
  //dof_map.print_dof_constraints(); 	  	

//   system.matrix->print();

   
 

#ifdef DEBUG
  std:: cout<< "matrix is done \n";  
   
  // system.rhs->print();


  std:: cout << "Active dofs number   " << system.n_active_dofs()   	<< "\n";

  std:: cout << "Total dofs number   " << system.n_dofs()   	<< "\n";

  std:: cout << "Constraint dofs number " <<  system.n_constrained_dofs()   	<< "\n";
#endif
      


}




//-----------------------------------------------------------------//
 bool Macrostrain::belongs_to_substrate(unsigned int n, const Elem* elem )
 {



   if (grown_on_substrate)
     {
     
       const Node* nd = elem->get_node(n);

       std::set <const Node*> :: iterator it = substrate_points.find( nd );

       if ( it != substrate_points.end() )
	 return(true);
       else
	 return (false);      

        
     }
  

   else
     {

       //there is always one node that is fixed. 
       // std :: cerr << elem->node(n) << "  " << fixed_node_number << "\n";
       //return(elem->node(n) == fixed_node_number);

       if (elem->node(n) == fixed_node1 ) 
	 {
	  
	   return(true);
	   
	 }
       else 	
	 return(false);
     }

 
 }

//-----------------------------------------------------------------//
 inline int Macrostrain::delta (int i, int j)
{
  return ((i==j) ? 1 : 0);
}
//-----------------------------------------------------------------//
void Macrostrain::refer_objects()
{


  static_this = this; 


}
//-----------------------------------------------------------------//
Mesh* Macrostrain::get_mesh()
{
  return( &(equation_systems->get_mesh()) );
}

//-----------------------------------------------------------------//
void Macrostrain::do_solve()

{
  

  parse_options();


  SimulationEnvironment& si = get_environment();   

 

  //------------------------------------------------------
  NumericVector<Number>& old_solution = 
    my_system->get_vector("old solution");
  

 

  //------------------------------------------------------

 
  Mesh& mesh = equation_systems->get_mesh();

 

 
  initialize_eps0_list();
  initialize_el_number_map();

  define_fixed_nodes();

  make_nodes_periodic();

  //init_u_node(); //not necessary
  
  init_substrate();



  refer_objects();

  if (grown_on_substrate) si.get_boundary_nodes (substrate_name, substrate_points);

  set_up_additional_dofs();

  my_system->solution->zero();

  apply_periodic_bc();

  apply_antirotation_constraints();

  
   

  my_system->solve();

 
  
  old_solution = * (my_system->solution);
 
  update_substrate();


  if (intermediate_output)
    {
      
      if (output_type=="GMV") GMVIO (mesh).write_equation_systems ("displacement_field.dat.000", *equation_systems);
      if (output_type=="tecplot") TecplotIO_cell(mesh,false).write_equation_systems ("displacement_field.dat.000", *equation_systems);

      output_strain("strain.dat.000");
      output_add_strain_variables("add_var.000");
    }


  if (dim > 1) mesh.write("mesh0.ucd");
  //-----------------------------------------------------
 
 


  //refinement loop----------------------------------------------------------------
  MeshRefinement mesh_refinement(mesh);

  for (unsigned int r_step = 1; r_step <= max_r_steps; ++r_step)
    {
      std::cerr << "\nRefining the mesh... (Step" << r_step << ")\n" << std::endl;
      

      

     
      ErrorVector error;
      
      KellyErrorEstimator error_estimator;

     
      error_estimator.estimate_error (*my_system,error);
		
      
      
      mesh_refinement.flag_elements_by_error_fraction (error,
						       refine_fraction,
						       coarsen_fraction,
						       max_ref_level);
      

   
	
      // This call actually refines and coarsens the flagged
      // elements.
      if (uniform_refinement == 1)
	mesh_refinement.uniformly_refine(1);
      else
        mesh_refinement.refine_and_coarsen_elements();

     
      equation_systems->reinit();

      old_solution = *(my_system->solution);
      
      old_solution.close();
   
      initialize_eps0_list();

      initialize_el_number_map();
      
      set_up_additional_dofs();
      
      init_substrate();

      define_fixed_nodes();

      if (grown_on_substrate) si.get_boundary_nodes (substrate_name, substrate_points);

      refer_objects();

      make_nodes_periodic();
      
      apply_periodic_bc();

      apply_antirotation_constraints();
      
      mesh.print_info();
      
      my_system->solution->zero();
      
      my_system->solve();
      

      std::cout << "Norm of the difference  " << norm_of_difference( old_solution, *(my_system->solution) ) 
		<< "  after step number " << r_step << "\n";
     
      if (intermediate_output)
      {      
	  
	std::ostringstream os;
	os << "displacement_field.dat.00" << r_step;


	if (output_type=="GMV")  GMVIO (mesh).write_equation_systems (os.str(), *equation_systems);
	if (output_type=="tecplot")   TecplotIO_cell(mesh,false).write_equation_systems(os.str(), *equation_systems);

	std::ostringstream os_mesh;
	os_mesh << "mesh" << r_step << ".ucd";
	if (dim > 1) mesh.write(os_mesh.str());
	  
	std::ostringstream os1;
	os1 <<"strain.dat.00" << r_step;
	output_strain(os1.str());	 
	  
	std::ostringstream os2;
	os2 <<"add_var.00" << r_step;
	  
	output_add_strain_variables(os2.str());
      }

      std::cout << "\n" ;
      std::cout << "Final Mesh after  " <<  max_r_steps <<" refinements  steps   " <<  "\n" ;
      mesh.print_info();
      std::cerr << "Grid refinement is done \n";
    }


 

  //cerr << atom_structure_filename << "\n";

  if (calculate_atom_displacements)
    { 
      read_atom_structure(atom_structure_filename);
      
      std::ostringstream disp_file;
      disp_file << atom_displacements_filename << ".out"  ;
     
   

      write_atom_displacements( disp_file.str());
    }
  //------------------------------------------------------------------------------------
  //geometry relaxation

 

 

  for (unsigned int geom_it = 1 ; geom_it <= max_shape_steps; geom_it++)
    {
      
      if (geom_it > 1) update_u_node();
      
      equation_systems->print_info();
      //------move nodes------------------------------------------
	
      update_eps0_list();


      refer_objects();

      move_nodes();
      
      update_substrate();

      
      //---------------------------------------------------------
     
      //-------solve---------------------------------------------
      my_system->solution->zero();
      
      //apply_periodic_bc();
      
      equation_systems->print_info();
      
      my_system->solve();
      if (intermediate_output)
	{

	  //------write-----------------------------------------------
	  if (dim > 1) mesh.write("mesh0.ucd");
	  
	  std::ostringstream os;
	  os << "displacement_field.dat.00" << geom_it + max_r_steps;
	  if (output_type=="GMV")  GMVIO (mesh).write_equation_systems (os.str(), *equation_systems);
	  if (output_type=="tecplot") TecplotIO_cell(mesh,false).write_equation_systems (os.str(), *equation_systems);

	  std::ostringstream os_mesh;
	  os_mesh << "mesh" << geom_it + max_r_steps<< ".ucd";
	  if (dim > 1) mesh.write(os_mesh.str());
      
	  std::ostringstream os1;
	  os1 <<"strain.dat.00" << geom_it;
	  output_strain(os1.str() );

	  
	  std::ostringstream os2;
	  os2 <<"add_var.00" << geom_it + max_r_steps;

	  output_add_strain_variables(os2.str());
	  
	  if (calculate_atom_displacements)
	    {
	      std::ostringstream disp_file;
	      disp_file << atom_displacements_filename << geom_it <<".out";

	     
	      write_atom_displacements( disp_file.str() );
	    }

	}
            //---------------------------------------------------

     
      

    }//end of shape loop
    

 

  if (max_shape_steps >= 1) update_u_node(); //temporary

  update_u_node();
 
  //------write-------------------------------------------------------------------------------------//
  //--  output of the final result
  if (intermediate_output)
  {
    if (dim > 1) mesh.write("mesh0.ucd");

	  
    std::ostringstream os;
    os << "displacement_field.dat" ;

    if (output_type=="GMV") GMVIO (mesh).write_equation_systems (os.str(), *equation_systems);
    if (output_type=="tecplot")   TecplotIO_cell(mesh,false).write_equation_systems (os.str(), *equation_systems);
    
    std::ostringstream os_mesh;
    os_mesh << "mesh"<< ".ucd";
    if (dim > 1) mesh.write(os_mesh.str());
      
    std::ostringstream os1;
    os1 <<"strain.dat" ;
    output_strain(os1.str() );
    
      
    std::ostringstream os2;
    os2 <<"add_var";
      
    output_add_strain_variables(os2.str());
      
   
  }


  if (calculate_atom_displacements)
  {
    std::ostringstream disp_file;
    disp_file << atom_displacements_filename <<".out";


   

    write_atom_displacements(disp_file.str() );
  }

  //--------------------------------------------------------------------------------------------------//
  calculate_result_elem_strain_map();


 
 
  //--------------------------------------------------------------------------------------------------//
 }


//-----------------------------------------------------------------//
void Macrostrain::update_eps0_list()
{
  //calculate eps_new = eps_old + 1/2(du/dx + du/dx) 

  const Mesh& mesh = equation_systems->get_mesh();
 
  LinearImplicitSystem& system = *my_system;

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  DofMap& dof_map = system.get_dof_map();
  
  unsigned int uvar[3] ;
  
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
 

 
  

  FEType fe_type = dof_map.variable_type(0);
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type, true));//no scaling here!
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  
  std::vector<Point> point_vec(1);

  std::vector<unsigned int> dof_indices_component1;
  std::vector<unsigned int> dof_indices_component2;
  
  MeshBase::const_element_iterator el  = mesh.active_elements_begin();
  MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  unsigned int el_number = 0;
  Tensor2Sym eps1;
  for ( ; el != end_el ; ++el) 
    {
     
      const Elem* elem = *el;

      point_vec[0]  = FEInterface::inverse_map(dim, fe_type, elem, elem->centroid());
      fe->reinit (elem, &point_vec);

     
      for (int i = 1; i <=3 ; i++)
       for (int j = 1; j <=i; j++)
	 {
	   double du_i_over_dx_j = 0;  
	    
	  
           dof_map.dof_indices (elem, dof_indices_component1, uvar[i-1]);
	   dof_map.dof_indices (elem, dof_indices_component2, uvar[j-1]);

	   const unsigned int n_u_dofs = dof_indices_component1.size(); 
	   
	   for (unsigned int p1=0; p1<n_u_dofs; p1++)
	     {   
	       
	       if (j<=dim) du_i_over_dx_j += 0.5 *  dphi[p1][0](j-1) * (*solution)(dof_indices_component1[p1]);
	       
		
	       if (i<=dim) du_i_over_dx_j += 0.5 *  dphi[p1][0](i-1) * (*solution)(dof_indices_component2[p1]); 
	       

	     }
	    
	   eps1(i,j) = du_i_over_dx_j;
	 }
     
      eps0_of_elem[el_number] += eps1;

     
      el_number++;
    }
  
}

//----------------------------------------------------------------------//
//--------------------------------------------------------------------
void Macrostrain::initialize_eps0_list()
{

  //initialize vector with lattice matching strain

  const Mesh& mesh = equation_systems->get_mesh();

  const unsigned int N_elem = mesh.n_active_elem();

  eps0_of_elem.resize( N_elem, Tensor2Sym(0) );




}
//------------------------------------------------------------------------
void Macrostrain::initialize_el_number_map()
{
  const Mesh& mesh = equation_systems->get_mesh();
  const unsigned int N_elem = mesh.n_active_elem();

  elem_numbers.clear();
  unsigned int el_number = 0;

  Elem* elem;

  MeshBase::const_element_iterator el  = mesh.active_elements_begin();
  MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  
  for ( ; el != end_el ; ++el) 
    {
      elem = *el;
      
      elem_numbers.insert(pair<Elem*,unsigned int> (elem, el_number));
     

      el_number++;

    }
  
}
//--------------------------------------------------------------------------
//------------------------------------------------------------------------

void Macrostrain::make_nodes_periodic()
{
  const double pos_tol = 1e-10;
  const Mesh& mesh = equation_systems->get_mesh();
  nodes_periodic.clear();

  const Node* node_fix = &mesh.node(fixed_node1);

  

  for (unsigned dir = 0; dir <=dim-1; dir++)
    {//directions
      std::vector < const Node*> temp_vec;
      temp_vec.clear();
    
      if (periodicity[dir]) 
	{
	  //----------------------------------------------------------------
	  //check if the fixed_point 1 is good
	  if(!grown_on_substrate)
	    if ( ( std::abs( (*node_fix)(dir) - max_coord[dir]) < pos_tol)  ||
		 ( std::abs( (*node_fix)(dir) - min_coord[dir]) < pos_tol)  )
	    {
	      cerr << "dir  " << dir << "\n";
	      
	      cerr << (*node_fix)(dir) << "   " << max_coord[dir] << "\n";
	      
	      cerr << "Fixed node 1 is on a boundary that periodic boundary conditions are applied for.\n\n This should be changed.\n";
			 exit(1);
			 
	    }
	  //------------------------------------------------------------------
	  for (unsigned int n = 0; n < mesh.n_nodes(); n++) // Loop over all the nodes
	    {
	      const Node* node1 = &mesh.node(n);
	      if (node1->active())
		{		
		  if ( std::abs( (*node1)(dir) - min_coord[dir]) < pos_tol)  temp_vec.push_back(node1);
		}
	    }
	}
      nodes_periodic.push_back(temp_vec);
    }
}

//-------------------------------------------------------------------------


void  Macrostrain::apply_periodic_bc()
{

  // It is a good idea to make sure we are assembling
  // the proper system.
 

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Periodic bc. Assembly",false);


  
  // Get a constant reference to the mesh object.
  const Mesh& mesh = equation_systems->get_mesh();

  // The dimension that we are running
  //const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem& system = *my_system;

  unsigned int uvar[3] ;
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }

  unsigned int system_number=system.number();
  
  DofMap& dof_map = system.get_dof_map();
  
  FEType fe_type = dof_map.variable_type(uvar[0]);
  
 
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type, true));  //no scaling here
   

  

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

 
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_component;
  
  const double pos_tol = 1e-10;
  const double func_tol = 1e-10;
 
  //dof_map.print_dof_constraints();  

  for (int i = 0; i < mesh.mesh_dimension(); i++) //Loop over all the mesh directions
    { 
      if  (periodicity[i]) //Check if the periodic b.c. are applied along the direction i
	
	{
	 
	  std :: vector <const Node*>& vec =  nodes_periodic[i];

	  for (unsigned int n = 0; n < vec.size(); n++) // Loop over all the nodes
	    {
	      const Node* node1 = vec[n];
	      
	    
		 
	      for (unsigned int var_index = 0 ; var_index  <= 2;  var_index ++)
		{//let us find dof for it-----------------
		  

		 // const Node& node = mesh.node(n);
		  
		  const unsigned int  n_dof = node1->dof_number(system_number,uvar[var_index],0);
		  
		  //dof is found-------------------------------
		 
		  
		  if (! dof_map.is_constrained_dof(n_dof) ) //only if the dof is not constrained do the job
		    {
		      //let us make a  point that lies at the opposite side
		      Point point2(*node1);
		
		      point2(i) = point2(i) + max_coord[i] - min_coord[i];
		  
		      
		      //corresponding point is created
		      
		      
		      //let us find an element this point belongs to and calculate the constraints
		      //the most coarse element first
		      unsigned int refinement_level = 0; 
		      MeshBase::const_element_iterator el3  = mesh.level_elements_begin(refinement_level);
		      MeshBase::const_element_iterator end_el3 = mesh.level_elements_end(refinement_level);
		      
		      Elem*  elem1;
		      for ( ; ( (el3 != end_el3) ) ; ++el3)  
			{
			  Elem* elem = *el3;
			  if (element_on_boundary(elem))
			    {
			      if (elem->contains_point(point2))
				{
				  elem1 = elem;
				  
				  break;
				  
				}
			    }
			}
		      
		      //children of the  most coarse element 
		      while ( !( elem1->active() ) )
			{
			  
			  for (unsigned int i=0 ; i < elem1->n_children() ; i++)
			    {
			      Elem* 	child = elem1->child(i);
			      if (element_on_boundary(child))
				{
				  if (child->contains_point(point2))
				    {
				      elem1 = child;
				      break;
				    }
				}
			    }
			}
			
		      

		      
		      //active elem1 contains the opposite  point, we can constrain it now
		      
		      DofConstraintRow constraint;
		      constraint.clear();
		      //dof_map.dof_indices (elem1, dof_indices);
		      dof_map.dof_indices (elem1, dof_indices_component, uvar[var_index]);
		      
		      std::vector<Point> point2_vec(1);
		      
		      point2_vec[0] = point2;
		      
		      std::vector<Point> point2_ref_vec(1);
			  
			  
		      FEInterface::inverse_map (elem1->dim(), fe_type , elem1,  point2_vec,  point2_ref_vec)  ;
		      
		      fe->reinit (elem1, &point2_ref_vec);

		      Point point_temp = point2_ref_vec[0];
		     
		      
		      for (int i1 = 0; i1 < phi.size(); i1++)
			{
			
			  if ( std::abs(phi[i1][0]) >  func_tol )  
			    {
			     
			      constraint[dof_indices_component[i1]] = phi[i1][0];
			    }
			}
		       

		      dof_map.add_constraint_row (n_dof,  constraint); 
		      
		    }
		     
		    
		}
	    }
	  
	  
	}
    }  
  //  std::cout << '+++++++++++++++++\n';
  // dof_map.print_dof_constraints();  
  //  std::cout << '+++++++++++++++++\n';

}


//---------------------------------------------------------------------------------//
void Macrostrain::apply_antirotation_constraints()
{
  /*
This method is to apply constraints that are necessary for a freestanding structure simulation
The constrants are the following:

             *p2
            /
           /
          /
         /______________
      p1*               *p3



     1)Point p1 is fixed (considered as a substrate point in a another method)
     2)Displacement vector of the point p2 is parallel to the [p1,p2] vector
     3)Displacement vector of the point p3 belongs to (p1,p2,p3) plane    

     Notes:
     a) 1D case does not require the antirotation constraints 2) and 3)
     b) 2D case does not require the antirotation constraint 3)
     c) Dimension of the problem is defined as:
     Geometric dimension - number_of_periodic_directions


  */


 const Mesh& mesh = equation_systems->get_mesh();
 // Get a reference to the LinearImplicitSystem we are solving
 LinearImplicitSystem& system = *my_system;

 unsigned int uvar[3] ;
 for (unsigned int i = 0; i<= 3 - 1; i++) 
   {
      uvar[i] = system.variable_number(uname_vec[i]);
   }
 
 unsigned int system_number=system.number();
 
 DofMap& dof_map = system.get_dof_map();


  
 FEType fe_type = dof_map.variable_type(uvar[0]);
  
 
 //AutoPtr<FEBase> fe (build_finite_element(dim, fe_type, true));

 DofConstraintRow constraint; 




 if ( !grown_on_substrate )
   {//substrate system is not treated
     if (dim > 1)
       {//1D system does not need a treatment 


	 if (fixed_node1 == fixed_node2) 
	   {
	     
	     cerr << "Error: fixed_node1 == fixed_node2  " << fixed_node1 <<"  "
		  << fixed_node2 <<"\n";  
	     cerr << fixed_node1 <<"  "<< fixed_node2 <<"\n";   
	     exit(1);
	   }

	 const Node& node1= mesh.node(fixed_node1);
	 const Node& node2= mesh.node(fixed_node2);

        



	 vector<double> nd1_to_nd2(3);
	 for (int i = 0; i < 3; i++) nd1_to_nd2[i] = node2(i) - node1(i);
	 
	 unsigned int dir1;
	 unsigned int dir2;
	 unsigned int dir3;
	 
	 if (dim == 2)
	   {

	     if (periodicity[0] || periodicity[1])
	       {
		 cout << "Since you have periodic boundary conditions the antirotation constraints are ignored\n";
	       }
	     else

	       {

		 if  (abs(nd1_to_nd2[0]) > 1e-5) 
		   {
		     dir1 = 0;
		     dir2 = 1;
		   }
		 else
		   { 
		     dir1 = 1;
		     dir2 = 0;
		   }

		 const unsigned int  n_dof_constrained = node2.dof_number(system_number,uvar[dir2],0);
		 const unsigned int  n_dof_used = node2.dof_number(system_number,uvar[dir1],0);
		 DofConstraintRow constraint;
		 constraint.clear();
		 constraint[n_dof_used] = nd1_to_nd2[dir2]/nd1_to_nd2[dir1];
		 dof_map.add_constraint_row(n_dof_constrained,  constraint);
		 
	       }
	   }
	 if (dim == 3)
	   {//dim ==3
	     
	    //calculate number of periodic bc
	     unsigned int n_per = 0;
	     for (unsigned int i = 0; i < 3; i++)
	       {
		 if (periodicity[i]) n_per++;
	       }
	     

	     if (n_per > 1)
	       {
	
		 cout << "Since you have periodic more than boundary conditions the antirotation constraints are ignored\n";
	       }
	     else
	       {//---

		 if ((nd1_to_nd2[0]*periodicity[0] + nd1_to_nd2[1]*periodicity[1] + nd1_to_nd2[2]*periodicity[2])/
		     ( sqrt(nd1_to_nd2[0]* nd1_to_nd2[0]+ nd1_to_nd2[1]*nd1_to_nd2[1] + nd1_to_nd2[2] * nd1_to_nd2[2] )) < 1e-4)

		   {
		   
		     if  (abs(nd1_to_nd2[0]) > 1e-5) 
		       {
			 dir1 = 0;
			 dir2 = 1;
			 dir3 = 2;
		       }
		     
		     if  (abs(nd1_to_nd2[1]) > 1e-5) 
		       {
			 dir1 = 1;
			 dir2 = 0;
			 dir3 = 2;
		       }
		 
		     if  (abs(nd1_to_nd2[2]) > 1e-5) 
		       {
			 dir1 = 2;
			 dir2 = 0;
			 dir3 = 1;
		       }
		     
		    
		     const unsigned int  n_dof_constrained1 = node2.dof_number(system_number,uvar[dir2],0);
		     const unsigned int  n_dof_constrained2 = node2.dof_number(system_number,uvar[dir3],0);
		     const unsigned int  n_dof_used = node2.dof_number(system_number,uvar[dir1],0);
		     
		     
		     constraint.clear();
		     constraint[n_dof_used] = nd1_to_nd2[dir2]/nd1_to_nd2[dir1];
		     dof_map.add_constraint_row(n_dof_constrained1,  constraint);
		     
		     constraint.clear();
		     constraint[n_dof_used] = nd1_to_nd2[dir3]/nd1_to_nd2[dir1];
		     dof_map.add_constraint_row(n_dof_constrained2,  constraint);
		   }
		 else
		   {
		     cerr << "Direction between fixed point 1 and 2 can not has to be  perpendicular to the  parallel direction\n";
		     exit(1);
		   }
		 
		 //---------------------------------------------------------------
		 //in a 3D case we need additional constrain
		 if (n_per == 0)
		   {
		     if ((fixed_node1 == fixed_node3)|| (fixed_node2 == fixed_node3)) 
		       {
		     
			 cerr << "Error: fixed_node1_temp == fixed_node3_temp  or \n";
			 cerr << "Error: fixed_node2_temp == fixed_node3_temp   \n";
			 cerr << fixed_node1 <<"  "<< fixed_node2 << "  "<< fixed_node2  << "\n";   
			 exit(1);
		       }
		 
		     const Node& node3= mesh.node(fixed_node3);
		     vector<double> nd1_to_nd3(3);
		     for (int i = 0; i < 3; i++) nd1_to_nd3[i] = node3(i) - node1(i);
	     
		 
		     vector<double> ortogon_to_plane(3); //ortogon_to_plane(3) = vect_prod( nd1_to_nd2, nd1_to_nd3)
		 
		     ortogon_to_plane[0] =  nd1_to_nd2[1]*nd1_to_nd3[2] - nd1_to_nd2[2]*nd1_to_nd3[1];
		     ortogon_to_plane[1] = -nd1_to_nd2[0]*nd1_to_nd3[2] + nd1_to_nd2[2]*nd1_to_nd3[0];
		     ortogon_to_plane[2] =  nd1_to_nd2[0]*nd1_to_nd3[1] - nd1_to_nd2[1]*nd1_to_nd3[0];


		     if (abs(ortogon_to_plane[0])>1e-5)
		       {
			 dir1 = 0;
			 dir2 = 1;
			 dir3 = 2;
		 
		       }
		     else
		       {
			 if (abs(ortogon_to_plane[1])>1e-5)
			   {
			     dir1 = 1;
			     dir2 = 2;
			     dir3 = 0;
			   }
			 else
			   {
			     dir1 = 2;
			     dir2 = 1;
			     dir3 = 0;
			   }
		       }
	     
	     

		     const unsigned int  n_dof_constrained = node3.dof_number(system_number,uvar[dir1],0);
		     const unsigned int  n_dof_used1 = node3.dof_number(system_number,uvar[dir2],0);
		     const unsigned int  n_dof_used2 = node3.dof_number(system_number,uvar[dir3],0);
	     
		     constraint.clear();
		     constraint[n_dof_used1] = - ortogon_to_plane[dir2]/ortogon_to_plane[dir1];
		     constraint[n_dof_used2] = - ortogon_to_plane[dir3]/ortogon_to_plane[dir1];
		     dof_map.add_constraint_row(n_dof_constrained,  constraint);
		     
		     //--------------------------------------------------------------
		   }
		 else
		   {
		     cout << "Since you have 1 periodic boundary conditions only 2 antiritation constraints are added\n";
		   }

	       }

	   }

	   

       }

   }

 

}
//--------------------------------------------------------------------------------//
void Macrostrain::prepare_strain_data_for_output( std::vector<std::string>& eps_names, std::vector<double>& eps_data ) 
{
  char num_i[2];
  char num_j[2];
  string eps_ij;
  
  
 
  double a_substrate[3];  substrate_crystal -> get_lat_const(a_substrate);

  unsigned int index = 0;

  const Mesh& mesh = equation_systems->get_mesh();

  //const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = *my_system  ;

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  DofMap& dof_map = system.get_dof_map();

  unsigned int uvar[3] ;
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
 
  

  FEType fe_type = dof_map.variable_type(0);
 
  
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type, true)); //no scaling here

 

  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();


  unsigned int Number_of_elements = mesh.n_active_elem();

  eps_data.resize(Number_of_elements*6);
  eps_names.resize(6);




  std::vector<Point> point_vec(1);

  

  std::vector<unsigned int> dof_indices_component1;
  

  std::vector<unsigned int> dof_indices_component2;

  Tensor2Sym eps0;

  // double der_tranf[3][3];

 

  for (int i = 1; i <=3 ; i++)
    for (int j = 1; j <=i; j++)
      {
	sprintf( num_i, "%i",i);
	sprintf( num_j, "%i",j);
	

	eps_ij = "eps_" + string(num_i) + string(num_j);

	eps_names[index] = eps_ij ;



	MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

	unsigned int elem_number = 0;

	for ( ; el != end_el ; ++el) 
	  {
	    const Elem* elem = *el;

	    Point center=elem->centroid();

	   
	    
	    eps0 = eps0_of_elem[elem_number] + calculate_eps_lat_matching( elem ); //previous iterations and lattice matching
	
	    Point center1 = FEInterface::inverse_map(dim, fe_type, elem, center);
	    point_vec[0] = center1;
	    fe->reinit (elem, &point_vec);
	    
	   
	    double du_i_over_dx_j = 0;
	 	    
	  
	    dof_map.dof_indices (elem, dof_indices_component1, uvar[i-1]);
	    dof_map.dof_indices (elem, dof_indices_component2, uvar[j-1]);

	    const unsigned int n_u_dofs = dof_indices_component1.size(); 
	  
	    for (unsigned int p1=0; p1<n_u_dofs; p1++)
	    {   
	      
	      if (j<= dim) 
	      {
		du_i_over_dx_j += 0.5 *  dphi[p1][0](j-1) * (*solution)(dof_indices_component1[p1]);
	      }
	      
	      if (i<= dim)
	      {
		du_i_over_dx_j += 0.5 *  dphi[p1][0](i-1) * (*solution)(dof_indices_component2[p1]); 
	      }
	     
	

	    }
	    
	    
	    double eps_value = eps0(i,j) + du_i_over_dx_j ;  
	    
	  
	    eps_data[index + elem_number * 6  ] = eps_value; //that's a correct order of variables
	     
	    elem_number++;
	  }

	index++;
      }
}

//--------------------------------------------------------------------------------//
//write out strain tensor component------------------------     ---------------
void Macrostrain::output_strain(std::string filename )
{

  std::vector<std::string> eps_names;

  std::vector<double> eps_data;

  prepare_strain_data_for_output(  eps_names,  eps_data );

  const Mesh& mesh = equation_systems->get_mesh();

  if (output_type == "GMV")     GMVIO_cell(mesh).write_ascii_cell_data(filename, eps_data, eps_names);

  if (output_type == "tecplot") TecplotIO_cell(mesh,false).write_cell_data(filename,eps_data,eps_names);

}
//---------------------------------------------------------------------------

void Macrostrain::prepare_polarization_data_for_output( std::vector<std::string>& polariz_names, std::vector<double>& polariz_data )
{

  char num_i[2];
  const Mesh& mesh = equation_systems->get_mesh();

  unsigned int Number_of_elements = mesh.n_active_elem();

  polariz_data.resize(Number_of_elements*3);

  polariz_names.resize(3);

  polariz_names[0] = "Px";
  polariz_names[1] = "Py";
  polariz_names[2] = "Pz";


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  unsigned int el_number = 0;

  Tensor1 polariz_vec;

  for ( ; el != end_el ; ++el) 
    {

      Elem* elem = *el;

      polariz_vec =  get_piezopolarization( elem );

    
      polariz_data[0 + el_number*3] = polariz_vec (1);
      polariz_data[1 + el_number*3] = polariz_vec (2);
      polariz_data[2 + el_number*3] = polariz_vec (3);

      el_number++;

    }

}

//---------------------------------------------------------------------------



//writes piezopolarization in GMV ot tecplot format
void Macrostrain::output_piezo(std :: string filename)
{
  
  const Mesh& mesh = equation_systems->get_mesh();

  

  std::vector<double> polariz_data;

  std::vector<std::string> polariz_names;
 
  
  prepare_polarization_data_for_output( polariz_names,  polariz_data );

  
  if (output_type == "GMV")  GMVIO_cell(mesh).write_ascii_cell_data(filename, polariz_data, polariz_names);
 
  if (output_type == "tecplot") TecplotIO_cell(mesh,false).write_cell_data(filename, polariz_data, polariz_names);
}

//---------------------------------------------------------------------------

void Macrostrain::move_nodes()
{
  const Mesh& mesh = equation_systems->get_mesh();

  //const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = *my_system;

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  // DofMap& dof_map = system.get_dof_map();

  unsigned int uvar[3] ;
  
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
  
  //--------------lattice matching correction---------------------------
  
  
  Tensor2Sym lat_matching_transformation(1);


  for (unsigned int i = 0; i <  number_of_add_var; i++)
    {      
      if (add_var[i].lat_cons)
	{
	  unsigned int index = add_var[i].index1;
	  
	  double a_new = (*solution)(add_var[i].dof_number);

	  lat_matching_transformation(index,index) +=(a_new -  substrate_lat_const[index - 1])/substrate_lat_const[index - 1];
	}
      else
	{
	  unsigned int index1 = add_var[i].index1;
	  unsigned int index2 = add_var[i].index2;
	  double shear = (*solution)(add_var[i].dof_number);

	  lat_matching_transformation(index1,index2) = shear;

	}
    }


  const unsigned int system_number = system.number();
  MeshBase::const_node_iterator       nd     = mesh.active_nodes_begin();
  const MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();

  Tensor1 r0(0); //fixed node coordinates
  Tensor1 r ;

  const Node& node_fix = mesh.node(fixed_node1);  

  for (unsigned int i = 0; i < dim; i++) r0(i + 1) = node_fix(i);



 
  for ( ; ( (nd != nd_end) ) ; ++nd)
    {
      //-----------------------------------------------------------------------------
      //
      Node* node1 = *nd;

      r = Tensor1(0);
      
      for (unsigned int i = 0; i < dim; i++)   r(i + 1) = (*node1)(i); //current node position
	  
      if (number_of_add_var !=0)
	{
	  for (unsigned int i = 0; i < dim; i++)   r(i + 1) -= u_node[node1][i]; //position without displacement

	  r = lat_matching_transformation * (r - r0) + r0; //transform cells according to the lattice matching transformation

	  for (unsigned int i = 0; i < dim; i++) r(i + 1) += u_node[node1][i]; //add back displacements
	}

      for (unsigned int i = 0; i < dim; i++) 
	{
	  const unsigned int  n_dof = node1->dof_number(system_number,uvar[i],0);
	  
	  r(i + 1) +=  (*solution)(n_dof);
	} 

      //-----------------------------------------------------
      //moove new node
       Point p;
       p(0) = r(1);
       p(1) = r(2);
       p(2) = r(3);
       
       *node1 = p;


     

    } 

 


  //--------------------------------------------------------------------

  std :: cout << "Nodes are moved. \n";

}
//-------------------------------------------------------------------------------------//

void Macrostrain::calculate_result_elem_strain_map()
{
  result_strain.clear();
  const Mesh& mesh = equation_systems->get_mesh();
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  for ( ; el != end_el ; ++el) 
    {
      const Elem* elem = *el;
      Tensor2Sym strain_result = get_strain(elem, true); //strain in CRYSTAL system
      result_strain.insert(pair <const Elem*, Tensor2Sym>  (elem, strain_result));
    }

}
//-------------------------------------------------------------------------------------//
Tensor2Sym Macrostrain::get_strain_crystal(const Elem* elem, const Point& quadratur_point )
{
  Tensor2Sym eps(0);
  map <const Elem*, Tensor2Sym> :: iterator it;
  it = result_strain.find(elem);
  //--------------------------------------------------------------------------------------------------
  //if the element elem is active for the strain simulation, it must be included in the map-----------
  if (it != result_strain.end() )
    {
      eps =   it->second;
    }
  else
    { //if the element is not included, we have to check his children or parents
      //-------------------------------------------------------------------------
      //1) may be it has a parent the belongs to the  result_strain map
     
      const Elem* el1 = elem->parent();
      
      bool out = false;
      bool found = false;

      while ( !out )
	{
	  if ( el1 != NULL )
	    {
	      it = result_strain.find(el1);
	      if ( it != result_strain.end() )
		{
		  eps = it -> second;

		  out = true;
		  found = true;
		}
	      else
		{
		  el1 = el1->parent();
		}
	    }
	  else
	    {
	      out = true;
	    }
	}
      //--------------------------------------------------------  
      if (!found)
	{//2) may be it has a child that belongs to the result_strain map
	  std::vector< const Elem * > active_children;
	  elem -> active_family_tree ( active_children, true);
	  unsigned int n = active_children.size();
      
	  for (unsigned int i1 = 0; ( i1 < n || found ); i1++)
	    {
	      it  =  result_strain.find(active_children[i1]);
	      if (it !=  result_strain.end() )
		{
		  //we have to check if this child contains a quadrature point q
		  if (active_children[i1]->contains_point(quadratur_point))
		    {
		      eps = it -> second;
		      found = true;
		    }
		  
		}
	    }
	  
      
	}
    }
  return(eps);
 
}

//-------------------------------------------------------------------------------------//
Tensor2Sym Macrostrain::get_strain(const Elem* elem, bool crystal_system )
{

  Tensor2Sym eps0;
  Tensor2Sym eps;
 

  LinearImplicitSystem& system = *my_system;

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  DofMap& dof_map = system.get_dof_map();


  unsigned int uvar[3] ;
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
  


  FEType fe_type = dof_map.variable_type(0);
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type, true)); //no scaling here
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  
  std::vector<Point> point_vec(1);


  map<const Elem*, unsigned int> :: iterator el_numb_it;

  el_numb_it = elem_numbers.find(elem);


 
  
  const unsigned int elem_number = el_numb_it->second;

 
  

  

  eps0 = eps0_of_elem[elem_number];//this is shape deformation from previous iterations

    
 

  point_vec[0]  = FEInterface::inverse_map(dim, fe_type, elem, elem->centroid());

 

  fe->reinit (elem, &point_vec);

  std::vector<unsigned int> dof_indices_component1;
  

  std::vector<unsigned int> dof_indices_component2;

  //-----------------------------------------------------------------
  //here we calculate 1/2(du/dx + du/dx) and add it to the eps0
  for (int i = 1; i <=3 ; i++)
    for (int j = 1; j <=i; j++)
      {
	dof_map.dof_indices (elem, dof_indices_component1, uvar[i-1]);
	dof_map.dof_indices (elem, dof_indices_component2, uvar[j-1]);

	const unsigned int n_u_dofs = dof_indices_component1.size(); 

	double du_i_over_dx_j = 0;
	for (unsigned int p1=0; p1<n_u_dofs; p1++)
	  {   
	    
	    if (j<=dim) du_i_over_dx_j += 0.5 *  dphi[p1][0](j-1) * (*solution)(dof_indices_component1[p1]);
	    
	    
	    if (i<=dim) du_i_over_dx_j += 0.5 *  dphi[p1][0](i-1) * (*solution)(dof_indices_component2[p1]); 
	

            //if ((j == 1) && (i == 1)) cerr <<  dof_indices_component1[0] << "\n";
	    
	  }
	    
	
	  
	eps(i,j) = eps0(i,j) + du_i_over_dx_j ; 
	
      } 


 
  //-----------------------------------------------------------------
  //we have to add lattice matching deformation
  eps += calculate_eps_lat_matching(elem);
  //------------------------------------------------------------------

  if (crystal_system)
    {//convert to crystal system

   

      ID subdomain = elem->subdomain_id();
         
      const Material* mat = _device->get_material(subdomain);

      const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

      Tensor2Gen RotM = (crystal_el->RotMatrix).transpose();//get rotation matrix
      
      Tensor2Gen eps1 = (RotM*eps)*RotM.transpose();  //transform
      
      eps  = sym(eps1); //result has to be symmetric

      assert (norm(eps - eps1) < 1e-6); //is it really symmetric

      return(eps); //return strain tensor in crystal system
    }
  else
    {
      //output in simulation system
      return(eps);
    }

}

//==============================================================================//
Tensor1 Macrostrain::get_built_in_polarization(const Elem* el, const Point& quadratur_point )
{
  //---------------calculate strain in crystal system---------------------------


  Tensor2Sym strain_cr = get_strain_crystal( el, quadratur_point);

 
  //----------------calculate polarization---------------------------------

  // std::map< unsigned int, Piezoelectricity*>::iterator piezo_it =
  //  piezo_parameters.find( material) ;

  

  ID subdomain = el->subdomain_id();
      

      
  const Material* mat = _device->get_material(subdomain);

  const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

  MacrostrainModel* macrostrain_model =  dynamic_cast<MacrostrainModel*>(   mat ->get_model(get_id())     );


  Tensor1 polariz(0);      

  if (macrostrain_model != NULL)
  {

    polariz = ( macrostrain_model->get_piezo() )-> get_polariz_cryst(strain_cr);
  


    // Tensor1 polariz = (piezo_it -> second)->get_polariz_cryst(strain_cr); //crystal system

    // std::map< unsigned int, Macrostrain::strain_param>::iterator str_it =
    //  strain_parameters.find( material) ;

    polariz =( crystal_el->RotMatrix) * polariz; //calculation system

  }
  return(polariz);

}


//-------------------------------------------------------------------------------------------/
Tensor1 Macrostrain::get_piezopolarization(const Elem* el)
{
  //---------------calculate strain in crystal system-----------------------------


   Tensor2Sym strain_cr= get_strain( el, true);

  //---------------get material number -------------------------------------------
   /*
     map<const Elem*, unsigned int> :: iterator el_numb_it;

     el_numb_it = elem_numbers.find(el);
  
     const unsigned int elem_number = el_numb_it->second;

     const unsigned int material = material_of_elem[elem_number]; //get material number

     //----------------calculate polarization-----------------------------------------

     std::map< unsigned int, Piezoelectricity*>::iterator piezo_it =
     piezo_parameters.find( material) ;

   */

   ID subdomain = el->subdomain_id();
      

      
   const Material* mat = _device->get_material(subdomain);

   const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

   MacrostrainModel* macrostrain_model =  dynamic_cast<MacrostrainModel*>(   mat ->get_model(get_id())     );
      

   Tensor1 polariz = ( macrostrain_model->get_piezo() )-> get_polariz_cryst(strain_cr); //crystal system

   

   polariz =(crystal_el->RotMatrix) * polariz; //calculation system

   return(polariz);


}
//-------------------------------------------------------------------------------------------/
void Macrostrain::define_additional_variables()
{
  number_of_add_var = 0;
  eps0_var_log = 0;

  if (!grown_on_substrate)
    {

     
      //----------------------------------------------------------------
      //lattice constants first ----------------------------------------
      if (periodicity[0]) 
	{
	  number_of_add_var++;
	  add_variable add_var1;
	  add_var1.name = "ax";
	  add_var1.lat_cons = true;
	  add_var1.index1 = 1;
	  add_var1.index2 = 1;

	  
	  add_var.push_back(add_var1);
	  eps0_var_log(1,1) = 1;
	  

	}
      
  
      if (periodicity[1] || dim < 2) 
	{
	  number_of_add_var++;
	  add_variable add_var1;
	  add_var1.name = "ay";
	  add_var1.lat_cons = true;
	  add_var1.index1 = 2;
	  add_var1.index2 = 2;

	  add_var.push_back(add_var1);
	  eps0_var_log(2,2) = 1;
	}
      
      if (periodicity[2] || dim < 3) 
	{
	  number_of_add_var++;
	  add_variable add_var1;
	  add_var1.name = "az";
	  add_var1.lat_cons = true;
	  add_var1.index1 = 3;
	  add_var1.index2 = 3;

	  add_var.push_back(add_var1);
	  eps0_var_log(3,3) = 1;
	}
      //---------------------------------------------------------------------------//
      //---shears------------------------------------------------------------------//
/*
      if (dim == 1)
	{
	  number_of_add_var++;
	  add_variable add_var1;
	  add_var1.name = "eps_xy";
	  add_var1.lat_cons = false;
	  add_var1.index1 = 2;
	  add_var1.index2 = 1;
	  add_var.push_back(add_var1);
	  eps0_var_log(2,1) = 1;
	}
*/
      

      if ( (eps0_var_log(1,1) == 1) && ( eps0_var_log(2,2) == 1 ) ) 
	{
	  number_of_add_var++;
	  add_variable add_var1;
	  add_var1.name = "eps_xy";
	  add_var1.lat_cons = false;
	  add_var1.index1 = 2;
	  add_var1.index2 = 1;
	  add_var.push_back(add_var1);
	  eps0_var_log(2,1) = 1;
	}
      
      
      if ( (eps0_var_log(1,1) == 1) && ( eps0_var_log(3,3) == 1 ) ) 
	{
	  number_of_add_var++;
	  add_variable add_var1;
	  add_var1.name = "eps_xz";
	  add_var1.lat_cons = false;
	  add_var1.index1 = 3;
	  add_var1.index2 = 1;

	  add_var.push_back(add_var1);
	  eps0_var_log(3,1) = 1;
	}
      
      if ( (eps0_var_log(2,2) == 1) && ( eps0_var_log(3,3) == 1 ) ) 
	{
	  number_of_add_var++;
	  add_variable add_var1;
	  add_var1.name = "eps_yz";
	  add_var1.lat_cons = false;
	  add_var1.index1 = 3;
	  add_var1.index2 = 2;
	  add_var.push_back(add_var1);
	  eps0_var_log(3,2) = 1;
	}

	
    }


 
  
}


  //---------------------------------------------------------------------------//


void  Macrostrain::set_up_additional_dofs()
{


  const Mesh& mesh =  equation_systems->get_mesh();

  LinearImplicitSystem& system = *my_system;

  DofMap& dof_map = system.get_dof_map();

  std::vector<unsigned int> dof_indices_component;

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  unsigned int el_number = 0;

  unsigned int var_fict ;

  if ( add_var.size() != 0) var_fict = system.variable_number("fict");


  for ( ; (el != end_el && el_number <  add_var.size() )  ; ++el  ) 
    {
      Elem* elem = *el;

      add_var[el_number].element = elem;
      
      dof_map.dof_indices (elem, dof_indices_component, var_fict);

      add_var[el_number].dof_number = dof_indices_component[0];

      el_number++; 
    }

  add_dofs_vector.clear();
  for (unsigned int i = 0; i <  add_var.size(); i++)
    {
      add_dofs_vector.push_back( add_var[i].dof_number )   ;
    } 


}
//-------------------------------------------------------------------------------------------/

Tensor2Sym Macrostrain::calculate_eps_lat_matching(const Elem* elem)
{
  //constant part of the lattice matching tensor

  ID subdomain = elem->subdomain_id();
      

      
  const Material* mat = _device->get_material(subdomain);

  const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());
  

  Tensor2Sym eps0 = crystal_el->get_const_eps0(substrate_lat_const, eps0_var_log);
  //------
  //variable part:
  for (unsigned int i = 0; i < number_of_add_var; i++)
    {
      unsigned int dof_number = add_dofs_vector[i];

      double coeff = ( *(my_system->solution) )( dof_number ); 

      eps0 += coeff * crystal_el->get_var_eps0( add_var[i].name );
    } 

  return(eps0);

}
//-------------------------------------------------------------------------------------------/
void Macrostrain::init_substrate()
{
  substrate_shear = Tensor2Sym(0);

 
  
  substrate_crystal->get_lat_const(substrate_lat_const);
}
//--------------------------------------------------------------------------------------------/

void Macrostrain::update_substrate()
{

  LinearImplicitSystem& system = *my_system;

  AutoPtr<NumericVector<Number> >& solution = system.solution;


  for (unsigned int i = 0; i <  number_of_add_var; i++)
    {
      
      if (add_var[i].lat_cons)
	{
	  unsigned int index = add_var[i].index1;
	  
	  substrate_lat_const[index - 1] = (*solution)(add_var[i].dof_number);
	

	}
      else
	{
	  unsigned int index1 = add_var[i].index1;
	  unsigned int index2 = add_var[i].index2;
	  double shear = (*solution)(add_var[i].dof_number);

	  substrate_shear(index1,index2) = shear;

	}
    }

}
//--------------------------------------------------------------------------------------------/
void Macrostrain::init_u_node()
{
  std :: vector <double> single_node(3);

  for (unsigned int i = 0 ; i < 3 ; i++) single_node[i] = 0.0;
 
  u_node.clear();
  
  const Mesh& mesh = equation_systems->get_mesh();
  MeshBase::const_element_iterator el  = mesh.active_elements_begin();
  MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el ; ++el) 
  { 
    const Elem* elem = *el;

    const unsigned int num_nodes = elem->n_nodes();

    for (short i = 0; i < num_nodes; i++)
    {
      const Node* nd = elem->get_node(i);
      if (u_node.find(nd) == u_node.end())
	u_node.insert( pair<const Node*,vector<double> > (nd,single_node)  ); 
    }
     
  }
  


}
//-------------------------------------------------------------------------------------------/
void Macrostrain::update_u_node()
{
  const Mesh& mesh = equation_systems->get_mesh();

  LinearImplicitSystem& system = *my_system;

  AutoPtr<NumericVector<Number> >& solution = system.solution;
 
  const unsigned int system_number = system.number();

  unsigned int uvar[3] ;
  
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }


 

  map<const Node*, vector<double> > temp;

  MeshBase::const_element_iterator el  = mesh.active_elements_begin();
  MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el ; ++el) 
  { 
    const Elem* elem = *el;

    const unsigned int num_nodes = elem->n_nodes();

    for (short i = 0; i < num_nodes; i++)
    {
      const Node* nd = elem->get_node(i);

      if (temp.find(nd) == temp.end())
      {
	vector <double> du(3);
	for (unsigned int i = 0; i < 3; i++) //<3 , not < dim (necessary for atoms!) 
	{
	  const unsigned int  n_dof = nd->dof_number(system_number,uvar[i],0);
	  du[i] = (*solution)(n_dof);
	}
	temp.insert( pair<const Node*,vector<double> > (nd,du)  ); 
      }

     
    }
     
  }


  map<const Node*, vector<double> >::iterator it = u_node.begin();
  map<const Node*, vector<double> >::iterator it_end = u_node.end();

  for( ; it != it_end ; )
  {
    const Node* nd = it->first;

    for (unsigned int i = 0; i < 3; i++) (it->second)[i] +=temp[nd][i]; 
  }


 

}
//-------------------------------------------------------------------------------------------/
void Macrostrain::output_add_strain_variables(string filename)
{
  std::ofstream out;

  if (number_of_add_var !=0) 
    {
      LinearImplicitSystem& system = *my_system;

      AutoPtr<NumericVector<Number> >& solution = system.solution;
      

      std::ofstream out (filename.c_str());

      assert (out.good());

      //std::cerr << " number_of_add_var  = " << number_of_add_var << "\n"; 


      for (unsigned int i = 0 ; i <  number_of_add_var; i++ )
	{
	  out << add_var[i].name  << "            "  << (*solution)(add_var[i].dof_number) << "\n";
	  
	  
	  
	}

    }

}
//-------------------------------------------------------------------------------------------/
unsigned int Macrostrain::get_number_of_the_fixed_node(Point point)
{
  const Mesh& mesh =  equation_systems->get_mesh();

  Elem*  elem1;
  bool  found;
  //-------------------------------------
  // the most coarse elements first
  unsigned int refinement_level = 0; 
  MeshBase::const_element_iterator el3  = mesh.level_elements_begin(refinement_level);
  MeshBase::const_element_iterator end_el3 = mesh.level_elements_end(refinement_level);
  
  
  found = false;

  //-------------------------------------  
  for ( ; ( (el3 != end_el3) ) ; ++el3)  
    {
      Elem* elem = *el3;
    
      if (elem->contains_point(point))
	{
	  elem1 = elem;
	  found = true;
	  break;
	  
	}
      
    }
  //----------------------------------------
  assert(found);
  //children of the  most coarse element 
  while ( !( elem1->active() ) )
    {
      
      for (unsigned int i=0 ; i < elem1->n_children() ; i++)
	{
	  Elem* 	child = elem1->child(i);
	  
	  if (child->contains_point(point))
	    {
	      elem1 = child;
	      break;
	    }
	    
	}
    } 

  
  // elem1 contains node
  vector <double> distances;
  for (unsigned int i = 0; i < elem1->n_nodes(); i++)
    {
      Point point1 = elem1->point(i);

      double d = std::sqrt( (point1(0) - point(0))*(point1(0) - point(0)) + (point1(1) - point(1)) *(point1(1) - point(1))  + 
			    (point1(2) - point(2))*(point1(2) - point(2)) );

      distances.push_back(d);

    }

  vector<double>::const_iterator it = min_element(distances.begin(), distances.end());

  //----------------------------------------------------------------  

  for (unsigned int i = 0; i < elem1->n_nodes(); i++)
    {
      if (*it == distances[i]) 
	{
	  return( elem1->node(i) );
	}
    }


}

//-------------------------------------------------------------------------------------------/
void Macrostrain::read_atom_structure(const std::string filename)
{
  string  string_from_file;

  std::ifstream atoms_file;

  atoms_file.open(filename.c_str());
  //----------------------------------------------------------------


  //----------------------------------------------------------------

  if (!atoms_file.good())
    {
      cerr << "Error: file with atom coordinates  "<< filename << "   can not be opened";
      error();	
    }

  //-----------------------------------------------------------
  //determination of number of atoms 
  unsigned int N_atoms = 0;

  while (getline(atoms_file, string_from_file ))  
    {
     
      N_atoms++;     
    }

  
  atoms_file.close(); 
  //------------------------------------------------------------
  atom  atom_from_file;

  atom_from_file.mat_number = 0;
  atom_from_file.type = 0;
  atom_from_file.relative_point = Point(0.0, 0.0, 0.0);
  atom_from_file.element = NULL; 
  
  atom_structure.resize(N_atoms,atom_from_file);

  cout << "Number of atoms in " << filename << "  file  " <<  N_atoms;
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  //---reading of atoms and processing---------------------------------


  //--------mesh related objects----------------------------------------
  const Mesh& mesh =  equation_systems->get_mesh();

  LinearImplicitSystem& system = *my_system;

  DofMap& dof_map = system.get_dof_map();

  FEType fe_type = dof_map.variable_type(0);
  




  //--------------------------------------------------------------------
  std::ifstream atoms_file1(filename.c_str());
 
  for (unsigned int i = 0; i < N_atoms; i++)
    {//loop over atoms
      getline(atoms_file1, string_from_file );
      
     

      istringstream input_string(string_from_file);
      
      unsigned int n ;
      int t ;
      int mat;
      double x;
      double y;
      double z;
      input_string >>mat >> t >> x >> y >> z ;
       
      vector<double> coordinate(3);
      coordinate[0] = x; coordinate[1] = y; coordinate[2] = z;
      
    

      //   cerr << x << "    "<< y <<"   "<< z <<"\n";
      


      atom_structure[i].mat_number = mat;
      atom_structure[i].type  = t;
      
      //-----------------------------------------------------------
      //determination if the atom belongs to the similation domain
      //-----------------------------------------------------------
      Point point2;

      for (unsigned int i1 = 0; i1 < dim; i1++)
	{
	  point2(i1) = coordinate[i1];
	
	}

      //-------------------------------------------------------------
      //find element that contains the point 

      unsigned int refinement_level = 0; 
      MeshBase::const_element_iterator el3  = mesh.level_elements_begin(refinement_level);
      MeshBase::const_element_iterator end_el3 = mesh.level_elements_end(refinement_level);
		      
      Elem*  elem1;
      
      bool   found = false;

      for ( ; ( (el3 != end_el3) ) ; ++el3)  
	{
	  Elem* elem = *el3;
	 
	  if  (may_belong_to_element(elem,  point2))
	    {
	      if (elem->contains_point(point2))
		{
		  elem1 = elem;
		  found = true;
		  break;
		}
	      
	    }
	    
	}
		      
      if (!found) 
	{
	  //atom is not found and will be ignored
	  cerr << "WARNING: atom does not belong to the macroscopic domain\n";
	  cerr << "The atom number  " << i << "  will be ignored\n";  
	  cerr << x <<"   "<< y <<"   " << z <<"   "<<  mat <<"   "<< t <<"\n"; 
	    
	}
      else
	{ //atom is found and will be processed
	  //children of the  most coarse element 
	  while ( !( elem1->active() ) )
	    {
	      
	      for (unsigned int i=0 ; i < elem1->n_children() ; i++)
		{
		  Elem* 	child = elem1->child(i);
		  if  (may_belong_to_element(elem1,  point2))
		    {
		      if (child->contains_point(point2))
			{
			  elem1 = child;
			  break;
			}
		    }
		  
		}
	    }
		      
	
	  atom_structure[i].relative_point = FEInterface::inverse_map(dim, fe_type, elem1, point2);
	 
	  atom_structure[i].element = elem1; 
	  //-------------------------------------------------------------------
	}
	  
	  
    }

  atoms_file1.close(); 
  

}

//-------------------------------------------------------------------------------------------/
void  Macrostrain::write_atom_potential()
{

 

  std::ofstream potential_file;

  if (poisson_equation != NULL)
  {

    ID pot_ID = poisson_equation->get_variable_id("ElPotential");


    potential_file.open( atom_potential_filename.c_str() );
    if (!potential_file.good())
    {
      cerr << "Error: file with atom potentials can not be opened\n";
      cerr <<  atom_potential_filename.c_str() << "\n";
      error();	
    }

    unsigned int Number_of_atom = atom_structure.size();
    vector<Point> point_vec(1);
    double potential_value;

    for (unsigned int i = 0; i < Number_of_atom ; i++)
    {//atoms loop
      if (atom_structure[i].element != NULL)
      {

	point_vec[0] =  atom_structure[i].relative_point ;
	vector<double> values;
	poisson_equation->get_solution(atom_structure[i].element, point_vec, pot_ID, values);
	potential_value = values[0];
      }

      potential_file << setw(20) <<   setprecision(12) << potential_value << "\n";


    }

  }



}



//-------------------------------------------------------------------------------------------/
void  Macrostrain::write_atom_displacements(const std::string filename)
{//-------------------------------------------------------------------
 


 //file opening
  string  string_from_file;

  std::ifstream atom_in_file;
  std::ofstream displacement_file;

  
 
  

  displacement_file.open(filename.c_str());

  if (!displacement_file.good())
  {
    cerr << "Error: file with atom displacements can not be opened\n";
    cerr << filename.c_str() << "\n";
    error();	
  }


 

  //--------------------------------------------------------------------
 if (atom_output_type=="uptight")
   {
     //for uptight we have output displacements.
     //therefore we need also the initial positions
    
     atom_in_file.open(atom_structure_filename.c_str());

     if ( ! atom_in_file.good () )
       {
	 cerr << "Error: file with atom positions can not be opened\n";
	 cerr << atom_structure_filename.c_str() << "\n";
	 error();
       }
    
   }
  //-------------------------------------------------------------------
  const Mesh& mesh =  equation_systems->get_mesh();


  LinearImplicitSystem& system = *my_system;

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  DofMap& dof_map = system.get_dof_map();

  std::vector<unsigned int> dof_indices_component1;

  unsigned int uvar[3] ;
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
 
  

  FEType fe_type = dof_map.variable_type(0);
 
  
  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type, true)); //no scaling here

 
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  const std::vector<Point>& q_point = fe->get_xyz();

  vector<Point> point_vec(1);

  //----------------------------------------------------------------------
  double substrate_lat_const_initial[3];

 
  
  substrate_crystal->get_lat_const(substrate_lat_const_initial);
  //----------------------------------------------------------------------
  

  unsigned int Number_of_atom = atom_structure.size();

  for (unsigned int i = 0; i < Number_of_atom ; i++)
    {//atoms loop
      if (atom_structure[i].element != NULL)
	{//atom belongs to the simulation domain
	 
	  vector <double> new_pos_of_atom(3,0.0);
	  

	  

	  point_vec[0] =  atom_structure[i].relative_point ;
	  
	  


	  fe->reinit(atom_structure[i].element, &point_vec);


	  /*----------------------------------------------
	    first we moove atom because the grid moves 
	    in 3D this is enough to get into account the external strain
	    ---------------------------------------------*/
	  for (short coord = 0; coord < dim; coord++)
	    new_pos_of_atom[coord] = q_point[0](coord);

	 /*----------------------------------------------
	   in 2D and 1D there are displacements in the direction perpendcular to the simulation space
	   if the reference lattice is fixed, this is enough
	   ------------------------------------------------*/
	  if (dim < 3)
	    {//dim = 1,2
	      vector <double> u_vector(3,0.0);	      
	      for (unsigned int nd=0; nd < atom_structure[i].element->n_nodes(); nd++)
		{
		  for (short coord = dim; coord < 3; coord++)
		    {
		     
 
		      u_vector[coord] += u_node[atom_structure[i].element->get_node(nd)][coord] * phi[nd][0]; 
		    }
		}

	      
	      for (short i1 = 1; i1 < 3; i1++)  new_pos_of_atom[i1] += u_vector[i];

	      /*------------
		If the reference lattice in a parallel space may change
		then we need additional things 
	      ---------------*/
	      if (!grown_on_substrate) /*otherwise the parallel space is fixed*/
		{ //free standing
		  //----------------------------------------------------------------------------------
		  //creation of the transformation tensor
		  Tensor2Sym substr_deform(1);
		 
		  for (short i1 = dim + 1; i1 <= 3; i1++)
		    {
		      substr_deform(i1,i1)  += ( (substrate_lat_const[i1 - 1] - substrate_lat_const_initial[i1 - 1])/
			 substrate_lat_const_initial[i1 - 1] ) ;
		      for (short j1 = dim + 1; j1 < i1; j1++)
			{
			  substr_deform(i1,j1) +=   substrate_shear(i1,j1);
			}	

		    } 
		  //--------------------------------------------------------------------------------------
		  //application of the transformation
		  Tensor1 pos_vector_math;
		  for (short i1 =1; i1 <= 3; i1++)		    pos_vector_math(i1) =  new_pos_of_atom[i1 - 1];
		  
		  pos_vector_math = substr_deform*pos_vector_math;

		  for (short i1 =1; i1 <= 3; i1++)		   new_pos_of_atom[i1 - 1] =  pos_vector_math(i1);
		  //---------------------------------------------------------------------------------------

		}
	      
	     


	    }

	  //-------------------------------------------------------

	

	 
	  //-------------------------------------------------------------------
	  //output of new coordinates
	  
	  if (atom_output_type=="povray")
	    {
	      displacement_file <<  setw(20) <<  atom_structure[i].mat_number << ",";
	      displacement_file <<  setw(20) <<  atom_structure[i].type << ",";
	  
	  
	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[0]<< ","  ;
	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[1]<< ","  ;
	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[2]<< ","  ;

	      if (output_strain_on_atoms)
		{
		  Tensor2Sym epsilon = get_strain(atom_structure[i].element);
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(1,1) << "," ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(2,2) << "," ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,3) << "," ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(2,1) << "," ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,1) << "," ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,2) << "," ;
		}
	    }


	  if (atom_output_type=="uptight")
	    {


	      getline(atom_in_file, string_from_file );

	   
	      istringstream input_string(string_from_file);
	      

	      int t ;
	      int mat;
	      double x;
	      double y;
	      double z;

	      
	      
	      input_string >>mat >> t >> x >> y >> z ;
	  
	    

	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[0] - x << ","  ;
	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[1] - y << ","  ;
	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[2] - z << ","  ;

	      if (output_strain_on_atoms)
		{
		  Tensor2Sym epsilon = get_strain(atom_structure[i].element);
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(1,1)  ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(2,2)  ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,3)  ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(2,1)  ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,1)  ;
		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,2)  ;
		}


	  

	    }

	  displacement_file <<  '\n';
	  
	
	}//


    }//end of atoms loop
    
  
    



}
//-------------------------------------------------------------------------------------




//-------------------------------------------------------------------------------------------/
bool Macrostrain::may_belong_to_element(const Elem* element, Point& point)
{
  const unsigned int n = element->n_nodes();
  double  min_x;double  min_y; double  min_z;
  double  max_x;double  max_y; double  max_z;
  
  Point vertex = element->point(0);
  min_x = vertex(0); min_y = vertex(1); min_z = vertex(2);
  max_x = min_x;  max_y = min_y; max_z = min_z;

  for (unsigned int i = 1 ; i < n ; i++)
    {
      vertex = element->point(i);
      double x = vertex(0);
      double y = vertex(1);
      double z = vertex(2);

      if (min_x > x) min_x = x; if (min_y > y) min_y = y; if (min_z > z) min_z = z;

      if (max_x < x) max_x = x; if (max_y < y) max_y = y; if (max_z < z) max_z = z;

      
      
    }

  if ( (point(0) > max_x) ||  (point(0) < min_x) ||
       (point(1) > max_y) ||  (point(1) < min_y) ||
       (point(2) > max_z) ||  (point(2) < min_z) ) 

    {    return(false);}

  else
    {    return(true) ; }

    

}


//-------------------------------------------------------------------------------------------/

//-------------------------------------------------------------------------------------------/
unsigned int Macrostrain::find_nearest_node(Point& point)
{
  //finds a node number nearest to the point
  const Mesh& mesh =  equation_systems->get_mesh();
  unsigned int num_nodes = mesh.n_nodes();
  double distance;

  const Node& nd = mesh.node(0);
  distance =        (nd(0) - point(0))*(nd(0) - point(0)) +
    (nd(1) - point(1))*(nd(1) - point(1)) +
    (nd(2) - point(2))*(nd(2) - point(2));
			
  unsigned int result = 0;
	      
  for (unsigned int i = 1 ; i < num_nodes; i++)
    {
      const Node& nd = mesh.node(i);
      if (nd.active()) 
	{
	  const double distance1 =        (nd(0) - point(0))*(nd(0) - point(0))+
	    (nd(1) - point(1))*(nd(1) - point(1)) +
	    (nd(2) - point(2))*(nd(2) - point(2));
	  
	  if (distance1 < distance)
	    {
	      result = i;
	      distance = distance1;
	    }

	}
    }

  return(result);

}
//-------------------------------------------------------------------------------------------/
double Macrostrain::norm_of_difference(NumericVector<Number>& solution1, NumericVector<Number>& solution2)
{

  //the idea is to exclude from comparison the additional variables 

  MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

  DofMap& dof_map = my_system->get_dof_map();
  vector<unsigned int> dof_indices;

  double  norm = 0;

  for ( ; el != end_el ; ++el) 
    {//el

      const Elem* elem = *el;

      for (short alpha = 0; alpha < 3; alpha++)
	{
	  dof_map.dof_indices (elem, dof_indices, alpha);

	  short n = dof_indices.size();

	  for (short i = 0; i < n; i++)
	    {
	      double t = abs( solution1(dof_indices[i]) - solution2(dof_indices[i]) );

	      if (t > norm) norm = t;
	    }

	}

    }

  return norm;

}

//-------------------------------------------------------------------------------------------/
Macrostrain::~Macrostrain()
{

 

  //equation_systems->delete_system(system_name);
}


//-------------------------------------------------------------------------------------------//


Macrostrain::Macrostrain(void )
{
  poisson_equation = NULL;

  my_solver = NULL;
}
