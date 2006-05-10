#include "macrostrain.h"

//-----------------------------------------------------------------//
//Static objects must be "mentioned" here! 

std::string    Macrostrain:: uname_vec[3];
unsigned int   Macrostrain:: dim;
bool           Macrostrain:: grown_on_substrate;

unsigned  int  Macrostrain:: substr_mat;


Tensor2Sym     Macrostrain:: eps0_var_log;
unsigned int   Macrostrain:: number_of_add_var_static;

std::vector<stiffness>*       Macrostrain:: C_tensor_temp;
std::vector<rotated_crystal>* Macrostrain:: crystal_temp;
std::vector <int>*            Macrostrain:: material_of_elem_temp;
std::vector <Tensor2Sym>*     Macrostrain:: eps0_of_elem_temp;
std::vector<Macrostrain::add_variable>*    Macrostrain:: add_var_temp; 
std::vector<unsigned int>     Macrostrain:: add_dofs_vector;
std::vector<unsigned int>*    Macrostrain:: zero_set_dofs_temp;
std::map<const Elem*, map <unsigned int, double> >* Macrostrain::boundary_cond_elem_temp;

std::set <unsigned int>* Macrostrain::substrate_nodes_temp;

double Macrostrain::substrate_lat_const[3];

Tensor2Sym Macrostrain::substrate_shear;

unsigned int Macrostrain::fixed_node1_temp;
unsigned int Macrostrain::fixed_node2_temp;
unsigned int Macrostrain::fixed_node3_temp;


//-----------------------------------------------------------------//
Macrostrain::Macrostrain(const options& opt,   Mesh&  mesh )
{
  
 
    max_r_steps = opt.max_r_steps;
    uniform_refinement = opt.uniform_refinement;
    refine_fraction = opt.refine_fraction;
    coarsen_fraction = opt.coarsen_fraction;
    max_ref_level = opt.max_ref_level;
    tolerance  = opt.tolerance  ;
    max_shape_steps = opt.max_shape_steps;
 
    mesh_input_file = opt.mesh_input_file;
    grown_on_substrate = opt.grown_on_substrate;
    calculate_atom_displacements = opt.calculate_atom_displacements;
    atom_structure_filename = opt.atom_structure_filename;
    atom_displacements_filename = opt.atom_displacements_filename;
  

    for (unsigned int i =0; i <=2; i ++ ) periodicity[i] = opt.periodicity[i];
    substr_mat = opt.substr_mat;

      
    dim = mesh.mesh_dimension();
    uname_vec[0]="ux";
    uname_vec[1]="uy";
    uname_vec[2]="uz";

    //define additional parameters-------------------------------------------

    define_additional_variables();


    // Create an equation systems object.
    equation_systems = new EquationSystems(mesh);
    
    // Declare the Poisson system and its variables.
    // The Poisson system is another example of a steady system.
    equation_systems->add_system<LinearImplicitSystem> ("Strain");


    
      
    equation_systems->parameters.set<Real>("linear solver tolerance") = tolerance; 
    // Adds the variable "u" to "Poisson".  "u"
    // will be approximated using second-order approximation.
  
    dim = mesh.mesh_dimension();

    //---------------------------------------------------------------------------------------
    //add normal variables
	
    for (unsigned int i = 0; i <  3 ; i++)  
      {  
	equation_systems->get_system("Strain").add_variable(uname_vec[i], FIRST);
      }
    
 
    //---------------------------------------------------------------------------------------
    //add aditional varables 


    if (number_of_add_var != 0)
      {
	FEType fe_type(CONSTANT,MONOMIAL);
	equation_systems->get_system("Strain").add_variable("fict", fe_type);
      }
    //-----------------------------------------------------------------------------------
      

    // Give the system a pointer to the matrix assembly
    // function.  This will be called when needed by the
    // library.

   



    equation_systems->get_system("Strain").attach_assemble_function (assemble_strain_matrix);
  
    //---------------------------------------------------------------------------------------
 
    for (unsigned int i=0 ; i < 3; i++) 
      {
	fixed_point1(i) = opt.fixed_point1[i];
	fixed_point2(i) = opt.fixed_point2[i];
	fixed_point3(i) = opt.fixed_point3[i];
      }

    //-------------------------------------------------------------------//

    //define  max and min coordinates 
    
    unsigned int num_nodes = mesh.n_nodes();
    const Node& nd = mesh.node(0);
    for (unsigned i = 0; i < 3; i++)
      {
	min_coord[i] = nd(i);
	max_coord[i] = nd(i);
      }

    for (unsigned i = 1; i < num_nodes; i++)
      {
	const Node& nd = mesh.node(i);
	for (unsigned i = 0; i < 3; i++)
	  {
	    if (min_coord[i] < nd(i)) min_coord[i] = nd(i);
	    if (max_coord[i] > nd(i)) max_coord[i] = nd(i);
	  
	  }

      }
   

  //-------------------------------------------------------------------//

  


    intermediate_output = opt.intermediate_output;

    output_strain_on_atoms = false;
    atom_output_type = "uptight";

    //------------------------------------------------------------------//
    if ((opt.output_type == "GMV")|| (opt.output_type == "tecplot"))
      {
	output_type = opt.output_type;
      }
    else
      {
	cerr << "Strange output type: " << opt.output_type << "\n";
	cerr << "We set it to GMV\n";
	output_type = "GMV";
    }
    
                     
}
//--------------------------------------------------------------------//

void Macrostrain::define_fixed_nodes()
{
  fixed_node1 = find_nearest_node(fixed_point1);
  fixed_node2 = find_nearest_node(fixed_point2);
  fixed_node3 = find_nearest_node(fixed_point3);
}

//-------------------------------------------------------------------//
void Macrostrain::define_substrate_bc(unsigned int substrate_bc_number_in)
{
  substrate_bc_number = substrate_bc_number_in;
}

//--------------------------------------------------------------------//
void Macrostrain::define_stress_value (const map <unsigned int, double> & stress_map_in)
{
  stress_values = stress_map_in;
}

//-------------------------------------------------------------------//
void Macrostrain::define_BC_map (const map <unsigned int , vector<unsigned int> > & bc_cond_in  )
{
  boundary_cond_nodes = bc_cond_in;
}
 


//-------------------------------------------------------------------//
void Macrostrain::assign_mesh_data(MeshData& mesh_data_in)
{
   meshdata = &mesh_data_in;
}


//-----------------------------------------------------------------//
void Macrostrain::define_strain_parameters(const std::vector<stiffness>&        C_tensor_in,
				const std::vector<rotated_crystal>&  crystal_in)
{
  C_tensor = C_tensor_in;
  crystal = crystal_in;
}


//------------------------------------------------------------------//


//------------------------------------------------------------------//

void Macrostrain::define_piezo_moduli(std::vector<Piezoelectricity>&  piezo_in)
{
  piezo = piezo_in; 
}
//------------------------------------------------------------------//
void Macrostrain::create_substate_nodes_set()
{
  if (grown_on_substrate)
    {
      //--------------------------------------------------------------
      //first we create a list of substrate nodes
      map <unsigned int , vector<unsigned int> >   :: iterator bc_pos;
  
      bc_pos =  boundary_cond_nodes.find(substrate_bc_number);
      if (bc_pos == boundary_cond_nodes.end()) 
	{
	  cerr << "Can not find boundary condition for substrate #  "<< substrate_bc_number  <<"  \n";
	  exit(1);
	}
  
      vector<unsigned int> sub_nodes = bc_pos->second;
      unsigned int num_nodes = sub_nodes.size();

      substrate_nodes.clear();
      
      for (unsigned int i = 0; i <  num_nodes; i++)
	{
	  substrate_nodes.insert(sub_nodes[i]);
	}
      
      //-------------------------------------------------------------
      // we create substrate_faces (2D/3D only)

      if (dim > 1)
	{
	  set<unsigned int>:: iterator  node_begin = substrate_nodes.begin();
	  set<unsigned int>:: iterator  node_end = substrate_nodes.end();

	  substrate_faces.clear();

	  const Mesh& mesh = equation_systems->get_mesh();

	  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
	  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
	  
	  for ( ; el != end_el ; ++el) 
	    {
	      
	      const Elem* elem = *el;
	      
	      
	      unsigned int n_sides = elem -> n_sides();
	      set<unsigned int> s_faces;
 

	      for (unsigned int s = 0; s < n_sides; s++)
		{
		  
		  
		  AutoPtr<Elem> side = elem->build_side(s);
		  bool found = true;
		  for (int i = 0; i < side->n_nodes(); i++)
		    {
		      if (substrate_nodes.find(side->node(i)) == node_end)       found = false;
		    }
		  if (found)
		    {
		      s_faces.insert(s);
		    }
		}

	    

	 

	      if (!(s_faces.empty() )) substrate_faces.insert(pair < const Elem*, set<unsigned int> > (elem, s_faces ))  ;

	    }
	  
	}
    }
}
//------------------------------------------------------------------//
void Macrostrain::update_substrate_nodes_set()
{
  if (grown_on_substrate)
    {//grown on substrate
      if (dim > 1) //only in this case we have to update
	{
	  map<const Elem*, std::set <unsigned int > >  new_map;
	  new_map.clear();

	  substrate_nodes.clear();
	  const Mesh& mesh = equation_systems->get_mesh();

	  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
	  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
      
	  for ( ; el != end_el ; ++el) 
	    { //--el_loop
	      const Elem* elem = *el;
	      //---does this element belong to the map?----
	      map<const Elem*, std::set <unsigned int > > :: iterator it;
	      it = substrate_faces.find(elem);
	      if (it != substrate_faces.end())
		{//----yes, it belongs to the map
		  set<unsigned int> el_sides = it -> second;
		  
		  new_map.insert(pair<const Elem*, std::set <unsigned int > > (elem, el_sides)); //copy it to the new map
		  set<unsigned int> :: iterator  side_it = el_sides.begin();
		  for ( ; side_it != el_sides.end() ; ++side_it)
		    {
		      unsigned int side_number = *side_it;
		      AutoPtr<Elem> side = elem->build_side(side_number);
		      for (int i = 0; i < side->n_nodes(); i++) substrate_nodes.insert(side->node(i));

		   
		      
		    }
		}
	      else 
		{//---no, it does not belong to the map
		  //let us check its parent 
		  const Elem* elem_parent = elem->parent();
		  map<const Elem*, std::set <unsigned int > > :: iterator it;
		  it = substrate_faces.find(elem_parent);
		  if (it !=  substrate_faces.end())
		    { //top parent belongs to a map
		      set <unsigned int> new_sides;
		      std::set <unsigned int > parent_el_sides;
		      parent_el_sides = it -> second;
		      //now we check if a child lies on a necessary side
		      
		      set <unsigned int, double> :: iterator par_it  =  parent_el_sides.begin();
		      set <unsigned int, double> :: iterator par_end =  parent_el_sides.end();
		      for ( ; par_it != par_end; ++par_it) 
			{
			  unsigned int side_number = *par_it;
			  /* I assume that side number i of an element and its parent are parallel faces */
			  /* I hope this is correct */
		       
			  if ((elem->neighbor(side_number) == NULL)) // on boundary
			    {
			      AutoPtr<Elem> side = elem->build_side(side_number);
			      new_sides.insert(side_number);
			      for (int i = 0; i < side->n_nodes(); i++) substrate_nodes.insert(side->node(i));


			    }
			  
			}
		      if (!(new_sides.empty()))  new_map.insert(pair< const Elem*, set <unsigned int>  > (elem, new_sides));
		      
		    }
		}
	      
	      
	    } 
	  
	  
	  substrate_faces = new_map;
	}
    }

}

//-----------------------------------------------------------------//
void Macrostrain::create_bondary_conditions_map()
{
  // Get a constant reference to the mesh object.
  const Mesh& mesh = equation_systems->get_mesh();

  boundary_cond_elem.clear();

 
  map <unsigned int, double>::iterator   stress_bc;
  map <unsigned int, double>::iterator   stress_bc_end = stress_values.end();
  
  stress_bc = stress_values.begin();
  
  for (       ;  stress_bc != stress_bc_end ; ++stress_bc)
    { //boundary conditions loop
      //----------------------------------------------------------------
      //we have to find if there is information in a mesh file for this stress applied

      unsigned int stress_number = stress_bc->first;

      map <unsigned int , vector<unsigned int> >   :: iterator bc_pos;

      bc_pos =  boundary_cond_nodes.find(stress_number);

      if (bc_pos == boundary_cond_nodes.end()) 
	{
	  cerr << "Can not find boundary condition #  "<< stress_number <<"  \n";
	  exit(1);
	}
      //---------------------------------------------------------------
      // there is such information, we can find elements and sides
      const vector<unsigned int> nodes_of_bc = bc_pos->second;
      
      vector<unsigned int> :: const_iterator n_it;
      const vector<unsigned int> :: const_iterator n_begin = nodes_of_bc.begin();
      const vector<unsigned int> :: const_iterator n_end   = nodes_of_bc.end();

      MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

      for ( ; el != end_el ; ++el) 
	{
	
	  const Elem* elem = *el;
	  map <unsigned int, double>  element_map;
	  
	  if (dim > 1)
	    {
	      unsigned int n_sides = elem -> n_sides();
	      


	      for (unsigned int s = 0; s < n_sides; s++)
		{
		  if (elem->neighbor(s) == NULL)
		    {
		      bool found = true;
		      AutoPtr<Elem> side = elem->build_side(s);
		      for (int i = 0; i < side->n_nodes(); i++)
			{
			  if (find(n_begin, n_end,side->node(i) ) == n_end)
			    { found = false;}		    
			}
		      if (found) 
			{
			  
			  element_map.insert(pair<unsigned int, double>(s,stress_bc->second) );
			}
		    }
		} 
	    }
	  else
	    {
	      for (int i = 0; i < 2; i++)
		{
		   if (find(n_begin, n_end, elem->node(i) ) != n_end) 
		     {
		       element_map.insert(pair<unsigned int, double>(i,stress_bc->second) );
		     }
		}
	    }

	  if (!(element_map.empty())) 
	    {
	     
	      boundary_cond_elem.insert(pair<const Elem*, map<unsigned int, double> > (elem,element_map));
	    
	    }
	 
	}
    }
  
 



 
}
//-----------------------------------------------------------------//
void Macrostrain::update_bondary_conditions_map()
{

  if (dim > 1)
    {
      map <const Elem*, map <unsigned int, double> > new_map;

      new_map.clear();

      const Mesh& mesh = equation_systems->get_mesh();

      MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
      for ( ; el != end_el ; ++el) 
	{ //--el_loop
	  const Elem* elem = *el;
	  if (elem -> on_boundary()) 
	    {//element is on boundary
	      //---does this element belong to the map?----
	      map<const Elem*, std::map <unsigned int , double > > :: iterator it_bc;
	      map <unsigned int, double>  element_map;
	      it_bc = boundary_cond_elem.find(elem);
	      if (it_bc !=  boundary_cond_elem.end())
		{//----yes, it belongs to the map
		  element_map = it_bc->second ;
		  new_map.insert(pair<const Elem*, map<unsigned int, double> > (elem,element_map));
		}
	      else
		{//---no, it does not belong to the map
		  
		  //let us check its parent 
		  const Elem* elem_parent = elem->parent();
		  map<const Elem*, std::map <unsigned int , double > > :: iterator it_bc;
		  
		  it_bc = boundary_cond_elem.find(elem_parent);
		  if (it_bc !=  boundary_cond_elem.end())
		    { //top parent belongs to a map
		      
		      map <unsigned int, double>  parent_element_map;
		      parent_element_map = it_bc->second ;
		      
		      
		      map <unsigned int, double>  element_map;
		      
		      element_map.clear();
		      
		      
		      //now we check if a child lies on a necessary side
		      map <unsigned int, double> :: iterator par_it  =  parent_element_map.begin();
		      map <unsigned int, double> :: iterator par_end =  parent_element_map.end();
		      
		      for ( ; par_it != par_end; ++par_it) 
			{
			  unsigned int side = par_it->first;
			  /* I assume that side number i of an element and its parent are parallel faces */
			  /* I hope this is correct */
			  
			  if ((elem->neighbor(side) == NULL)) // on boundary
			    {
			      const double stress =  par_it->second;
			      element_map.insert(pair<unsigned int, double>  (side,stress));
			      
			    }
			}
	     	      
		      if (! element_map.empty() ) new_map.insert(pair<const Elem*, map<unsigned int, double> > (elem,element_map));
		      
		    }  
		}
	    }
	}
      
      boundary_cond_elem = new_map;
    }
}




//-----------------------------------------------------------------//

void Macrostrain::assemble_strain_matrix(EquationSystems& es,
				     const std::string& system_name)

{ //

 
 

  // It is a good idea to make sure we are assembling
  // the proper system.
  assert (system_name == "Strain");

  int temp_i;
  temp_i = 0;

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Matrix Assembly",false);

 
  
  // Get a constant reference to the mesh object.
  const Mesh& mesh = es.get_mesh();

  // The dimension that we are running
   
  //dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem> ("Strain");

  
  unsigned int uvar[3] ;
  unsigned int var_fict;

  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
  

  if (number_of_add_var_static !=0 ) var_fict = system.variable_number("fict");

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  // const DofMap& dof_map = system.get_dof_map();
  DofMap& dof_map = system.get_dof_map();
  
 

  FEType fe_type = dof_map.variable_type(uvar[0]);
 
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
 


  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);
  
  // Tell the finite element object to use our quadrature rule.

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

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element iterators are a nice way to iterate through
  // all the elements, or all the elements that have some property.
  // There are many types of element iterators, but here we will
  // use the most basic type, the  const_elem_iterator.  The iterator
  //  el will iterate from the first to the last element.  The
  // iterator  end_el tells us when to stop.  It is smart to make
  // this one  const so that we don't accidentally mess it up!
//   const_elem_iterator           el (mesh.elements_begin());
//   const const_elem_iterator end_el (mesh.elements_end());

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  //  std::vector<int> mapping(dof_map.n_dofs());
  
  const   unsigned int Number_of_act_el = mesh.n_active_elem(); //number of fict DOFs

  const   unsigned int Number_of_DOFs = dof_map.n_dofs(); //number of total DOFs
  

  
  //dof_map.print_dof_constraints();   	
 
  // Loop over the elements.  Note that  ++el is preferred to
  // el++ since the latter requires an unnecessary temporary
  // object.

  

  Tensor1 vec1;
  Tensor1 vec2;

  Tensor2Sym eps_var;

  Tensor2Sym eps_const;

  Tensor2Sym C_kl;

  // double a_substrate[3];
 
  double lattice_factor;

  //(*crystal_temp)[substr_mat].get_lat_const(substrate_lat_const);
  // for (int i = 0; i <=2; i++) std::cerr << substrate_lat_const[i];

  
  unsigned int el_number = 0; 
  
  system.matrix->zero();






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
      const unsigned int n_dofs   = dof_indices.size(); //in fact, could be  dof_indices.size() - 1, fict is not used 

      fe->reinit  (elem);


      Ke_total.resize (n_dofs + number_of_add_var_static, n_dofs + number_of_add_var_static);
      Fe_total.resize (n_dofs + number_of_add_var_static);

      dof_indices_total.resize(n_dofs + number_of_add_var_static);
      
     
      

      const int material = (*material_of_elem_temp)[el_number];


      

      eps_const = (*crystal_temp)[material].get_const_eps0(substrate_lat_const, eps0_var_log) 
	+ (*eps0_of_elem_temp)[el_number] ;//+ substrate_shear;

      
      for (unsigned int qp=0; qp<qrule.n_points(); qp++) // Scalar product integration
	{//qp
	  //----------------------------------------------------------//
	  //master equation:                                          //
          //                                                          //
	  //     d/dx_i  ( C_ijkl (du_k/dx_l + eps0_kl)) =0           //
	  //                                                          //
          //                                                          //
	  //----------------------------------------------------------//
	  for (unsigned int j = 0; j<=2; j++)
	    {//loop over j: master equation discretization
 	     
	      dof_map.dof_indices (elem, dof_indices_component, uvar[j]);
	      const unsigned int n_u_dofs = dof_indices_component.size(); 
	      Fe_sub.reposition (uvar[j]*n_u_dofs, n_u_dofs);
	 
	          
	      

	      for (unsigned int p1=0; p1<n_u_dofs; p1++)
		{//p1



		  if (!belongs_to_substrate(p1, elem))
		    {//---------not a substrate point, so we calculte the RHS for master equation system
     
		      vec1 = 0;
		      for (int i = 1; i<=dim; i++) vec1(i) = dphi[p1][qp](i-1);
		      
		      //-------------eps0 part------------------
		      for (unsigned int k = 0; k <= 2; k++)
			{
			  vec2 = 0;
			  for (int i = 1; i <=3; i++ ) 
			    {	     	
			      if (k+1 > i)
				vec2(i) = eps_const(k+1,i);
			      else
				vec2(i) = eps_const(i,k+1);
			    }
			  
			 
			  Fe_sub(p1) -= JxW[qp]*(vec1 * ( (*C_tensor_temp)[material].get_subtensor(j+1,k+1)* vec2 ));
			} 

		      //------------external normal stress--------
		     
		      

		      map<const Elem*, std::map <unsigned int , double > > :: iterator it_bc;
		      std::map <unsigned int, double>  element_map;
		     
		      it_bc = boundary_cond_elem_temp->find(elem);
		    
		      if (it_bc != boundary_cond_elem_temp->end()) 
			{ //if this element has applied stress
			 
			 
			 
			  element_map = it_bc->second ;

			  if (dim > 1)
			    {//2D, 3D case 

			      map <unsigned int, double> :: iterator side_it  =  element_map.begin();
			      map <unsigned int, double> :: iterator side_end =  element_map.end();
			      for ( ; side_it != side_end; ++side_it) 
				{//loop over stressed sides
				  unsigned int side = side_it->first;
			     	
				  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
				  
				  const std::vector<Real>& JxW_face = fe_face->get_JxW();
				
				  const std::vector<Point >& qface_point = fe_face->get_xyz();
				  
				  const std::vector<Point> & normal = fe_face->get_normals();
				  
				  double stress_value = side_it->second ;


				  fe_face->reinit(elem, side);

				 
				  for (unsigned int qp=0; qp<qface.n_points(); qp++)
				    {				  
				      Fe_sub(p1) += ((JxW_face[qp] * phi_face[p1][qp]) * stress_value) * normal[qp](j);
				    } 
				  
				    
				}
			    }
			  else
			    { //1D case 
			      if (j == 0) //the only non-zero component of the normal
				{
				  map <unsigned int, double> :: iterator node_it  =  element_map.begin();
				  map <unsigned int, double> :: iterator node_end =  element_map.end();
				  for ( ; node_it != node_end; ++node_it) 
				    {//loop over stressed nodes
				      unsigned int node_number = node_it -> first;
				      if (p1 == node_number)
					{
					  double stress_value = node_it->second ;
					  double normal;
					  Point p = elem->point(node_number);
					  Point pc = elem->centroid();
					  if (p(0) > pc(0)) 
					    normal = 1.0;
					  else
					    normal = -1.0;
                                          cerr << "stress_value   " <<  stress_value << "\n";

					  Fe_sub(p1) += stress_value * normal;
					} 
				    }
				}
			    }
			  
			}
			
		      //---RHS of master equation is done---------------------------------------------------


		      //----- additional variables first, if any---------------------------------------------
		    
		      for (unsigned int i1 = 0; i1 < (*add_var_temp).size()  ; i1++)
			{	  			 
			
			  
			  Ke_u_add_sub.reposition (uvar[j]*n_u_dofs, n_dofs, n_u_dofs, (*add_var_temp).size());
			  
			
			  eps_var = (*crystal_temp)[material].get_var_eps0( (*add_var_temp)[i1].name );
			  
			  
			  
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
			   
			      Ke_u_add_sub(p1,i1) += JxW[qp]*(vec1 * ( (*C_tensor_temp)[material].get_subtensor(j+1,k+1) * vec2 ));
			    }
 
			} 
		      //------ additional variables done ---------------------
		    }
		  else
		    {
		      //substrate point
		      Fe_sub(p1) = 0.0;
		      
		    }




		}		
	    
	      
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
			 
			

			  vec1 = 0;
			  for (int i = 1; i<=dim; i++) vec1(i) = dphi[p1][qp](i-1) ;
			  
			  vec2 = 0;
			  for (int i = 1; i<=dim; i++) vec2(i) = dphi[p2][qp](i-1) ;
			  
			  scal_prod = vec1 * ( (*C_tensor_temp)[material].get_subtensor(j+1,k+1) * vec2);
			  
			  if (!belongs_to_substrate(p1, elem))
			    {
			      Ke_sub(p1,p2) += JxW[qp]*scal_prod;
			     
			    }
			  else
			    {
			      Ke_sub(p1,p2) = delta(p1,p2)*delta(j,k);
			    }

			}
		      
		    }
			
		
		}
		      
		      
	    }//end of loop over j - end of master equation
	  

	  //----------------------------------------------------------------------------
	  //superlattice equations
	  for (unsigned int eq_number =  0;   eq_number <  number_of_add_var_static; eq_number++)
	    {
	      
	      dof_map.dof_indices (elem, dof_indices_component, uvar[0]);
	      const unsigned int n_u_dofs = dof_indices_component.size();
 
	    
	      if ( (*add_var_temp)[eq_number].lat_cons )
		{
		    unsigned int lat_index = (*add_var_temp)[eq_number].index1;
		    double lat_const = (*crystal_temp)[material].lat_const_calc[lat_index - 1];
		    
		    lattice_factor = 1/lat_const;
		 }
	       else
		 {
		   lattice_factor = 1.0;
		 }
		  
		  unsigned int lat_index1 = (*add_var_temp)[eq_number].index1;
		  unsigned int lat_index2 = (*add_var_temp)[eq_number].index2;
		  
		 
		  C_kl = (*C_tensor_temp)[material].get_another_subtensor(lat_index1,lat_index2);
		  //----------------RHS------------------
		  Fe_add_sub.reposition(n_dofs + eq_number,1);
		  

		  Fe_add_sub(0) -=  JxW[qp] * doubleContraction(C_kl , eps_const ) *lattice_factor  ;
		 
		  //-------------------------------------

		  //-------------ux,uy,uz----------------
		  for (unsigned int k = 0; k<=2; k++)
		    {//loop over k
		      dof_map.dof_indices (elem, dof_indices_component, uvar[k]);
		      const unsigned int n_u_dofs = dof_indices_component.size(); 
		      Ke_add_u_sub.reposition (n_dofs, uvar[k]*n_u_dofs, (*add_var_temp).size() , n_u_dofs);

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
		 
		  for (unsigned int i1 = 0; i1 < (*add_var_temp).size()  ; i1++)
		    {
		      Ke_add_add_sub.reposition(n_dofs + eq_number,n_dofs + i1,1,1);
		      eps_var = (*crystal_temp)[material].get_var_eps0( (*add_var_temp)[i1].name );
		    
		      Ke_add_add_sub(0,0) += JxW[qp]  * doubleContraction(eps_var,C_kl) *  lattice_factor;
		   
		      
		    }
		    
		  

		  //------------------------------------------

		
	    
	    }


	  

	  if  (number_of_add_var_static != 0)
	    {
	     
	      if ( el_number >=  number_of_add_var_static )
		{
		  dof_map.dof_indices (elem, dof_indices_component, uvar[0]);
		  const unsigned int n_u_dofs = dof_indices_component.size(); 
		  Ke_sub.reposition(3*n_u_dofs,3*n_u_dofs, 1, 1);
		  Ke_sub(0,0) += 1.0;
		}
		
	      
	      	      
	    }

	  
	  // end of superlattice equations
	  //--------------------------------------------------------------------------

	  
	  
	}

    
     
       
       for (unsigned i =0 ; i < n_dofs; i++)
	 dof_indices_total[i] = dof_indices[i];
	 
       
       for (unsigned i =0 ; i <number_of_add_var_static ; i++)
	 dof_indices_total[i+n_dofs] =  add_dofs_vector[i];
	
       
       

      
      dof_map.constrain_element_matrix_and_vector(Ke_total, Fe_total, dof_indices_total);
    

      system.matrix->add_matrix (Ke_total, dof_indices_total);
      system.rhs->add_vector    (Fe_total, dof_indices_total);
      
      // constraint is necessary for Ke_add!!
      

   

      el_number++;
    }
       
      

 
  
  // system.matrix->print();

 
  // system.rhs->print();

  //-----------------------------------------------------------------------
  //Application of periodicity constraints



  //-----------------------------------------------------------------------
  //dof_map.print_dof_constraints(); 	  	

  //  system.matrix->print();

   
 


  std:: cout<< "matrix is done \n";  
   
  // system.rhs->print();


  std:: cout << "Active dofs number   " << system.n_active_dofs()   	<< "\n";

  std:: cout << "Total dofs number   " << system.n_dofs()   	<< "\n";

  std:: cout << "Constraint dofs number " <<  system.n_constrained_dofs()   	<< "\n";
      


}




//-----------------------------------------------------------------//
 bool Macrostrain::belongs_to_substrate(unsigned int n, const Elem* elem )
 {



   if (grown_on_substrate)
     {
     
       set<unsigned int>::iterator substr_nodes_it;



       substr_nodes_it = substrate_nodes_temp->find( elem->node(n) );

       if (substr_nodes_it != substrate_nodes_temp->end())
	 return(true);
       else
	 return(false);
       


        
     }
  

   else
     {
       if (elem->node(n) == fixed_node1_temp ) return(true);
       //there is always one node that is fixed. 
       // std :: cerr << elem->node(n) << "  " << fixed_node_number << "\n";
       //return(elem->node(n) == fixed_node_number);
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
//------------------------------------------------------
  //referencing of data objects
  C_tensor_temp = &C_tensor;

  crystal_temp =  &crystal ;

  material_of_elem_temp  = &material_of_elem;
 
  eps0_of_elem_temp = &eps0_of_elem;

  number_of_add_var_static = number_of_add_var;

  add_var_temp = &add_var; 
  
  fixed_node1_temp = fixed_node1;
  fixed_node2_temp = fixed_node2;
  fixed_node3_temp = fixed_node3;

  zero_set_dofs_temp = &zero_set_dofs;
  

  boundary_cond_elem_temp = &boundary_cond_elem;

  substrate_nodes_temp = &substrate_nodes;
  
  //------------------------------------------------------
}
//-----------------------------------------------------------------//
Mesh* Macrostrain::get_mesh()
{
  return( &(equation_systems->get_mesh()) );
}

//-----------------------------------------------------------------//
void Macrostrain::solve()

{
  
 
  
  NumericVector<Number>& old_solution = 
    equation_systems->get_system("Strain").add_vector("old solution");
  
  // Initialize the data structures for the equation system.
  equation_systems->init();	

 
  

 
  Mesh& mesh = equation_systems->get_mesh();

 

  assemble_material_list();
  initialize_eps0_list();
  initialize_el_number_map();

  define_fixed_nodes();

  make_nodes_periodic();

  //init_u_node(); //not necessary
  
  init_substrate();

  create_substate_nodes_set();
  
  
  create_bondary_conditions_map();


  

  refer_objects();

 

  set_up_additional_dofs();

  equation_systems->get_system("Strain").solution->zero();

  apply_periodic_bc();

  apply_antirotation_constraints();

  
   

  equation_systems->get_system("Strain").solve();

 
  
  old_solution = *equation_systems->get_system("Strain").solution;
 
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
      
      // The \p ErrorVector is a particular \p StatisticsVector
      // for computing error information on a finite element mesh.
      ErrorVector error;
      
      // The \p ErrorEstimator class interrogates a finite element
      // solution and assigns to each element a positive error value.
      // This value is used for deciding which elements to refine
      // and which to coarsen.
      KellyErrorEstimator error_estimator;

      // Compute the error for each active element using the provided
      // \p flux_jump indicator.  Note in general you will need to
      // provide an error estimator specifically designed for your
      // application.
      error_estimator.estimate_error (equation_systems->get_system("Strain"),
				      error);
		
      // This takes the error in \p error and decides which elements
      // will be coarsened or refined.  Any element within 20% of the
      // maximum error on any element will be refined, and any
      // element within 10% of the minimum error on any element might
      // be coarsened. Note that the elements flagged for refinement
      // will be refined, but those flagged for coarsening _might_ be
      // coarsened.
      
      mesh_refinement.flag_elements_by_error_fraction (error,
							 refine_fraction,
							 coarsen_fraction,
							 max_ref_level);
	

      std :: cerr << " refine_fraction  " << refine_fraction <<"\n";
	
      // This call actually refines and coarsens the flagged
      // elements.
      if (uniform_refinement == 1)
	mesh_refinement.uniformly_refine(1);
      else
	mesh_refinement.refine_and_coarsen_elements();





      // This call reinitializes the \p EquationSystems object for
      // the newly refined mesh.  One of the steps in the
      // reinitialization is projecting the \p solution,
      // \p old_solution, etc... vectors from the old mesh to
      // the current one.
      equation_systems->reinit();

      old_solution = *equation_systems->get_system("Strain").solution;
      

      assemble_material_list();
      initialize_eps0_list();
      initialize_el_number_map();
      
      set_up_additional_dofs();
      
     
      init_substrate();


      define_fixed_nodes();

      update_substrate_nodes_set();
  
      update_bondary_conditions_map();

      refer_objects();


      make_nodes_periodic();
      
      apply_periodic_bc();

      apply_antirotation_constraints();
      
      mesh.print_info();
      
      equation_systems->get_system("Strain").solution->zero();
      
      equation_systems->get_system("Strain").solve();
      
      old_solution.add(-1.0, 
		       *equation_systems->get_system("Strain").solution);

      old_solution.close();

      std::cout << "Norm of the difference  " <<  old_solution.linfty_norm() << "\n";
     
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
    }

  std::cout << "\n" ;
  std::cout << "Final Mesh after  " <<  max_r_steps <<" refinements  steps   " <<  "\n" ;
  mesh.print_info();


  std::cerr << "Grid refinement is done \n";

  if (intermediate_output)
    {
      output_materials("materials_nondeformed.dat");
    }

  cerr << atom_structure_filename << "\n";

  if (calculate_atom_displacements)
    { 
      read_atom_structure(atom_structure_filename);
      
      std::ostringstream disp_file;
      disp_file << atom_displacements_filename << ".out"  ;
     

      write_atom_displacements(disp_file.str());
    }
  //------------------------------------------------------------------------------------
  //geometry relaxation

  init_u_node();

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
      equation_systems->get_system("Strain").solution->zero();
      
      //apply_periodic_bc();
      
      equation_systems->print_info();
      
      equation_systems->get_system("Strain").solve();
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
	      write_atom_displacements(disp_file.str());
	    }

	}
            //---------------------------------------------------

     
      

    }//end of shape loop
    
  std::cout << "SHAPE IS DONE\n";

  std::cout << "OUTPUT...\n";

  update_u_node();
   //------write-----------------------------------------------
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

  if (calculate_atom_displacements)
    {
      std::ostringstream disp_file;
      disp_file << atom_displacements_filename <<".out";
      write_atom_displacements(disp_file.str());
    }
  

  //------------------------------------------------------------------------------------------//

  
}

//-------------------------------------------------------------------------------------//
void Macrostrain::assemble_material_list()
{
  Mesh& mesh = equation_systems->get_mesh();
 

  const unsigned int N_elem = mesh.n_active_elem();

  material_of_elem.resize( N_elem );

 
  
  MeshBase::const_element_iterator el  = mesh.active_elements_begin();
  MeshBase::const_element_iterator end_el = mesh.active_elements_end();

    

  unsigned int el_number = 0;
  
  for ( ; el != end_el ; ++el) 
    {
      const Elem* elem = *el;
	  unsigned int  mat;
	  
	  mat  = (unsigned int) (*meshdata)(elem->top_parent(),0);
	  
	  
	  
	  material_of_elem[el_number] = mat - 1 ;
	  
	  el_number++;

    }
     
}

//-----------------------------------------------------------------//
void Macrostrain::update_eps0_list()
{
  //calculate eps_new = eps_old + 1/2(du/dx + du/dx) 

  const Mesh& mesh = equation_systems->get_mesh();
  //const unsigned int dim = mesh.mesh_dimension();
  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  DofMap& dof_map = system.get_dof_map();
  
  unsigned int uvar[3] ;
  
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
 

 
  

  FEType fe_type = dof_map.variable_type(0);
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
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

  /*



  MeshBase::const_element_iterator el  = mesh.active_elements_begin();
  MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  
  double a_substrate[3];  crystal[substr_mat].get_lat_const(a_substrate);

  const Tensor2Sym eps0(0);
 
  unsigned int el_number = 0;
 

  for ( ; el != end_el ; ++el) 
    {
      eps0_of_elem[el_number] = eps0;

      el_number++;
    }

  */


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

	  if ( ( std::abs( (*node_fix)(dir) - max_coord[dir]) < pos_tol)  ||
	       ( std::abs( (*node_fix)(dir) - min_coord[dir]) < pos_tol)  )
	    {
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
  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

  unsigned int uvar[3] ;
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }

  unsigned int system_number=system.number();
  
  DofMap& dof_map = system.get_dof_map();
  
  FEType fe_type = dof_map.variable_type(uvar[0]);
  
 
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
   

  

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
			  if (elem->on_boundary())
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
			      if (child->on_boundary())
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
  // std::cout << '+++++++++++++++++\n';
  //dof_map.print_dof_constraints();  
  // std::cout << '+++++++++++++++++\n';

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
 LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

 unsigned int uvar[3] ;
 for (unsigned int i = 0; i<= 3 - 1; i++) 
   {
      uvar[i] = system.variable_number(uname_vec[i]);
   }
 
 unsigned int system_number=system.number();
 
 DofMap& dof_map = system.get_dof_map();


  
 FEType fe_type = dof_map.variable_type(uvar[0]);
  
 
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

 DofConstraintRow constraint; 

 if ( !grown_on_substrate )
   {//substrate system is not treated
     if (dim > 1)
       {//1D system does not need a treatment 
	 if (fixed_node1_temp == fixed_node2_temp) 
	   {
	     
	     cerr << "Error: fixed_node1_temp == fixed_node2_temp  " << fixed_node1_temp <<"  "
		  << fixed_node2_temp <<"\n";  
	     cerr << fixed_node1 <<"  "<< fixed_node2 <<"\n";   
	     exit(1);
	   }

	 const Node& node1= mesh.node(fixed_node1_temp);
	 const Node& node2= mesh.node(fixed_node2_temp);

        



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
		     if ((fixed_node1_temp == fixed_node3_temp)|| (fixed_node2_temp == fixed_node3_temp)) 
		       {
		     
			 cerr << "Error: fixed_node1_temp == fixed_node3_temp  or \n";
			 cerr << "Error: fixed_node2_temp == fixed_node3_temp   \n";
			 cerr << fixed_node1 <<"  "<< fixed_node2 << "  "<< fixed_node2  << "\n";   
			 exit(1);
		       }
		 
		     const Node& node3= mesh.node(fixed_node3_temp);
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
//write out strain tensor component------------------------     ---------------
void Macrostrain::output_strain(std::string filename )
{

 
 
  char num_i[2];
  char num_j[2];
  string eps_ij;
  
  double a_substrate[3];  crystal[substr_mat].get_lat_const(a_substrate);

  unsigned int index = 0;

  const Mesh& mesh = equation_systems->get_mesh();

  //const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  DofMap& dof_map = system.get_dof_map();

  unsigned int uvar[3] ;
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
 
  

  FEType fe_type = dof_map.variable_type(0);
 
  
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

 

  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();


  unsigned int Number_of_elements = mesh.n_active_elem();
  std::vector<Number> eps_data(Number_of_elements*6);
  std::vector<std::string> eps_names(6);
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

	    const int material = material_of_elem[elem_number];
	    
	    eps0 = eps0_of_elem[elem_number] + calculate_eps_lat_matching( material); //previous iterations and lattice matching
	
	    Point center1 = FEInterface::inverse_map(dim, fe_type, elem, center);
	    point_vec[0] = center1;
	    fe->reinit (elem, &point_vec);
	    
	   
	    double du_i_over_dx_j = 0;
	 	    
	  
	    dof_map.dof_indices (elem, dof_indices_component1, uvar[i-1]);
	    dof_map.dof_indices (elem, dof_indices_component2, uvar[j-1]);

	    const unsigned int n_u_dofs = dof_indices_component1.size(); 
	  
	    for (unsigned int p1=0; p1<n_u_dofs; p1++)
	      {   

		if (j<=dim) du_i_over_dx_j += 0.5 *  dphi[p1][0](j-1) * (*solution)(dof_indices_component1[p1]);
	
		
		if (i<=dim) du_i_over_dx_j += 0.5 *  dphi[p1][0](i-1) * (*solution)(dof_indices_component2[p1]); 
   

	      }
	    
	    
	  
	    double eps_value = eps0(i,j) + du_i_over_dx_j ;  
	    
	  
	    eps_data[index + elem_number * 6  ] = eps_value; //that's a correct order of variables
	     
	    elem_number++;
	  }

	index++;
      }

  //std :: cout << filename << "\n";
  if (output_type == "GMV")     GMVIO_cell(mesh).write_ascii_cell_data(filename, eps_data, eps_names);

  if (output_type == "tecplot") TecplotIO_cell(mesh,false).write_cell_data(filename,eps_data,eps_names);

}
//---------------------------------------------------------------------------
//writes piezopolarization in GMV format
void Macrostrain::output_piezo(std :: string filename)
{
  char num_i[2];
  const Mesh& mesh = equation_systems->get_mesh();

  unsigned int Number_of_elements = mesh.n_active_elem();

  std::vector<Number> polariz_data(Number_of_elements*3);

  std::vector<std::string> polariz_names(3);
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
  

  if (output_type == "GMV")  GMVIO_cell(mesh).write_ascii_cell_data(filename, polariz_data, polariz_names);
  if (output_type == "tecplot") TecplotIO_cell(mesh,false).write_cell_data(filename, polariz_data, polariz_names);
}

//---------------------------------------------------------------------------

void Macrostrain::move_nodes()
{
  const Mesh& mesh = equation_systems->get_mesh();

  //const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

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

  const Node& node_fix = mesh.node(fixed_node1_temp);  

  for (unsigned int i = 0; i < dim; i++) r0(i + 1) = node_fix(i);

  unsigned int node_number = 0;
  for ( ; ( (nd != nd_end) ) ; ++nd)
    {
      //-----------------------------------------------------------------------------
      //
      Node* node1 = *nd;

      r = Tensor1(0);
      
      for (unsigned int i = 0; i < dim; i++)   r(i + 1) = (*node1)(i); //current node position
	  
      if (number_of_add_var !=0)
	{
	  for (unsigned int i = 0; i < dim; i++)   r(i + 1) -= u_node[node_number][i]; //position without displacement

	  r = lat_matching_transformation * (r - r0) + r0; //transform cells according to the lattice matching transformation

	  for (unsigned int i = 0; i < dim; i++) r(i + 1) += u_node[node_number][i]; //add back displacements
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


      node_number++;

    } 




  //--------------------------------------------------------------------

  std :: cout << "Nodes are moved. \n";

}
//-------------------------------------------------------------------------------------//


Tensor2Sym Macrostrain::get_strain(const Elem* elem, bool crystal_system )
{

  Tensor2Sym eps0;
  Tensor2Sym eps;
 

  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  DofMap& dof_map = system.get_dof_map();


  unsigned int uvar[3] ;
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
  


  FEType fe_type = dof_map.variable_type(0);
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
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
	    
	    
	  }
	    
	    
	  
	eps(i,j) = eps0(i,j) + du_i_over_dx_j ; 
	
      }

  //-----------------------------------------------------------------
  //we have to add lattice matching deformation
  eps += calculate_eps_lat_matching(material_of_elem[elem_number]);
  //------------------------------------------------------------------
  if (crystal_system)
    {//convert to crystal system
      const int material = material_of_elem[elem_number]; //get material number
     
      Tensor2Gen RotM = (crystal[material].RotMatrix).transpose();//get rotation matrix
      
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

//-------------------------------------------------------------------------------------------/
Tensor1 Macrostrain::get_piezopolarization(const Elem* el)
{
  //---------------calculate strain in crystal system---------------------------
  Tensor2Sym strain_cr= get_strain( el, true);

  //---------------get material number -----------------------------------------

  map<const Elem*, unsigned int> :: iterator el_numb_it;

  el_numb_it = elem_numbers.find(el);
  
  const unsigned int elem_number = el_numb_it->second;

  const int material = material_of_elem[elem_number]; //get material number
  //----------------calculate polarization---------------------------------

  Tensor1 polariz = piezo[material].get_polariz_cryst(strain_cr); //crystal system

  

  polariz =(crystal[material].RotMatrix) * polariz; //calculation system

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

  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

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

Tensor2Sym Macrostrain::calculate_eps_lat_matching(unsigned int material)
{
  //constant part of the lattice matching tensor
  Tensor2Sym eps0 = (*crystal_temp)[material].get_const_eps0(substrate_lat_const, eps0_var_log);
  //------
  //variable part:
  for (unsigned int i = 0; i < number_of_add_var; i++)
    {
      unsigned int dof_number = add_dofs_vector[i];

      double coeff = (*equation_systems->get_system("Strain").solution)( dof_number ); 

      eps0 += coeff * crystal[material].get_var_eps0( add_var[i].name );
    } 

  return(eps0);

}
//-------------------------------------------------------------------------------------------/
void Macrostrain::init_substrate()
{
  substrate_shear = Tensor2Sym(0);
  crystal[substr_mat].get_lat_const(substrate_lat_const);
}
//--------------------------------------------------------------------------------------------/

void Macrostrain::update_substrate()
{

  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

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

  const Mesh& mesh = equation_systems->get_mesh();

  u_node.resize(mesh.n_nodes());

  for (unsigned int i = 0; i< mesh.n_nodes(); i++ )  u_node[i] = single_node;

  active_node_number.clear();

  Node* nd1;
  MeshBase::const_node_iterator nd  = mesh.active_nodes_begin();
  MeshBase::const_node_iterator end_nd = mesh.active_nodes_end();

  unsigned int node_number = 0;

  for ( ; nd != end_nd ; ++nd) 
    {
      nd1 = *nd;
      
      active_node_number.insert( pair<Node*,unsigned int> (nd1,node_number)  );
     
      node_number++;
    }

}
//-------------------------------------------------------------------------------------------/
void Macrostrain::update_u_node()
{
  const Mesh& mesh = equation_systems->get_mesh();

  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

  AutoPtr<NumericVector<Number> >& solution = system.solution;
 
  const unsigned int system_number = system.number();

  unsigned int uvar[3] ;
  
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }


  MeshBase::const_node_iterator       nd     = mesh.active_nodes_begin();
  const MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();

  unsigned int node_number = 0;

  for ( ; ( (nd != nd_end) ) ; ++nd)
    {
      Node* node1 = *nd;
      for (unsigned int i = 0; i < 3; i++) //<3 , not < dim (necessary for atoms!) 
	{
	  const unsigned int  n_dof = node1->dof_number(system_number,uvar[i],0);
	  
	  u_node[node_number][i] +=  (*solution)(n_dof);
	}


      node_number++;

    }
}
//-------------------------------------------------------------------------------------------/
void Macrostrain::output_add_strain_variables(string filename)
{
  std::ofstream out;

  if (number_of_add_var !=0) 
    {
      LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

      AutoPtr<NumericVector<Number> >& solution = system.solution;
      

      std::ofstream out (filename.c_str());

      assert (out.good());

      std::cerr << " number_of_add_var  = " << number_of_add_var << "\n"; 


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

  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

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


  LinearImplicitSystem& system = equation_systems->get_system<LinearImplicitSystem> ("Strain");

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  DofMap& dof_map = system.get_dof_map();

  std::vector<unsigned int> dof_indices_component1;

  unsigned int uvar[3] ;
  for (unsigned int i = 0; i<= 3 - 1; i++) 
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }
 
  

  FEType fe_type = dof_map.variable_type(0);
 
  
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

 
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  const std::vector<Point>& q_point = fe->get_xyz();

  vector<Point> point_vec(1);

  //----------------------------------------------------------------------
  double substrate_lat_const_initial[3];
  crystal[substr_mat].get_lat_const(substrate_lat_const_initial);
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
		      map<Node*, unsigned int> :: iterator active_node_it;
		      unsigned int node_number = (active_node_number.find(atom_structure[i].element->get_node(nd))) -> second;
		      u_vector[coord] += u_node[node_number][coord] * phi[nd][0]; 
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
void Macrostrain::output_materials(std :: string filename)
{

  const Mesh& mesh = equation_systems->get_mesh();

 
  std::vector<std::string> mat_name(1);
  mat_name[0] = "material";

  unsigned int N =  material_of_elem.size();

  std::vector<double> mat_double(N);

  for (unsigned i = 0; i < N; i++) mat_double[i] = material_of_elem[i];

  if (output_type == "GMV")  GMVIO_cell(mesh).write_ascii_cell_data(filename, mat_double, mat_name);
  if (output_type == "tecplot") TecplotIO_cell(mesh,false).write_cell_data(filename,mat_double, mat_name);

}

//-------------------------------------------------------------------------------------------/
Macrostrain::~Macrostrain()
{
  delete equation_systems;
}
