// $Id$

#include "Macrostrain.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "BoundaryProperties.h"
#include "Boundary.h"
#include "MacrostrainBoundaryProperties.h"
#include "MacrostrainPressure.h"
#include "MacrostrainSubstrate.h"
#include "SimulationOptions.h"
#include "TiberLinearSystem.h"
#include "MeshUtils.h"

#include "TiberLinearSolver.h"
#include "TiberPetscLinearSolver.h"

//Needed for straining atomistic structure
#include "AtomisticStructure.h"

#include "Messages.h"
#include "Specie.h"


using namespace std;
//-----------------------------------------------------------------//



Macrostrain* Macrostrain::static_this;




void Macrostrain::do_setup_solution_variables(void)
{
  declare_solution(Strain, TENSOR, CELL, "");
  declare_solution(Stress, TENSOR, CELL, "GPa");
  declare_solution(Displacement, VECTOR, NODES, "m");
  declare_solution(vonMises, REAL, CELL, "GPa");
  declare_solution(PiezoPolarization, VECTOR, CELL, "C/m^2");
  declare_solution(StrainCrystal, TENSOR, CELL, "");
  declare_solution(StressCrystal, TENSOR, CELL, "GPa");
  declare_solution(EnergyDensity, REAL, CELL, "J/m^3");
}



void
Macrostrain::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& points)
{

  ID subdomain = elem->subdomain_id();
  MacrostrainModel* macrostrain_model =
      dynamic_cast<MacrostrainModel*>(get_physical_model(subdomain));

  const RotatedCrystal& crystal = macrostrain_model->get_material()->get_rotated_crystal();
  const Tensor2Gen& RotM = crystal.RotMatrix;

  // this is in crystal system
  const Tensor2Sym& str_cryst = result_strain[elem];
  Tensor2Sym elemstr = sym(RotM * (str_cryst * RotM.transpose()));


  if (values.count(Displacement))
  {
    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem& system = *my_system;
    const NumericVector<Number>& solution = *(system.solution);

    double x0 = get_scaling().get_length_scaling();

        unsigned int uvar[3] ;

    for (unsigned int i = 0; i < 3; i++)
      uvar[i] = system.variable_number(uname_vec[i]);

    DofMap& dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(uvar[0]);
    AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));
    const vector<vector<double> >& phi = fe->get_phi();

    vector<unsigned int> dof_indices_x;
    vector<unsigned int> dof_indices_y;
    vector<unsigned int> dof_indices_z;
    dof_map.dof_indices(elem, dof_indices_x, uvar[0]);
    dof_map.dof_indices(elem, dof_indices_y, uvar[1]);
    dof_map.dof_indices(elem, dof_indices_z, uvar[2]);

    fe->reinit(elem, &points);

    const unsigned int n_points = points.size();

    for (unsigned int i = 0; i < n_points; i++)
    {
      double ux = 0;
      double uy = 0;
      double uz = 0;

      for (unsigned int j = 0; j < dof_indices_x.size(); j++)
      {
        ux += solution(dof_indices_x[j]) * phi[j][i];
        uy += solution(dof_indices_y[j]) * phi[j][i];
        uz += solution(dof_indices_z[j]) * phi[j][i];
      }

      values[Displacement][3 * i]   = ux * x0;
      values[Displacement][3 * i + 1] = uy * x0;
      values[Displacement][3 * i + 2] = uz * x0;
    }
  }



  if (values.count(Strain))
  {
    values[Strain][0] = elemstr(1,1);
    values[Strain][1] = elemstr(2,2);
    values[Strain][2] = elemstr(3,3);
    values[Strain][3] = elemstr(2,1);
    values[Strain][4] = elemstr(3,2);
    values[Strain][5] = elemstr(3,1);
  }

  bool do_stress = values.count(Stress);
  bool do_vonmises = values.count(vonMises);
  bool do_stress_cryst = values.count(StressCrystal);
  bool do_energy = values.count(EnergyDensity);
  Tensor2Sym elemstress;

  if (do_stress || do_vonmises || do_stress_cryst || do_energy)
  {
    //Elasticity in the calculation system
    Stiffness* C_tensor_el = macrostrain_model->get_stiffness();
    const Tensor4DSym& C_calc =  C_tensor_el->C_calc;

    // add converse piezo stress
    macrostrain_model->get_converse_piezo_stress(elemstress, elem);
    elemstress = elemstr * C_calc - elemstress;

    if (do_stress)
    {
      values[Stress][0] = elemstress(1,1);
      values[Stress][1] = elemstress(2,2);
      values[Stress][2] = elemstress(3,3);
      values[Stress][3] = elemstress(2,1);
      values[Stress][4] = elemstress(3,2);
      values[Stress][5] = elemstress(3,1);
    }

    if (do_vonmises)
    {
      double e11 = elemstress(1,1);
      double e22 = elemstress(2,2);
      double e33 = elemstress(3,3);
      double e21 = elemstress(2,1);
      double e31 = elemstress(3,1);
      double e32 = elemstress(3,2);
      double I1 = e11 + e22 + e33;
      double I2 = e11*e22 + e22*e33 + e11*e33 - e21*e21 - e31*e31 - e32*e32;
      double J2 = I1 * I1 / 3.0 - I2;
      values[vonMises][0] = sqrt(3 * J2);
    }

    if (do_energy)
    {
      double energy = doubleContraction(elemstress, elemstr);
      // stress is in GPa
      values[EnergyDensity][0] = 0.5 * energy * 1e9;
    }

    if (do_stress_cryst)
    {
      elemstress = sym(RotM.transpose() * (elemstress * RotM));
      values[StressCrystal][0] = elemstress(1,1);
      values[StressCrystal][1] = elemstress(2,2);
      values[StressCrystal][2] = elemstress(3,3);
      values[StressCrystal][3] = elemstress(2,1);
      values[StressCrystal][4] = elemstress(3,2);
      values[StressCrystal][5] = elemstress(3,1);
    }
  }


  if (values.count(PiezoPolarization))
  {
    Tensor1 piezo = get_piezopolarization(elem);

    values[PiezoPolarization][0] = piezo(1);
    values[PiezoPolarization][1] = piezo(2);
    values[PiezoPolarization][2] = piezo(3);
  }


  if (values.count(StrainCrystal))
  {
    values[StrainCrystal][0] = str_cryst(1,1);
    values[StrainCrystal][1] = str_cryst(2,2);
    values[StrainCrystal][2] = str_cryst(3,3);
    values[StrainCrystal][3] = str_cryst(2,1);
    values[StrainCrystal][4] = str_cryst(3,2);
    values[StrainCrystal][5] = str_cryst(3,1);
  }
}


//-------------------------------------------------------------------------//



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






//-----------------------------------------------------------------//
void Macrostrain::parse_options( )
{

 const ModelOptions& opt = get_options();




 max_r_steps = opt.get_option("refinement_steps", 0);

 uniform_refinement = opt.get_option("uniform_refinement", false);


 refine_fraction = opt.get_option("refine_fraction", 0.25);

 coarsen_fraction = opt.get_option("coarsen_fraction", 0.0);

 max_ref_level = opt.get_option("max_refinement_level",10);
 max_shape_steps = opt.get_option("number_shape_steps",0);


 calculate_atom_displacements = opt.find_option("strain_atomistic_structure");
 structure_to_be_strained = opt.get_option("strain_atomistic_structure", "none");


 AtomisticStructure* as =
   get_environment().get_device().get_atomistic_structure(structure_to_be_strained);
 if (as == NULL) calculate_atom_displacements = false;

 if (calculate_atom_displacements && max_shape_steps == 0)
 {
   Messages::warning("Atomic strain requires number_shape_steps > 0");
   Messages::warning("number_shape_steps will be set to 1");
   max_shape_steps = 1;
 }

 internal_strain = opt.get_option("internal_strain_correction",true);


 fix_all_fixed_points = opt.get_option("fix_all_fixed_points", false);

 _preallocate = opt.get_option("preallocate_matrix", _preallocate);

 // assert(periodicity_x == false);



 // set default options for solver
 ModelOptions& solver_opts = get_solver_options();


 if (solver_opts.get_option("method", "") == "")
   solver_opts["method"] = (dim == 1) ? "gmres" : "bcgsl";

 if (solver_opts.get_option("preconditioner", "") == "")
   solver_opts["preconditioner"] = (dim == 1) ? "ilu" : "jacobi";





 if (!grown_on_substrate)
 {

   apply_antirotation = opt.get_option("apply_antirotation", true);


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
  //if (!get_options().find_option("constant_strain"))
  {
  StrainSimulation::do_init();

  const ModelOptions& options = get_options();


  equation_systems = & (get_equation_systems());


  SimulationEnvironment& si = get_environment();




  get_scaling().set_length_scaling(  _device->get_mesh_units() );


  uname_vec[0]="ux";
  uname_vec[1]="uy";
  uname_vec[2]="uz";

  mesh = &(equation_systems->get_mesh());

  dim = mesh->mesh_dimension();

  substrate_name = "";
  grown_on_substrate = false;
  substrate_crystal = NULL;

  SimulationEnvironment::BoundaryIterator bit(si.boundaries_begin());
  SimulationEnvironment::BoundaryIterator bend(si.boundaries_end());
  for ( ; bit != bend; ++bit)
  {
    Boundary* bd = *bit;
    MacrostrainBoundaryProperties* bp =
        static_cast<MacrostrainBoundaryProperties*>(bd->get_boundary_properties(get_id()));

    if (bp->get_type() == "substrate")
    {
      MacrostrainSubstrate* sub = static_cast<MacrostrainSubstrate*>(bp);
      substrate_name = bd->get_name();
      grown_on_substrate = true;
      Material* mat = sub->get_material();

      if (mat == NULL)
        throw InitFailedException("Macrostrain: Substrate \'" + substrate_name +
            "\' has no material.");

      substrate_crystal =  &(mat->get_rotated_crystal());
      break;
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
      if (   MeshUtils::may_belong_to_element(elem,  p) )
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

    if (elem1 == NULL)
      throw InitFailedException("Macrostrain: reference point is  wrong (probably, it does not belong to the simulation domain)");

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
  //my_system = TiberLinearSystem::create(*equation_systems,
  //    get_equation_system_name(), get_solver_options());
  create_equation_system("linear");
  my_system = &get_equation_system<TiberLinearSystem>();

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


  }
  //------init is done---------------------------------------------------------------------//

  parse_options();
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




  int verbose = SimulationOptions::verbose();


  if (verbose > 2)
    Messages::info("Starting the matrix assembly");


  int temp_i;
  temp_i = 0;

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Matrix Assembly",false);



  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

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

    // cerr << el_number << "\n";

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
      //!first we calculate volume partg
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

	    //-------------eps0 part and converse piezo volumic part

	    for (int i = 1; i<=dim; i++)
	    {
	      double stress;

	      if (i > j+1)
		stress = stress_converse_piezo(i, j+1);
	      else
		stress = stress_converse_piezo(j+1, i);

	      Fe_sub(p1) -= JxW[qp] * dphi[p1][qp](i-1) *
		( doubleContraction(C_tensor_el->get_another_subtensor(i,j+1), eps_const) - stress );
	    }

	  }



	  //! may be there is external pressure or extended device, so we have to calculate a surface part
	  //
	  //
	  //
	  for (unsigned int side=0; side<num_sides; side++)
	  {//side loop
	    Boundary* bd = si.get_boundary(ElementSide(elem,side));


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



		  const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();


		  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();

		  const std::vector<Real>& JxW_face = fe_face->get_JxW();

		  const std::vector<Point >& qface_point = fe_face->get_xyz();

		  const std::vector<Point> & normal = fe_face->get_normals();

		  fe_face->reinit(elem, side);

		  for (unsigned int qp=0; qp< qface.n_points(); qp++)
		  {



		    Tensor2Gen t = doubleContraction(C_tensor_el->C_calc, eps_const) - stress_converse_piezo;

		    Tensor1 v;
		    v(1) = normal[qp](0);  v(2) = normal[qp](1);  v(3) = normal[qp](2);

		    Tensor1 v1 = t*v;


		    Fe_sub(p1) += JxW_face[qp] * phi_face[p1][qp] * v1(j+1);




		  }

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
	      //! Here we impose the substrate point (or fixed point) is not coupled at all.
	      if ( belongs_to_substrate(p2, elem) && (p1 != p2) )
	      	  Ke_sub(p1,p2) = 0;

	    }

	  }

	}


	//if there are extended boundary condition we need surface integration!
	for (unsigned int side=0; side<num_sides; side++)
	{//side loop
	  Boundary* bd = si.get_boundary(ElementSide(elem,side));
	  if (bd != NULL && (bd->get_boundary_properties( get_id() ) != NULL )  )
	    if (dynamic_cast<MacrostrainBoundaryProperties*>( bd->get_boundary_properties( get_id() ))->get_type() == "extended")
	    {



	      const std::vector<std::vector<RealGradient> >& dphi_face = fe_face->get_dphi();

	      const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();

	      const std::vector<Real>& JxW_face = fe_face->get_JxW();

	      const std::vector<Point >& qface_point = fe_face->get_xyz();

	      const std::vector<Point> & normal = fe_face->get_normals();

	      fe_face->reinit(elem, side);




	      for (unsigned int p1=0; p1<n_u_dofs; p1++)
	      {
		if (!belongs_to_substrate(p1, elem))
		  for (unsigned int p2=0; p2<n_u_dofs; p2++)
		  {


		    for (unsigned int qp=0; qp<qface.n_points(); qp++)
		    {

		      vec2 = 0;
		      for (int i = 0; i < dim; i++) vec2(i+1) = dphi_face[p2][qp](i) ;


		      adjust_derivatives(vec2,  normal[qp]);


		      vec1 = 0;
		      for (int i = 0; i < dim; i++) vec1(i+1) = phi_face[p1][qp]* normal[qp](i);


		      const double scal_prod =  vec1 * ( C_tensor_el->get_subtensor(j+1,k+1) *vec2);

		      Ke_sub(p1,p2) -= JxW_face[qp] * scal_prod ;



		    }
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


	    if ( !belongs_to_substrate(p1, elem))
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
  static int counter = 0;


  std::ostringstream os;
  os << "rhs" << counter << ".m";

  //system.matrix->close();

  //system.matrix->print_matlab(os.str().c_str());

  system.rhs->close();

  system.rhs->print_matlab(os.str().c_str());

  counter++;
}

*/
   //-----------------------------------------------------------------------
   //Application of periodicity constraints



   //-----------------------------------------------------------------------
   //dof_map.print_dof_constraints();

   //   system.matrix->print();

   // system.rhs->print();

   if (verbose > 2) Messages::info("Matrix assembly done");

   if (verbose > 4)
   {
     ostringstream os;
     os << "Active DOFs      " << system.n_active_dofs() << Messages::endl;
     os << "Total DOFs       " << system.n_dofs() << Messages::endl;
     os << "Constrained DOFs " <<  system.n_constrained_dofs();
     Messages::info(os.str());
   }


/*
   {

     for (int i = 0; i < system.matrix->n(); i++)
     {
       cerr << (*system.matrix)(i,i) << "\n";
     }

   }
*/
}




//-----------------------------------------------------------------//
 bool Macrostrain::belongs_to_substrate(unsigned int n, const Elem* elem )
 {

   // return(false);
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


     if (elem->node(n) == fixed_node1 && apply_antirotation)
       return(true);
     else if ((dim > 1) && (fix_all_fixed_points && (elem->node(n) == fixed_node2)))
       return true;
     else if ((dim > 2) && (fix_all_fixed_points && (elem->node(n) == fixed_node3)))
       return true;
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
MeshBase* Macrostrain::get_mesh()
{
  return( &(equation_systems->get_mesh()) );
}

//-----------------------------------------------------------------//

//-----------------------------------------------------------------//
void Macrostrain::do_solve()

{

  int verbose = SimulationOptions::verbose();


  parse_options();


  SimulationEnvironment& si = get_environment();



  //------------------------------------------------------
  NumericVector<Number>& old_solution =
    my_system->get_vector("old solution");




  //------------------------------------------------------


  MeshBase& mesh = equation_systems->get_mesh();




  initialize_eps0_list();

  initialize_el_number_map();

  define_fixed_nodes();

  make_nodes_periodic();

  init_u_node();

  init_substrate();



  refer_objects();

  if (grown_on_substrate) si.get_boundary_nodes (substrate_name, substrate_points);

  set_up_additional_dofs();

  if (_first_run)
  {
    
    my_system->solution->zero();
    
    apply_periodic_bc();
    
    
    if (verbose > 2) cout << "apply_antirotation_constraints ... " << flush;
    
      if (!fix_all_fixed_points)
        if(apply_antirotation) apply_antirotation_constraints();
      
      if (verbose > 2) cout << "done \n" << flush;
      
  }




  if (!_is_reallocated && _preallocate)
  {
    reallocate_matrix();
    _is_reallocated = true;
  }

  my_system->set_options(get_solver_options());

  //if (verbose > 2) cout << "Assemble and solve the linear system ...\n" << flush ;

  my_system->solve();

  //if (verbose > 2) cout << "The linear system is solved... \n" << flush ;



  old_solution = * (my_system->solution);

  update_substrate();


  if (intermediate_output)
  {

    if (output_type=="GMV") GMVIO (mesh).write_equation_systems ("displacement_field.dat.000", *equation_systems);
    if (output_type=="tecplot") TecplotIO_cell(mesh,false).
      write_equation_systems ("displacement_field.dat.000", *equation_systems);
    
    output_strain("strain.dat.000");
    output_add_strain_variables("add_var.000");
  }


  //if (dim > 1) mesh.write("mesh0.ucd");
  //-----------------------------------------------------




  //refinement loop----------------------------------------------------------------
  MeshRefinement mesh_refinement(mesh);

  for (unsigned int r_step = 1; r_step <= max_r_steps; ++r_step)
  {
    if (verbose > 0)  cout << "\nRefining the mesh... (Step" << r_step << ")\n" << flush;
    

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
    
    if (apply_antirotation) apply_antirotation_constraints();
    
    //mesh.print_info();
    
    my_system->solution->zero();
    
    
    
    my_system->solve();
    
    
    if (verbose > 1)
      std::cout << "Norm of the difference  " << norm_of_difference( old_solution, *(my_system->solution) )
                << "  after step number " << r_step << "\n" << flush;
    
    if (intermediate_output)
    {
      
      std::ostringstream os;
      os << "displacement_field.dat.00" << r_step;
      
      
      if (output_type=="GMV")  
        GMVIO (mesh).write_equation_systems (os.str(), *equation_systems);
      if (output_type=="tecplot")   
        TecplotIO_cell(mesh,false).write_equation_systems(os.str(), *equation_systems);
      
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
    
    if (verbose > 1)
    {
      
      //std::cout << "\n" ;
      //std::cout << "Final Mesh after  " <<  max_r_steps <<" refinements  steps   " <<  "\n" ;
      //mesh.print_info();
      //std::cout << flush;
    }
    
    if (verbose > 0) Messages::info("Grid refinement finished");
    
  }
  
  // Setup for atom displacemets: compute relative coords of unstrained atoms
  //  with respect to the original mesh
  if (calculate_atom_displacements)
  {
    
    //Initialise relative points
    AtomisticStructure* as = NULL;
    as = get_environment().get_device().get_atomistic_structure(structure_to_be_strained);
    
    LinearImplicitSystem& system = *my_system;

    DofMap& dof_map = system.get_dof_map();

    FEType fe_type = dof_map.variable_type(0);


    AutoPtr<FEBase> fe (build_finite_element(dim, fe_type, true));

    _atom_relative_points.resize(as->get_structure_atoms().size());
    for (unsigned int i = 0; i < as->get_structure_atoms().size(); i++)
    {
      if (as->get_structure_atoms()[i].get_elem() != NULL)
      {
        Point tmp_point;
        tmp_point(0) = as->get_structure_atoms()[i].get_position()(1) / as->get_scale();
        tmp_point(1) = as->get_structure_atoms()[i].get_position()(2) / as->get_scale();
        tmp_point(2) = as->get_structure_atoms()[i].get_position()(3) / as->get_scale();

        //get atom relative point
        _atom_relative_points[i] =  
          FEInterface::inverse_map(dim, fe_type, as->get_structure_atoms()[i].get_elem(), tmp_point);
      }
    }
  }


  //------------------------------------------------------------------------------------
  // geometry relaxation (mesh deformation)

  for (unsigned int geom_it = 1 ; geom_it <= max_shape_steps; geom_it++)
  {
    if (verbose > 1)  cout << "\n Geometry relaxation of the mesh... (Step" << geom_it << ")\n" << flush;
    
    
    if (verbose > 2) Messages::info("Update nodes... ", false);
    
    update_u_node();
    
    if (verbose > 2) Messages::info("done");
    
    
    //equation_systems->print_info();
    //------move nodes------------------------------------------
    
    update_eps0_list();
    
    
    refer_objects();
    
    move_nodes();
    
    update_substrate();
    
    
    //---------------------------------------------------------
    
    //-------solve---------------------------------------------
    my_system->solution->zero();
    
    //apply_periodic_bc();
    
    //equation_systems->print_info();
    
    
    if (verbose > 2) cout << "Will solve...   " << flush;
    
    my_system->solve();
    
    if (verbose  > 2) cout << "solved \n" << flush;
    
    
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
            
    }
    //---------------------------------------------------
    
    
  }//end of shape loop


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
    
    if (verbose > 0) Messages::info("Applying strain to atoms");

    apply_atom_displacements(structure_to_be_strained);
    
    if (internal_strain) 
    {
      if (verbose > 0) Messages::info("Applying internal strain");
      internal_strain_correction(structure_to_be_strained);
    }

    if (verbose > 0)
      Messages::info("Saving strained structure in 'strained.xyz'");

    get_environment().get_device().get_atomistic_structure(structure_to_be_strained)
      ->print_structure("strained.xyz");

  }

  //--------------------------------------------------------------------------------------------------//
  calculate_result_elem_strain_map();


  _first_run = false;


}


//-----------------------------------------------------------------//
void Macrostrain::update_eps0_list()
{
  //calculate eps_new = eps_old + 1/2(du/dx + du/dx)

  const MeshBase& mesh = equation_systems->get_mesh();

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

  const MeshBase& mesh = equation_systems->get_mesh();

  const unsigned int N_elem = mesh.n_active_elem();

  eps0_of_elem.resize( N_elem, Tensor2Sym(0) );




}
//------------------------------------------------------------------------
void Macrostrain::initialize_el_number_map()
{
  const MeshBase& mesh = equation_systems->get_mesh();
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
  const MeshBase& mesh = equation_systems->get_mesh();
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
  const MeshBase& mesh = equation_systems->get_mesh();

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


 const MeshBase& mesh = equation_systems->get_mesh();
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

  const MeshBase& mesh = equation_systems->get_mesh();

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
void Macrostrain::prepare_stress_data_for_output( std::vector<std::string>& stress_names, std::vector<double>& stress_data )
{



       //write names
       unsigned int index = 0;
       char num_i[2];
       char num_j[2];
       string S_ij;
       stress_names.resize(6);
       const MeshBase& mesh = equation_systems->get_mesh();
       unsigned int Number_of_elements = mesh.n_active_elem();
       stress_data.resize(Number_of_elements*6);
       for (int i = 1; i <=3 ; i++)
            {
              for (int j = 1; j <=i; j++)
              {
	         sprintf( num_i, "%i",i);
	         sprintf( num_j, "%i",j);

                 S_ij = "stress_" + string(num_i) + string(num_j);

	         stress_names[index] = S_ij ;
                index++;
             }

        }

        //Write stress


	MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

       	unsigned int elem_number = 0;
	for ( ; el != end_el ; ++el)
	  {
	    const Elem* elem = *el;
            ID subdomain = elem->subdomain_id();
            const Material* mat = _device->get_material(subdomain);
            MacrostrainModel* macrostrain_model;
            macrostrain_model = dynamic_cast<MacrostrainModel*>(   mat ->get_model(get_id())     );



            //strain in calculation system
	    Tensor2Sym strain_el = result_strain[elem];
            const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());
            Tensor2Gen RotM = (crystal_el->RotMatrix);//get rotation matrix
            strain_el = sym(RotM * (strain_el * RotM.transpose())); //calculation system

            //Elasticity in the calculation system
            Stiffness* C_tensor_el;
            C_tensor_el = macrostrain_model->get_stiffness();
            Tensor4DSym C_calc =  C_tensor_el->C_calc;

            Tensor2Sym stress_el = strain_el * C_calc;

            //Write stress

             index = 0;
             for (int i = 1; i <=3 ; i++)
             {
              for (int j = 1; j <=i; j++)
               {
                  double stress_value = stress_el(i,j);
                  stress_data[index + elem_number * 6  ] = stress_value; //that's a correct order of variables
                  index++;
               }
             }

                elem_number++;
             }

}

//--------------------------------------------------------------------------------//
//write out strain tensor component------------------------     ---------------
void Macrostrain::output_strain(std::string filename )
{

  std::vector<std::string> eps_names;

  std::vector<double> eps_data;

  prepare_strain_data_for_output(  eps_names,  eps_data );

  const MeshBase& mesh = equation_systems->get_mesh();

  if (output_type == "GMV")     GMVIO_cell(mesh).write_ascii_cell_data(filename, eps_data, eps_names);

  if (output_type == "tecplot") TecplotIO_cell(mesh,false).write_cell_data(filename,eps_data,eps_names);

}
//---------------------------------------------------------------------------

void Macrostrain::prepare_polarization_data_for_output( std::vector<std::string>& polariz_names, std::vector<double>& polariz_data )
{

  char num_i[2];
  const MeshBase& mesh = equation_systems->get_mesh();

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

  const MeshBase& mesh = equation_systems->get_mesh();



  std::vector<double> polariz_data;

  std::vector<std::string> polariz_names;


  prepare_polarization_data_for_output( polariz_names,  polariz_data );


  if (output_type == "GMV")  GMVIO_cell(mesh).write_ascii_cell_data(filename, polariz_data, polariz_names);

  if (output_type == "tecplot") TecplotIO_cell(mesh,false).write_cell_data(filename, polariz_data, polariz_names);
}

//---------------------------------------------------------------------------

void Macrostrain::move_nodes()
{
  const MeshBase& mesh = equation_systems->get_mesh();

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

  //std :: cout << "Nodes are moved. \n";

}
//-------------------------------------------------------------------------------------//

void Macrostrain::calculate_result_elem_strain_map()
{
  result_strain.clear();
  const MeshBase& mesh = equation_systems->get_mesh();
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

   if (get_options().find_option("constant_strain"))

   {
     std::vector<double> strain(6, 0.0);
     get_options().get_option("constant_strain", strain);
     //Note: upper part is given as Tensor is simmetric

     eps(1,1) = strain[0]; eps(2,1) = strain[2]; eps(2,2) = strain[3];
     eps(3,1) = strain[4]; eps(3,2) = strain[5]; eps(3,3) = strain[6];

   return eps;
   }


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
  Tensor2Sym eps(0);

  if (get_options().find_option("constant_strain"))

  {
    std::vector<double> strain(6, 0.0);
    get_options().get_option("constant_strain", strain);
    eps(1,1) = strain[0]; eps(2,1) = strain[2]; eps(2,2) = strain[3];
    eps(3,1) = strain[4]; eps(3,2) = strain[5]; eps(3,3) = strain[6];
  }

  else

  {


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

  }

   if (crystal_system)
    {//convert to crystal system



      ID subdomain = elem->subdomain_id();

      const Material* mat = _device->get_material(subdomain);

      const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

      Tensor2Gen RotM = (crystal_el->RotMatrix).transpose();//get rotation matrix

      Tensor2Gen eps1 = (RotM*eps)*RotM.transpose();  //transform

      eps  = sym(eps1); //result has to be symmetric

      assert (::norm(eps - eps1) < 1e-6); //is it really symmetric

      for (unsigned int i = 0; i++; i<3)
      {
        for (unsigned int j = 0; j++; j<3)
        {
          std::cerr << "crystal_system_eps " << j + 1 << i + 1 << " " << eps(j + 1, i + 1);
        }
      }
      
      return(eps); //return strain tensor in crystal system
    }
   else
   {
     for (unsigned int i = 0; i++; i<3)
     {
       for (unsigned int j = 0; j++; j<3)
       {
         std::cerr << "calc_system_eps " << j + 1 << i + 1 << " " << eps(j + 1, i + 1);
       }
     }
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


  const MeshBase& mesh =  equation_systems->get_mesh();

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

  const MeshBase& mesh = equation_systems->get_mesh();
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
  const MeshBase& mesh = equation_systems->get_mesh();

  LinearImplicitSystem& system = *my_system;

  AutoPtr<NumericVector<Number> >& solution = system.solution;

  const unsigned int system_number = system.number();

  unsigned int uvar[3] ;

  for (unsigned int i = 0; i<= 3 - 1; i++)
    {
      uvar[i] = system.variable_number(uname_vec[i]);
    }



  vector <double> du(3);

  map<const Node*, vector<double> > temp;
  {
    MeshBase::const_element_iterator el  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

    for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;

      const unsigned int num_nodes = elem->n_nodes();

      for (short i = 0; i < num_nodes; i++)
      {
	const Node* nd = elem->get_node(i);

	if (temp.find(nd) == temp.end())
	{

	  for (unsigned int i = 0; i < 3; i++) //<3 , not < dim (necessary for atoms!)
	  {
	    const unsigned int  n_dof = nd->dof_number(system_number,uvar[i],0);
	    du[i] = (*solution)(n_dof);
	  }
	  temp.insert( pair<const Node*,vector<double> > (nd,du)  );
	}


      }

    }

  }



  map<const Node*, vector<double> >::iterator it = u_node.begin();
  map<const Node*, vector<double> >::iterator it_end = u_node.end();

  for( ; it != it_end ; ++it)
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
  const MeshBase& mesh =  equation_systems->get_mesh();

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

void
Macrostrain::apply_atom_displacements(const std::string structure_name)
{//-------------------------------------------------------------------


  //Get pointer to structure
  AtomisticStructure* as = NULL;
  as = get_environment().get_device().get_atomistic_structure(structure_name);


  const MeshBase& mesh =  equation_systems->get_mesh();

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

  //std::cout << " q_point " << q_point[0] << "becomes ";

  vector<Point> point_vec(1);


  //----------------------------------------------------------------------
  double substrate_lat_const_initial[3];



  substrate_crystal->get_lat_const(substrate_lat_const_initial);
  //----------------------------------------------------------------------


  unsigned int Number_of_atom = as->get_structure_atoms().size();

  //We need to keep track of displacements to move hydrogen atoms in NULL element
  //(apply the same displacement as their neighbour)
  std::vector<Tensor1> u_atm(Number_of_atom);

  for (unsigned int i = 0; i < Number_of_atom ; i++)
    {//atoms loop
      if (as->get_structure_atoms()[i].get_elem() != NULL)
        {//atom belongs to the simulation domain

          vector <double> new_pos_of_atom(3,0.0);

          //Point tmp_point;  (Alex: apparently not used)
          //tmp_point(0) = as->get_structure_atoms()[i].get_position()(1) / as->get_scale();
          //tmp_point(1) = as->get_structure_atoms()[i].get_position()(2) / as->get_scale();
          //tmp_point(2) = as->get_structure_atoms()[i].get_position()(3) / as->get_scale();

          //get atom relative point
          point_vec[0] =  _atom_relative_points[i];

          fe->reinit(as->get_structure_atoms()[i].get_elem(), &point_vec);

          /*----------------------------------------------
            first we move atom because the grid moves
            in 3D this is enough to get into account the external strain
            ---------------------------------------------*/
          for (short coord = 0; coord < dim; coord++)
            {
            new_pos_of_atom[coord] = q_point[0](coord);

            }
          /*----------------------------------------------
           in 2D and 1D there are displacements in the direction perpendicular to the simulation space
           if the reference lattice is fixed, this is enough
           ------------------------------------------------*/
          if (dim < 3)
            {//dim = 1,2
              vector <double> u_vector(3,0.0);
              for (unsigned int nd=0; nd < as->get_structure_atoms()[i].get_elem()->n_nodes(); nd++)
                {
                  for (short coord = dim; coord < 3; coord++)
                    {
                      u_vector[coord] += u_node[as->get_structure_atoms()[i].get_elem()->get_node(nd)][coord] * phi[nd][0];
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
                  for (short i1 =1; i1 <= 3; i1++)                  pos_vector_math(i1) =  new_pos_of_atom[i1 - 1];

                  pos_vector_math = substr_deform*pos_vector_math;

                  for (short i1 =1; i1 <= 3; i1++)                 new_pos_of_atom[i1 - 1] =  pos_vector_math(i1);
                  //---------------------------------------------------------------------------------------

                }

            }

          //-------------------------------------------------------



          Tensor1 tmp(0);
          //std::cout << "setting position from " << as->get_structure_atoms()[i].get_position();
          tmp(1) = new_pos_of_atom[0] * as->get_scale();

          if (dim >= 2) tmp(2) = new_pos_of_atom[1] * as->get_scale();
	  else tmp(2) = as->get_structure_atoms()[i].get_position()(2);

	  if (dim == 3) tmp(3) = new_pos_of_atom[2] * as->get_scale();
	  else tmp(3) = as->get_structure_atoms()[i].get_position()(3);

          u_atm[i] = tmp - as->get_structure_atoms()[i].get_position();

          as->get_structure_atoms()[i].set_position(tmp);

        }//end of non NULL element atoms loop

    }//end of atoms loop


  std::vector< std::vector< unsigned int > > bond_map = as->get_bond_map();
  for (unsigned int i = 0; i < Number_of_atom ; i++)
      {//atoms loop
        if (as->get_structure_atoms()[i].get_elem() == NULL)
          {
          Tensor1 tmp(0);
          tmp = as->get_structure_atoms()[i].get_position();
          //Use first neighbour displacement
          if (bond_map[i].size() == 0) Messages::error("One atom has no neighbours!");
          as->get_structure_atoms()[i].set_position(tmp + u_atm[bond_map[i][0]]);
          }
      }

}
//-------------------------------------------------------------------------------------------/
// This is to shift the cations sublattice in order to make the tethraheron with equal bonds
// Only for nitrides works because of Specie::N
 
void 
Macrostrain::internal_strain_correction(const std::string structure_name)
{

  AtomisticStructure* as = get_environment().get_device().get_atomistic_structure(structure_name);

  std::vector<Atom>& structure = as->get_structure_atoms();

  std::set<ID> IDs = as->get_IDset();

  // create a map between region IDs and lattice constants (for later fast access)
  // (a0 and c0 are the relaxed lattice constants of each material or alloy)
  std::map<ID,double> a, c;
  std::map<ID,const Tensor2Gen*> RotM;

  for(std::set<ID>::iterator reg = IDs.begin(); reg != IDs.end(); reg++)
  {
    const Material* mat = as->get_device()->get_material( (*reg) );

    if (mat->get_structure() == "wz")
    {
      Database db = mat->get_database();
      db.set_section("lattice");
      a[*reg] = db.get("a",0.0);
      c[*reg] = db.get("c",0.0);

      const RotatedCrystal& cry =  mat->get_rotated_crystal();
      RotM[*reg] = &cry.RotMatrix;

      //std::cout << "RotM: "<< cry.RotMatrix(1,3) << " " << cry.RotMatrix(2,3)  << " " << cry.RotMatrix(3,3)  <<std::endl;
    }
    else
    {
      a[*reg] = 0.0;
      RotM[*reg] = NULL;
    }

  }

  // Main Loop on structure

  unsigned int Number_of_atoms = structure.size();

  for (unsigned int i = 0; i < Number_of_atoms ; i++)
  { 

    ID id = structure[i].get_region_ID();

    if (structure[i].get_specie() == Specie::N && RotM[id]!=NULL)
    {
      
      // take strain in crystal system
      Tensor2Sym eps = get_strain(structure[i].get_elem(),true);


      double exx = eps(1,1);
      double eyy = eps(2,2);
      double ezz = eps(3,3);

      double a0 = a[id];
      double c0 = c[id];     
      
      // compute internal displacement:
      // du = a0^2/(3*c0^2) * (exx+eyy-2ezz) * c0 

      Tensor1 du_cry(0), du(0), ro(0), r(0); 

      du_cry(3) = a0*a0/(3*c0) * (exx+eyy-2*ezz) * 10.0; //1 nm -> 10.0 Angstroms   
       
      du = *(RotM[id]) * du_cry;

      ro = structure[i].get_position(); 

      r = ro + du;

      //std::cout<< "ro: "<< ro(1) <<" "<< ro(2) <<" "<< ro(3) << std::endl;
      //std::cout<< "r: "<< r(1) <<" "<< r(2) <<" "<<r(3) << std::endl;
     
      structure[i].set_position(r);

    }

       
  }
    
}

//-------------------------------------------------------------------------------------------/
unsigned int Macrostrain::find_nearest_node(Point& point)
{
  //finds a node number nearest to the point
  const MeshBase& mesh =  equation_systems->get_mesh();
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

void Macrostrain::adjust_derivatives(Tensor1& deriv_vectors, const Point& normal)
{
  short index;

  if ( ( abs(normal(0)) >= abs(normal(1)) ) && ( abs(normal(0)) >= abs(normal(2) )) )
    index  =  0;
  else if ( ( abs(normal(1)) >= abs(normal(0)) ) && ( abs(normal(1)) >= abs(normal(2) ))  )
    index  =  1;
  else
    index  =  2;

  assert(abs(normal(index)) > 1e-10);

  deriv_vectors(index + 1) = 0;

  if ( index == 0 )
    deriv_vectors(1) = -normal(1)/normal(index)*deriv_vectors(2) - normal(2)/normal(index)* deriv_vectors(3);


  else if ( index == 1 )
    deriv_vectors(2) =  -normal(2)/normal(index)*deriv_vectors(3) - normal(0)/normal(index)* deriv_vectors(1);


  else if ( index == 2 )
   deriv_vectors(3) =  -normal(0)/normal(index)*deriv_vectors(1) - normal(1)/normal(index)* deriv_vectors(2);





}

//-------------------------------------------------------------------------------------------//
void Macrostrain::reallocate_matrix(void)
{
  {
    LinearImplicitSystem& system = *my_system;

    DofMap& dof_map = system.get_dof_map();

    std::vector< unsigned int > n_nz =	dof_map.get_n_nz();
    std::vector< unsigned int > n_oz =	dof_map.get_n_oz();


    int n = system.matrix->n();

    for (int i = 0; i < n; i++)
      n_nz[i] += add_dofs_vector.size() ;


    for (short j = 0 ; j <add_dofs_vector.size(); j++)
      n_nz[ add_dofs_vector[j] ] = n;





    //!system matrix
    Mat _matr =  (dynamic_cast<PetscMatrix<Number>* >(system.matrix) ) ->mat();

    int ierr;

    int non_zeros[n];

    for (int i = 0; i < n; i++)
      non_zeros[i] = n_nz[i];

    ierr = MatDestroy(_matr);

    ierr = MatCreateSeqAIJ (PETSC_COMM_WORLD, n, n, 0, non_zeros, &_matr);

    ierr = MatSetFromOptions (_matr);



  }

}
//-------------------------------------------------------------------------------------------//
Macrostrain::Macrostrain(const ModelOptions& options)
  : StrainSimulation(options),
    _is_reallocated(false),
    _preallocate(false),
    _first_run(true),
    poisson_equation(NULL),
    my_system(NULL),
    apply_antirotation(true)

{
}


//Tensor2Sym
//Macrostrain::get_stress_crystal(const Elem* el)
//{


//  Stiffness* C_tensor_el;

//  MacrostrainModel* macrostrain_model;

//   ID subdomain = el->subdomain_id();

//    const Material* mat = _device->get_material(subdomain);

//    const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

//    macrostrain_model = dynamic_cast<MacrostrainModel*>(   mat ->get_model(get_id())     );

//    C_tensor_el = macrostrain_model->get_stiffness();

//    Tensor4DSym C_calc =  C_tensor_el->C_calc;

//    Tensor2Sym strain =  get_strain_crystal(el, el->centroid());

//   Tensor2Sym stress = strain * C_calc;



//  return stress;
//}

// void  Macrostrain::write_atom_displacements(const std::string filename)
// {//-------------------------------------------------------------------



//  //file opening
//   string  string_from_file;

//   std::ifstream atom_in_file;
//   std::ofstream displacement_file;





//   displacement_file.open(filename.c_str());

//   if (!displacement_file.good())
//   {
//     cerr << "Error: file with atom displacements can not be opened\n";
//     cerr << filename.c_str() << "\n";
//   }




//   //--------------------------------------------------------------------
//  if (atom_output_type=="uptight")
//    {
//      //for uptight we have output displacements.
//      //therefore we need also the initial positions

//      atom_in_file.open(atom_structure_filename.c_str());

//      if ( ! atom_in_file.good () )
//        {
// 	 cerr << "Error: file with atom positions can not be opened\n";
// 	 cerr << atom_structure_filename.c_str() << "\n";
//        }

//    }
//   //-------------------------------------------------------------------
//   const MeshBase& mesh =  equation_systems->get_mesh();


//   LinearImplicitSystem& system = *my_system;

//   AutoPtr<NumericVector<Number> >& solution = system.solution;

//   DofMap& dof_map = system.get_dof_map();

//   std::vector<unsigned int> dof_indices_component1;

//   unsigned int uvar[3] ;
//   for (unsigned int i = 0; i<= 3 - 1; i++)
//     {
//       uvar[i] = system.variable_number(uname_vec[i]);
//     }



//   FEType fe_type = dof_map.variable_type(0);


//   AutoPtr<FEBase> fe (build_finite_element(dim, fe_type, true)); //no scaling here


//   const std::vector<std::vector<Real> >& phi = fe->get_phi();

//   const std::vector<Point>& q_point = fe->get_xyz();

//   vector<Point> point_vec(1);

//   //----------------------------------------------------------------------
//   double substrate_lat_const_initial[3];



//   substrate_crystal->get_lat_const(substrate_lat_const_initial);
//   //----------------------------------------------------------------------


//   unsigned int Number_of_atom = atom_structure.size();

//   for (unsigned int i = 0; i < Number_of_atom ; i++)
//     {//atoms loop
//       if (atom_structure[i].element != NULL)
// 	{//atom belongs to the simulation domain

// 	  vector <double> new_pos_of_atom(3,0.0);




// 	  point_vec[0] =  atom_structure[i].relative_point ;




// 	  fe->reinit(atom_structure[i].element, &point_vec);


// 	  /*----------------------------------------------
// 	    first we moove atom because the grid moves
// 	    in 3D this is enough to get into account the external strain
// 	    ---------------------------------------------*/
// 	  for (short coord = 0; coord < dim; coord++)
// 	    new_pos_of_atom[coord] = q_point[0](coord);

// 	 /*----------------------------------------------
// 	   in 2D and 1D there are displacements in the direction perpendcular to the simulation space
// 	   if the reference lattice is fixed, this is enough
// 	   ------------------------------------------------*/
// 	  if (dim < 3)
// 	    {//dim = 1,2
// 	      vector <double> u_vector(3,0.0);
// 	      for (unsigned int nd=0; nd < atom_structure[i].element->n_nodes(); nd++)
// 		{
// 		  for (short coord = dim; coord < 3; coord++)
// 		    {


// 		      u_vector[coord] += u_node[atom_structure[i].element->get_node(nd)][coord] * phi[nd][0];
// 		    }
// 		}


// 	      for (short i1 = 1; i1 < 3; i1++)  new_pos_of_atom[i1] += u_vector[i];

// 	      /*------------
// 		If the reference lattice in a parallel space may change
// 		then we need additional things
// 	      ---------------*/
// 	      if (!grown_on_substrate) /*otherwise the parallel space is fixed*/
// 		{ //free standing
// 		  //----------------------------------------------------------------------------------
// 		  //creation of the transformation tensor
// 		  Tensor2Sym substr_deform(1);

// 		  for (short i1 = dim + 1; i1 <= 3; i1++)
// 		    {
// 		      substr_deform(i1,i1)  += ( (substrate_lat_const[i1 - 1] - substrate_lat_const_initial[i1 - 1])/
// 			 substrate_lat_const_initial[i1 - 1] ) ;
// 		      for (short j1 = dim + 1; j1 < i1; j1++)
// 			{
// 			  substr_deform(i1,j1) +=   substrate_shear(i1,j1);
// 			}

// 		    }
// 		  //--------------------------------------------------------------------------------------
// 		  //application of the transformation
// 		  Tensor1 pos_vector_math;
// 		  for (short i1 =1; i1 <= 3; i1++)		    pos_vector_math(i1) =  new_pos_of_atom[i1 - 1];

// 		  pos_vector_math = substr_deform*pos_vector_math;

// 		  for (short i1 =1; i1 <= 3; i1++)		   new_pos_of_atom[i1 - 1] =  pos_vector_math(i1);
// 		  //---------------------------------------------------------------------------------------

// 		}




// 	    }

// 	  //-------------------------------------------------------




// 	  //-------------------------------------------------------------------
// 	  //output of new coordinates

// 	  if (atom_output_type=="povray")
// 	    {
// 	      displacement_file <<  setw(20) <<  atom_structure[i].mat_number << ",";
// 	      displacement_file <<  setw(20) <<  atom_structure[i].type << ",";


// 	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[0]<< ","  ;
// 	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[1]<< ","  ;
// 	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[2]<< ","  ;

// 	      if (output_strain_on_atoms)
// 		{
// 		  Tensor2Sym epsilon = get_strain(atom_structure[i].element);
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(1,1) << "," ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(2,2) << "," ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,3) << "," ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(2,1) << "," ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,1) << "," ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,2) << "," ;
// 		}
// 	    }


// 	  if (atom_output_type=="uptight")
// 	    {


// 	      getline(atom_in_file, string_from_file );


// 	      istringstream input_string(string_from_file);


// 	      int t ;
// 	      int mat;
// 	      double x;
// 	      double y;
// 	      double z;



// 	      input_string >>mat >> t >> x >> y >> z ;



// 	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[0] - x << ","  ;
// 	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[1] - y << ","  ;
// 	      displacement_file <<  setw(20) <<  setprecision(12) << new_pos_of_atom[2] - z << ","  ;

// 	      if (output_strain_on_atoms)
// 		{
// 		  Tensor2Sym epsilon = get_strain(atom_structure[i].element);
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(1,1)  ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(2,2)  ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,3)  ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(2,1)  ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,1)  ;
// 		  displacement_file <<  setw(20) <<  setprecision(12) << epsilon(3,2)  ;
// 		}




// 	    }

// 	  displacement_file <<  '\n';


// 	}//


//     }//end of atoms loop






// }
// 
// //-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------/
// void Macrostrain::read_atom_structure(const std::string filename)
// {
//   string  string_from_file;

//   std::ifstream atoms_file;

//   atoms_file.open(filename.c_str());
//   //----------------------------------------------------------------


//   //----------------------------------------------------------------

//   if (!atoms_file.good())
//     {
//       cerr << "Error: file with atom coordinates  "<< filename << "   can not be opened";
//     }

//   //-----------------------------------------------------------
//   //determination of number of atoms
//   unsigned int N_atoms = 0;

//   while (getline(atoms_file, string_from_file ))
//     {

//       N_atoms++;
//     }


//   atoms_file.close();
//   //------------------------------------------------------------
//   atom  atom_from_file;

//   atom_from_file.mat_number = 0;
//   atom_from_file.type = 0;
//   atom_from_file.relative_point = Point(0.0, 0.0, 0.0);
//   atom_from_file.element = NULL;

//   atom_structure.resize(N_atoms,atom_from_file);

//   cout << "Number of atoms in " << filename << "  file  " <<  N_atoms;
//   //-------------------------------------------------------------------
//   //-------------------------------------------------------------------
//   //---reading of atoms and processing---------------------------------


//   //--------mesh related objects----------------------------------------
//   const MeshBase& mesh =  equation_systems->get_mesh();

//   LinearImplicitSystem& system = *my_system;

//   DofMap& dof_map = system.get_dof_map();

//   FEType fe_type = dof_map.variable_type(0);





//   //--------------------------------------------------------------------
//   std::ifstream atoms_file1(filename.c_str());

//   for (unsigned int i = 0; i < N_atoms; i++)
//     {//loop over atoms
//       getline(atoms_file1, string_from_file );



//       istringstream input_string(string_from_file);

//       unsigned int n ;
//       int t ;
//       int mat;
//       double x;
//       double y;
//       double z;
//       input_string >>mat >> t >> x >> y >> z ;

//       vector<double> coordinate(3);
//       coordinate[0] = x; coordinate[1] = y; coordinate[2] = z;



//       //   cerr << x << "    "<< y <<"   "<< z <<"\n";



//       atom_structure[i].mat_number = mat;
//       atom_structure[i].type  = t;

//       //-----------------------------------------------------------
//       //determination if the atom belongs to the similation domain
//       //-----------------------------------------------------------
//       Point point2;

//       for (unsigned int i1 = 0; i1 < dim; i1++)
// 	{
// 	  point2(i1) = coordinate[i1];

// 	}

//       //-------------------------------------------------------------
//       //find element that contains the point

//       unsigned int refinement_level = 0;
//       MeshBase::const_element_iterator el3  = mesh.level_elements_begin(refinement_level);
//       MeshBase::const_element_iterator end_el3 = mesh.level_elements_end(refinement_level);

//       Elem*  elem1;

//       bool   found = false;

//       for ( ; ( (el3 != end_el3) ) ; ++el3)
// 	{
// 	  Elem* elem = *el3;

// 	  if  (MeshUtils::may_belong_to_element(elem,  point2))
// 	    {
// 	      if (elem->contains_point(point2))
// 		{
// 		  elem1 = elem;
// 		  found = true;
// 		  break;
// 		}

// 	    }

// 	}

//       if (!found)
// 	{
// 	  //atom is not found and will be ignored
// 	  cerr << "WARNING: atom does not belong to the macroscopic domain\n";
// 	  cerr << "The atom number  " << i << "  will be ignored\n";
// 	  cerr << x <<"   "<< y <<"   " << z <<"   "<<  mat <<"   "<< t <<"\n";

// 	}
//       else
// 	{ //atom is found and will be processed
// 	  //children of the  most coarse element
// 	  while ( !( elem1->active() ) )
// 	    {

// 	      for (unsigned int i=0 ; i < elem1->n_children() ; i++)
// 		{
// 		  Elem* 	child = elem1->child(i);
// 		  if  (MeshUtils::may_belong_to_element(elem1,  point2))
// 		    {
// 		      if (child->contains_point(point2))
// 			{
// 			  elem1 = child;
// 			  break;
// 			}
// 		    }

// 		}
// 	    }


// 	  atom_structure[i].relative_point = FEInterface::inverse_map(dim, fe_type, elem1, point2);

// 	  atom_structure[i].element = elem1;
// 	  //-------------------------------------------------------------------
// 	}


//     }

//   atoms_file1.close();


// }
