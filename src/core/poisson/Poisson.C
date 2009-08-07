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
#include "tensor_value.h"

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



  unsigned int n_vars = 0;
  legend.resize(0);

  int POT = -1;
  if (variables.count("ElPotential") ||
      variables.count("PoissonVariables") )
  {

    POT = n_vars;
    legend.push_back("Potential[V]");
    legend.resize(legend.size() + 1);
    n_vars ++;

  }


  int POL = -1;
  if (variables.count("Polarization") ||
      variables.count("PoissonVariables") )
  {

    POL = n_vars;
    legend.resize(legend.size() + dim + 1);

    switch (dim)
    {
    case 3:
      legend[POL + 2] = "P_z";
      n_vars++;
    case 2:
      legend[POL + 1] = "P_y";
      n_vars++;
      legend[POL + dim] = "modP";
      n_vars++;
    default:
      legend[POL] = "P_x";
      n_vars++;
    }
  }

  int EF  = -1;
  if (variables.count("ElectricField") ||
      variables.count("PoissonVariables") )
  {

    EF = n_vars;
    legend.resize(legend.size() + dim + 1);

    switch (dim)
    {
    case 3:
      legend[EF + 2] = "E_z";
      n_vars++;
    case 2:
      legend[EF + 1] = "E_y";
      n_vars++;
      legend[EF + dim] = "modE";
      n_vars++;
    default:
      legend[EF] = "E_x";
      n_vars++;
    }
  }

  const unsigned int nn  = mesh->n_nodes();
  legend.resize(n_vars);
  results.resize(nn * n_vars,0.0);


  // Get the number of elements that share each node.  We will
  // compute the average value at each node.
  std::vector<unsigned short int> node_conn(nn);
  {
  MeshBase::const_element_iterator it =
    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =
    mesh->active_local_elements_end();

  for ( ; it != end; ++it)
    for (unsigned int n=0; n<(*it)->n_nodes(); n++)
      node_conn[(*it)->node(n)]++;
  }


  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();

  std::vector<unsigned int> dof_indices;

  DofMap& dof_map = my_system->get_dof_map();

  //Build electrif field ID
  set<ID> E_ID;
  E_ID.insert(EX);
  E_ID.insert(EY);
  E_ID.insert(EZ);

  std::vector<std::map<ID, double> > solution;
  //----

  Tensor1 E(0);
  RealGradient _tot_pol(0);
  for ( ; it != end; ++it)
    {
      const Elem* elem = *it;
      dof_map.dof_indices (elem, dof_indices);

      init_poisson_model(elem);
      poisson_model->get_total_polarization(_tot_pol);

      std::vector<Point> p(1);
      p[0] = elem->centroid();
      get_solution_secure(elem,p,E_ID,solution);
      E(1) = solution[0].find(EX)->second;
      E(2) = solution[0].find(EY)->second;
      E(3) = solution[0].find(EZ)->second;

      for (ID n = 0; n < elem->n_nodes(); n++)
      {
	assert (node_conn[elem->node(n)] != 0);
	double conn = static_cast<double>(node_conn[elem->node(n)]);

	unsigned int id =  (elem->node(n) * n_vars);

	if (POT != -1)
	  results[POT + id]  =  (*(my_system->solution))(dof_indices[n]);


	if (POL != -1)
	{
	  switch (dim)
	  {
	  case 3:
	    results[id + POL + 2] += _tot_pol(2) * 1e4/conn;
	  case 2:
	    results[id + POL + 1] += _tot_pol(1) * 1e4/conn ;
	    results[id + POL + dim] = 0.0;
	  default:
	    results[id + POL ] += _tot_pol(0) * 1e4/conn;
	  }
	}


       	if (EF != -1)
	{


	  switch (dim)
	  {
	  case 3:
	    results[id + EF + 2] += E(3)/conn;
	  case 2:
	    results[id + EF + 1] += E(2)/conn ;
	    results[id + EF + dim] = norm(E);
	  default:
	    results[id + EF ] += E(1)/conn;
	  }
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

  RealTensor epsilon(0);
  RealGradient _tot_pol(0);

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

    init_poisson_model(elem);

    poisson_model->get_charge_density(q_point,charge_density);
    poisson_model->get_dielectric_constant(epsilon);
    poisson_model->get_total_polarization(_tot_pol);

    //    std::cout<<_tot_pol<<std::endl;

    for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)  //loop over quadrature points
	for (unsigned int p2=0; p2<n_dofs; p2++) //loop over basis functions
	  Ke(p1,p2) += JxW[qp] * dphi[p1][qp] * (epsilon * dphi[p2][qp]);

    //Fe construction----

    for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)  //loop over quadrature points
      {
	//	Fe(p1) += JxW[qp] * charge_density[qp] * phi[p1][qp];
	total_charge +=  JxW[qp] * charge_density[qp];
	Fe(p1) += JxW[qp] * (dphi[p1][qp] * _tot_pol);

      }


    const unsigned int num_sides = elem->n_sides();

    //Surface///
    for (unsigned int side = 0; side<num_sides; side++)
    {
      const ElementSide elside(elem->top_parent(), side);

      if (se.is_on_boundary(elside))
      {
	fe_face->reinit(elem,side);
	for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
	  for (unsigned int qp = 0; qp < qface.n_points(); qp++){}
	    // {Fe(p1) +=   JxW_face[qp] * (_tot_pol * normal[qp]) * phi_face[p1][qp];}

	for (unsigned int p1=0; p1<n_dofs; p1++) // loop over test function
	  for (unsigned int qp=0; qp<qrule.n_points(); qp++)  //loop over quadrature points
	    for (unsigned int p2=0; p2<n_dofs; p2++){} //loop over basis functions
	      //Ke(p1,p2) -= JxW_face[qp] * (epsilon * dphi_face[p2][qp]) * normal[qp] * phi_face[p1][qp];



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



  std::vector<Point> points(elem->n_nodes());

  for (unsigned n = 0 ; n< elem->n_nodes(); ++n)
  {
    points[n] = elem->point(n);
  }

  get_solution_secure(elem,points,ids,values);


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
  RealTensor epsilon(0);
  RealGradient normals(0);


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
    poisson_model->get_dielectric_constant(epsilon);
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


	  RealGradient D(0),E(0);
	   double flux = 0.0;

 	   E(0) = electric_field[qp].find(EX)->second;
 	   E(1) = electric_field[qp].find(EY)->second;
 	   E(2) = electric_field[qp].find(EZ)->second;


           D = epsilon * E;

	   for( unsigned int kt =0;kt<3;kt++)
	     normals(kt)=normal[qp](kt);


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
