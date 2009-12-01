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

#include "Messages.h"

using namespace std;

MacroHeatBalance* MacroHeatBalance::static_this;
Device* MacroHeatBalance::_device;
//-----------------------------------------------------------------//


void MacroHeatBalance::parse_options( )
{

  const ModelOptions& opts = SimulationInterface::get_options();

  myopts.integration_order = static_cast<libMeshEnums::Order>(
      opts.get_option("integration_order", 5));

  myopts.heat_scheme =  opts.get_option("heat_scheme","surface");

  //std::cout<<"Heat Scheme:  "<< myopts.heat_scheme<<std::endl;

  myopts.work_units = opts.get_option("Work_length_units", 1e-2);


  string qrule =_device->get_mesh().mesh_dimension() < 3 ? "trapez" : "gauss";

  //string qrule =mesh->mesh_dimension() < 3 ? "trapez" : "gauss";

  qrule = opts.get_option("quadrature_rule", qrule);
  if (qrule == "gauss")
    myopts.quadrature_type = QGAUSS;
  else if (qrule == "trapez")
    myopts.quadrature_type = QTRAP;
  else
    throw InitFailedException("Unknown quadrature rule");


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


  //------init is done---------------------------------------------------------------------//

}
//-------------------------------------------------------------------------------//
void  MacroHeatBalance::do_solve()
{


  parse_options();

  static_this = this;

  my_system->set_options(get_solver_options());
  my_system->solve();
   double power_emitted = calculate_power_emitted();
   double power_dissipated_rstf = calculate_power_dissipated_rstf();
  //double power_dissipated = calculate_power_dissipated();

   double check = 100 - std::abs((power_emitted - power_dissipated_rstf)/(power_emitted));
   if (power_emitted<1e-10)
     check = 0;


     //std::cout<<"Power Emitted"   <<"          "<<power_emitted<<std::endl;
    //std::cout<<"Power Dissipated"<<"       "   <<power_dissipated<<std::endl;
     //std::cout<<"Power Dissipated rstf"<<"  "   <<power_dissipated_rstf<<std::endl;


     if (SimulationOptions::verbose() > 1)
     {
       ostringstream os;
       os << "Energy Conservation:  " << check << " %";
       Messages::info(os.str());
     }

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

    //heat_model->set_element(elem);
    //heat_model->set_side(-1);

    heat_model->re_init();

    std::vector<double> heat_source;
    heat_model->get_total_heat_source(elem,q_point,heat_source);

    //std::vector<RealGradient> flux_power;
    //    heat_model->get_total_power_flux(q_point,flux_power);

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

	  Fe(p1) +=JxW[qp] * heat_source[qp] * phi[p1][qp];

      }//end Loop over quadrature points
    } // end loop over test function

    //Boundary conditions and source

    //The loop over element is the only loop that is surviving at this point

    const unsigned int num_sides = elem->n_sides();

    for (unsigned int side = 0; side<num_sides; side++)
    {

      const ElementSide elside(elem->top_parent(), side);

      //  heat_model->set_side(side);

      Boundary* bd = se.get_boundary(elside);

      if (bd != NULL)
      {
	if (bd->get_boundary_properties( get_id() ) != NULL )
	{

	  // heat_model->re_init();
	  fe_face->reinit(elem,side);

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
  // heat_model->set_element(elem);
  //heat_model->set_side(-1);
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


  //if (variables.find("HeatSource") != varend)
  //{
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

    // std::vector<std::set<ID> > source_index(nm);

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
	//	source_index[i].insert(leg->first);
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

    //heat_model->set_element(elem);
    //heat_model->set_side(-1);

    heat_model->re_init();
    heat_model->get_thermal_conductivity(kappa);

    if (HS != -1)
    {

      unsigned int k = 0;

      for (int i = 0; i < nm; i++)
      {
        std::vector<std::map<ID, double> > heat_sources;

        heat_model->get_heat_source_model(ids[i])->get_heat_sources(elem,_node,heat_sources);

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

	heat_model->get_total_heat_source(elem,_node,total_heat_source);
	results[id + HS + k] =  total_heat_source[0];

      }

    } //if (HS != -1)


    if (PF != -1)
    {
      unsigned int k = 0;

      std::vector< std::map< ID, double > > jq_solution;
      get_solution_secure(elem,_node,JQ_var,jq_solution);

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
	heat_model->get_heat_source_model(ids[i])->get_power_fluxes(elem,_node,power_flux);


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

     //  if (variables.count("TotalFlux")  ||
// 	  variables.count("thermal")    ||
// 	  variables.count("PowerFlux") )
//       {

// 	std::vector<RealGradient > total_power_flux;

// 	heat_model->get_total_power_flux(_node,total_power_flux);

// 	double Px_tot = total_power_flux[0](0) + Pqx;
// 	double Py_tot = total_power_flux[0](1) + Pqy;
// 	double Pz_tot = total_power_flux[0](2) + Pqz;

//  	switch (dim)
//  	{
//  	case 3:
//  	  results[id + W[k] + 2] = Pz_tot;
// 	case 2:
// 	  results[id + W[k] + 1] = Py_tot;
// 	  results[id + W[k] + dim] = sqrt(Px_tot * Px_tot + Py_tot * Py_tot + Pz_tot * Pz_tot);
// 	default:
// 	  results[id + W[k] ] = Px_tot;
// 	}

//       }

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

 double power = calculate_power_dissipated();

     values.resize(1,power);


  }
}



double
MacroHeatBalance::calculate_power_dissipated(void)
{

  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return 0;

  TiberLinearSystem& system = *my_system;

  const unsigned int  var = system.variable_number("T");

  const unsigned int dim = mesh->mesh_dimension();

  DofMap& dof_map =  system.get_dof_map();

  FEType fe_type = dof_map.variable_type(var);

  AutoPtr<FEBase> fe_face (build_finite_element(dim, fe_type));


  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order =  myopts.integration_order;

  AutoPtr<QBase> qface(QBase::build(
				    myopts.quadrature_type, dim-1, integration_order));


  fe_face->attach_quadrature_rule(qface.get());


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

  double  power_dissipated = 0.0;

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

	//heat_model->set_element(elem);
	//heat_model->set_side(s);

	heat_model->re_init();

	fe_face->reinit(elem, s);

	std::vector< std::map< ID, double > > jq_solution;

        get_solution_secure(elem,q_point,JQ_var,jq_solution);

	RealGradient P;

	for (unsigned int qp = 0; qp <  qface->n_points(); qp++)
	{

	  P(0) = jq_solution[qp].find(JQX)->second;
	  P(1) = jq_solution[qp].find(JQY)->second;
	  P(2) = jq_solution[qp].find(JQZ)->second;
	  if (dim> 1)
	    power_dissipated += JxW[qp] * P * face_normals[qp];
	  else
	    power_dissipated +=  P * face_normals[qp];

	}
      }
    } // end loop over elem sides
  } // end loop over elements

  return power_dissipated;
}



double
MacroHeatBalance::calculate_power_dissipated_rstf(void)
{
  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return 0;

  TiberLinearSystem& system = *my_system;

  const unsigned int  var = system.variable_number("T");

  const unsigned int dim = mesh->mesh_dimension();

  DofMap& dof_map =  system.get_dof_map();

  FEType fe_type = dof_map.variable_type(var);

  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type));


  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order =  myopts.integration_order;

  AutoPtr<QBase> q_rule(QBase::build(
				    myopts.quadrature_type, dim, integration_order));
  SimulationEnvironment&  env= get_environment();

  fe->attach_quadrature_rule(q_rule.get());

  std::vector<unsigned int> dof_indices;

  // Jacobian * quadrature weight at each integration point.
  const vector<Real>& JxW = fe->get_JxW();

  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  // element shape function gradients
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();


  // will contain the node ids if an element has boundary nodes
  vector<Boundary*> node_ids;

  MeshBase::const_element_iterator el =
                                  mesh->active_elements_begin();
  const MeshBase::const_element_iterator end_el =
                                  mesh->active_elements_end();



  double  power_dissipated = 0.0;
  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    const Elem* top_parent = (*el)->top_parent();

    fe->reinit(elem); //centroid
    ID subdomain = elem->subdomain_id();

    const Material* mat = _device->get_material(subdomain);

    HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );
    //heat_model->set_element(elem);
    //heat_model->set_side(-1);
    heat_model->re_init();
    std::vector<double> heat_source;
    heat_model->get_total_heat_source(elem,q_point,heat_source);

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices,var);

    bool has_node = false;
    node_ids.resize(elem->n_nodes());
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    {
      Boundary* bd = env.get_boundary(elem->get_node(n));
      node_ids[n] = bd;
      if (bd != NULL)
        has_node = true;
    }


    // if the element has no node on a boundary,
    // we can go to the next element
    if (!has_node)
      continue;


    // get DOF indices
    dof_map.dof_indices(elem, dof_indices, var);



     std::vector< std::map< ID, double > > jq_solution;
     get_solution_secure(elem,q_point,JQ_var,jq_solution);

    for (unsigned int qp = 0; qp < q_rule->n_points(); qp++)
    {

      unsigned int n_dofs = dof_indices.size();
      // get the solution values at the centroid


      RealGradient P;

      P(0) = jq_solution[qp].find(JQX)->second;
      P(1) = jq_solution[qp].find(JQY)->second;
      P(2) = jq_solution[qp].find(JQZ)->second;

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {

        Boundary* boundary = node_ids[n];
        if (boundary != NULL)
        {
	  power_dissipated  += JxW[qp] * (P * dphi[n][qp] + heat_source[qp] * phi[n][qp]);
          //power_dissipated  += P * dphi[n][qp];
        }

      }
    } // end loop over quadrature points
  } // end loop over elements

  return power_dissipated;
}



double
MacroHeatBalance::calculate_power_emitted(void)
{

  // we only do something if we are on processor 0
  if (libMesh::processor_id() != 0)
    return 0;

  TiberLinearSystem& system = *my_system;

  const unsigned int  var = system.variable_number("T");

  const unsigned int dim = mesh->mesh_dimension();

  DofMap& dof_map =  system.get_dof_map();

  FEType fe_type = dof_map.variable_type(var);

  AutoPtr<FEBase> fe (build_finite_element(dim, fe_type));


  libMeshEnums::Order integration_order;
  if (dim == 1)
    integration_order = libMeshEnums::CONSTANT;
  else
    integration_order =  myopts.integration_order;

  AutoPtr<QBase> qface(QBase::build(
				    myopts.quadrature_type, dim, integration_order));


  fe->attach_quadrature_rule(qface.get());



  const vector<Real>& JxW = fe->get_JxW();

  // physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz();

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();

  std::vector<unsigned int> dof_indices;


  MeshBase::const_element_iterator el =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =     mesh->active_local_elements_end();

  const SimulationEnvironment& env = get_environment();

  double PowerEmitted = 0.0;

  for ( ; el != end_el ; ++el)
  {

    const Elem* elem = *el;

    fe->reinit(elem);

    ID subdomain = elem->subdomain_id();

    // get DOF indices
    dof_map.dof_indices(elem, dof_indices,var);

    const Material* mat = _device->get_material(subdomain);

    HeatModel* heat_model =  (  dynamic_cast<HeatModel*> (  mat -> get_model(get_id()) )  );

    //heat_model->set_element(elem);

    //heat_model->set_side(-1);

    heat_model->re_init();

    std::vector<double> heat_source;
    heat_model->get_total_heat_source(elem,q_point,heat_source);

    for (unsigned int qp = 0; qp <  qface->n_points(); qp++)
      PowerEmitted  += JxW[qp] * heat_source[qp];



  } // end loop over elements

  return  PowerEmitted;
}


void
MacroHeatBalance::do_print_info(void)
{

  string space("  ");
  //cout << space << "linear solver is: petsc" <<std::endl;

}
