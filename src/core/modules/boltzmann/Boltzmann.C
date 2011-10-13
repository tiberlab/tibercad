// $Id: ThermalBalance.C 2457 2011-03-06 23:52:12Z gromano $

#include "Boltzmann.h"
#include "BoltzmannModel.h"
#include "BoltzmannBoundaryModel.h"
#include "TiberLinearSystem.h"
#include "Messages.h"
#include "ModelOptions.h"
#include "equation_systems.h"
#include "dof_map.h"
#include "quadrature_gauss.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "SimulationOptions.h"
#include "fe_interface.h"

//TODO
//Eliminate the FE:interface when the point is at the centroid
// The first string is the class name of the object to be created,
// the second one is the name of the module as it should be referred
// in the input file (the Makefile defines MODULE_NAME, which can be used here).
TIBER_MODULE(Boltzmann, MODULE_NAME)
using namespace std;


Boltzmann*
Boltzmann::_this = NULL;



Boltzmann::Boltzmann(const ModelOptions& options) :
  SimulationInterface(options),
  dim(0),
  is_gray_solved(false),
  is_fourier_solved(false)
{
  // there's nothing to be done
}


Boltzmann*
Boltzmann::create(const ModelOptions& options)
{
  // we could use the options to create different implementations
  // or something like that.
  return new Boltzmann(options);
}



void
Boltzmann::do_init(void)
{
  parse_options();

  //get_scaling().set_calc_mesh_units(100.0 * get_scaling().get_calc_mesh_units());

  dim = get_mesh().mesh_dimension();

  do_init_fourier();

  //Get gray options
  get_gray_options();
  
  if (is_gray)
    do_init_gray();

  //Create Global Domain
  const MeshBase& mesh = get_mesh();
  {
    MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

    for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;
      GlobalDomain.insert(elem);
    }
  }

  //Create node connection
  const unsigned int nn  = mesh.n_nodes();
  node_conn.resize(nn);
  {
    vector<unsigned short int> node_conn_local(node_conn.size());


    MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

    for ( ; el != end_el; ++el)
      for (unsigned int n = 0; n < (*el)->n_nodes(); n++)
	node_conn_local[(*el)->node(n)]++;

    node_conn = node_conn_local;
  }



}

void
Boltzmann::get_gray_options(void)
{

  ModelOptions gray_opt = get_options();
  const MeshBase& mesh = get_mesh();

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  is_gray = false;

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);
    HeatTransportModel* htm = mod.get_heat_transport_model();

    if (htm->get_type() == HeatTransportModel::Gray)
    {
      if (is_gray == false)
      {
	gray_opt = htm->get_options();
	is_gray = true;
      }
    }
  }

  if (is_gray)
  {
    myopts.theta_slices = gray_opt.get_option("theta_slices",0);
    myopts.phi_slices = gray_opt.get_option("phi_slices",0);
    myopts.max_error = gray_opt.get_option("max_error",1e-3);
    myopts.max_iter = gray_opt.get_option("max_iter",1);
    myopts.diffusive =  gray_opt.get_option("diffusive_walls",false);
    myopts.partitioning = gray_opt.get_option("partitioning","manual");
    myopts.threshold_value =  gray_opt.get_option("temp_threshold",0.0);
    myopts.phi_zero = gray_opt.get_option("phi_zero",0);
  }



}


void
Boltzmann::do_partition(void)
{

  TiberLinearSystem& system = static_cast<TiberLinearSystem&>(
  						      get_equation_systems().get_system("fourier"));

  const MeshBase& mesh = get_mesh();
  DofMap& dof_map =  system.get_dof_map();
  const unsigned int tvar = system.variable_number("T");
  FEType fe_type = dof_map.variable_type(tvar);

 //------------------------------DOMAIN PARTITIONING---------------------------
 MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
 const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

 for ( ; el != end_el ; ++el)
 {
   const Elem* elem = *el;

   if (myopts.partitioning.compare("automatic") == 0)
   {

     //Get temperature from the first Fourier temperature
     std::vector<Point>  p(1);
     p[0] = elem->centroid();
     vector<Point> points(1);
     FEInterface::inverse_map(dim, fe_type, elem, p, points);
     std::map<ID, std::vector<double> > values;
     std::vector<double> value_vec(1);
     values[LatticeTemp] = value_vec;
     get_solution(elem,values,points,false);
     double value = values[LatticeTemp][0];
     //---------------------------------------
   
     if (value>myopts.threshold_value)
       GrayDomain.insert(elem);
     else
       FourierDomain.insert(elem);
   }
   else
   {

     BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);
     HeatTransportModel* htm = mod.get_heat_transport_model();

     if (htm->get_type() == HeatTransportModel::Gray)
       GrayDomain.insert(elem);
     else
       FourierDomain.insert(elem);
   }

 }



 //BTE/Fourier boundary side
 {
   set<const Elem*>::iterator el_f = GrayDomain.begin();
     const set<const Elem*>::iterator end_el_f = GrayDomain.end();

     for ( ; el_f != end_el_f ; ++el_f)
     {
       const Elem* elem = *el_f;
       ID neighbor = elem->n_neighbors();

       for (ID k = 0; k < neighbor; k ++)
       {
	 const Elem* elem_n = elem->neighbor(k);

	 if (FourierDomain.count(elem_n))
	 {
	   const ElementSide elside(elem->top_parent(),k);
	   BoundarySide.insert(elside);
	 }
       }

     }
   }


   {
    set<const Elem*>::iterator el_f = FourierDomain.begin();
    const set<const Elem*>::iterator end_el_f = FourierDomain.end();

    for ( ; el_f != end_el_f ; ++el_f)
    {
      const Elem* elem = *el_f;
      ID neighbor = elem->n_neighbors();

      for (ID k = 0; k < neighbor; k ++)
      {
	const Elem* elem_n = elem->neighbor(k);

	if (GrayDomain.count(elem_n))
	{
	  const ElementSide elside(elem->top_parent(),k);
	  BoundarySide.insert(elside);
	}
      }

    }
   }





}


NumericVector<double>&
Boltzmann::do_get_solution_vector(void)
{

  TiberLinearSystem* system = TiberLinearSystem::create(get_equation_systems(),
							"fourier", get_solver_options());

  system->get_solution_vector().close();
  return system->get_solution_vector();

}



void
Boltzmann::do_init_fourier(void)
{




TiberLinearSystem* system = TiberLinearSystem::create(get_equation_systems(),
					     "fourier", get_solver_options());


  system->add_variable("T", FIRST);
  system->attach_assemble_function(assemble_fourier);
  system->init();

  thermal_flux_nodal.resize(3);
  for (ID i = 0; i<3; i++ )
  {
    thermal_flux_nodal[i] = (system->solution)->clone().release();
    thermal_flux_nodal[i]->zero();
    thermal_flux_nodal[i]->close();
  }

 initial_energy = (system->solution)->clone().release();

}

void
Boltzmann::do_init_gray(void)
{

 if (is_gray)
 {
   
   TiberLinearSystem* system = TiberLinearSystem::create(get_equation_systems(),
							 "gray", get_solver_options());
   
   system->add_variable("T",CONSTANT,MONOMIAL);
   system->attach_assemble_function(assemble_gray);
   system->init();
   
   //N.B.:Hereafter the dof are elemental
   gray_sys_number = system->number();
   //Initialize thermal flux
   
   thermal_flux.resize(3);
   for (ID i = 0; i<3; i++ )
   {
     thermal_flux[i] = (system->solution)->clone().release();
     thermal_flux[i]->zero();
     thermal_flux[i]->close();
   }
   
   AngInt.dim = dim;
   AngInt.theta_slices = myopts.theta_slices;
   AngInt.phi_slices = myopts.phi_slices;
   AngInt.phi_zero = myopts.phi_zero;
   
   if (myopts.custom)
     AngInt.compute_custom_direction(myopts.cd);
   else
     AngInt.compute_directions();
   
   //Get Gray Options
   //Initialize direction solution

   sol_dir.resize(AngInt.n_slices);
   for (ID k = 0; k<AngInt.n_slices ; k++ )
   {
     sol_dir[k] = (system->solution)->clone().release();
     sol_dir[k]->zero();
   }
   
   //Initialize Equilibrium Energy
   equilibrium_energy = (system->solution)->clone().release();
   equilibrium_energy->add(SimulationOptions::temperature); //Just put some dummy value
   equilibrium_energy->close();
   
   if  (SimulationOptions::verbose() > 2)
     AngInt.print_info();

 }

}


//-------------------------------------------------------------------------//
Boltzmann::~Boltzmann()
{
  //Release pointers
  
  if (is_gray)
  {
    for (ID i = 0; i<3; i++ )
      delete thermal_flux[i]; //elemental (for the gray model)
    
    for (ID k = 0; k<AngInt.n_slices ; k++ )
      delete sol_dir[k];
    
    delete  equilibrium_energy;
  }
  
  for (ID i = 0; i<3; i++ )
    delete thermal_flux_nodal[i];
 

}
void
Boltzmann::parse_options(void)
{

 const ModelOptions& options = SimulationInterface::get_options();

 opts.ms_error = options.get_option("ms_error",1e-3);
 opts.ms_iter = options.get_option("ms_iter" ,1);
 opts.automatic_partitioning = options.get_option("automatic_partitioning",false);
 opts.threshold_value = options.get_option("threshold_value",0.0);
 opts.fourier_guess = options.get_option("fourier_guess",true);
 opts.do_fourier = options.get_option("do_fourier" ,true);
 myopts.diffusive = options.get_option("diffusive_walls",true);
 myopts.hot_contact = options.get_option("hot_contact","Hot");
 myopts.cold_contact = options.get_option("cold_contact","Cold");

//  myopts.max_error =  options.get_option("max_error",1e-3);
//  myopts.max_iter =  options.get_option("max_iter",1);
//  myopts.theta_slices = options.get_option("theta_slices",0);
//  myopts.phi_slices =    options.get_option("phi_slices",0);

  //Custom Dir
  myopts.custom = options.get_option("custom_dir",0);

  int n_dir = options.get_option("n_dir",1);

  myopts.cd.resize(n_dir);
  //  Custom Direction
  for (ID n = 0; n<n_dir; n++)
  {

    std::string str = static_cast<ostringstream*>( &(ostringstream() << n+1) )->str();
    std::string opt_vec = "dir_" + str;
    std::vector<double> cd(3);
    cd[0] = 0;
    cd[1] = 0;
    cd[2] = 1;
    options.get_option(opt_vec,cd);

    myopts.cd[n](0) = cd[0];
    myopts.cd[n](1) = cd[1];
    myopts.cd[n](2) = cd[2];

  }


}


void
Boltzmann::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(LatticeTemp, REAL, NODES, "K");
  declare_solution(FourierTemp, REAL, NODES, "K");
  declare_solution(MaxTemp, REAL, GLOBAL, "K");
  declare_solution(ThermalFlux, VECTOR, NODES, "W/m^2");
  declare_solution(ThermCond, VECTOR, NODES, "W/m K");
  declare_solution(HeatSource, REAL, NODES, "W/m^3");
  // declare_solution(EffectiveKappa, REAL, NODES, "W/cm^3");
  declare_solution(SolDir,VECTOR,CELL, "W/m^3");
  declare_solution(Partition,REAL,CELL, "");
  declare_solution(DomainTest,REAL,CELL, "");
  declare_solution(GRAY,REAL,CELL, "");
  // we can define aliases (but the association name -> id
  // has to be surjective)
  //add_alias("Jq", ThermalFlux);
  //add_alias("LatticeTemp", temperature);
  //add_alias("thermal", LatticeTemp);
  //add_alias("thermal", ThermCond);
  //add_alias("thermal", ThermalFlux);
  //add_alias("thermal", HeatSource);
  //add_alias("thermal", Partition);
  //add_alias("thermal", FourierTemp);
  //add_alias("thermal", EffectiveKappa);
  //  add_alias("ThermCond", ThermalConductivity);
}



void
Boltzmann::solve_fourier(void)
{

  //cout<<endl;
  //cout<<"      FOURIER..."<<endl;
  EquationSystems& es = get_equation_systems();

  TiberLinearSystem& system_fourier =
    es.get_system<TiberLinearSystem>("fourier");

  system_fourier.set_options(get_solver_options());
  system_fourier.solve();

  //Write the nodal flux
  const MeshBase& mesh = get_mesh();

  DofMap& dof_map_fourier = system_fourier.get_dof_map();
  std::vector<unsigned int> dof_indices_fourier;
  const NumericVector<double>& solution_fourier = *(system_fourier.solution);

  const unsigned int var = system_fourier.variable_number("T");
  FEType fe_type = dof_map_fourier.variable_type(var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  //--------------------------------------------
  for (ID d = 0; d< 3; d++)
  {
    thermal_flux_nodal[d]->close();
    thermal_flux_nodal[d]->zero();
  }
   //-----------------------------------------------------

   //Fill the temperature and the nodal flux
  set<const Elem*>::iterator el = Domain.begin();
  const set<const Elem*>::iterator end_el = Domain.end();
  
  double max_T = 0.0;
  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    dof_map_fourier.dof_indices (elem, dof_indices_fourier);
    
    std::vector<Point> p(1);
    p[0]=elem->centroid();
    vector<Point> points(1);
    FEInterface::inverse_map(dim, fe_type, elem, p, points);
    fe->reinit(elem, &points);
    
    BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);
    mod.calculate(elem,elem->centroid());
    const RealTensor& kappa = mod.get_total_thermal_conductivity();
    
    RealGradient heat_flux(0);
    for (ID alpha = 0; alpha<dof_indices_fourier.size() ;alpha ++)
      heat_flux -= solution_fourier(dof_indices_fourier[alpha]) * (kappa * dphi[alpha][0]);
    


    for (int n = 0; n < elem->n_nodes(); n++)
      for (ID d = 0; d< 3; d++)
	thermal_flux_nodal[d]->add(dof_indices_fourier[n],heat_flux(d)/node_conn[elem->node(n)]);
    
  }
  
  for (ID d = 0; d< 3; d++)
    thermal_flux_nodal[d]->close();
 

}


void
Boltzmann::from_nodal_to_cell()
{

  //Fourier System
  TiberLinearSystem& system_fourier = get_equation_systems().get_system<TiberLinearSystem>("fourier");
  DofMap& dof_map_fourier = system_fourier.get_dof_map();
  std::vector<unsigned int> dof_indices_fourier;
  const NumericVector<double>& solution_fourier = *(system_fourier.solution);

  const unsigned int var = system_fourier.variable_number("T");
  FEType fe_type = dof_map_fourier.variable_type(var);
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  //Gray System
  TiberLinearSystem& system_gray = get_equation_systems().get_system<TiberLinearSystem>("gray");
  DofMap& dof_map_gray = system_gray.get_dof_map();
  std::vector<unsigned int> dof_indices_gray;
  const NumericVector<double>& solution_gray = *(system_gray.solution);

  set<const Elem*>::iterator el = Domain.begin();
  const set<const Elem*>::iterator end_el = Domain.end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    dof_map_fourier.dof_indices (elem, dof_indices_fourier);
    dof_map_gray.dof_indices (elem, dof_indices_gray);

    BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);
    mod.calculate(elem,elem->centroid());
    const RealTensor& kappa = mod.get_total_thermal_conductivity();

    //Compute temperature at the centroid. Maybe here we can compute value at the centroid by using the JXW.
    std::vector<Point> p(1);
    p[0]=elem->centroid();
    vector<Point> points(1);
    FEInterface::inverse_map(dim, fe_type, elem, p, points);
    //---------------------------------------------------------

    fe->reinit(elem, &points);
    Real T = 0.0;
    RealGradient heat_flux(0);
    for (ID alpha = 0; alpha<dof_indices_fourier.size() ;alpha ++)
    {
      T += solution_fourier(dof_indices_fourier[alpha]) *  phi[alpha][0];
      heat_flux -= solution_fourier(dof_indices_fourier[alpha]) * (kappa * dphi[alpha][0]);
    }
    //----------------------------------

    equilibrium_energy->set(dof_indices_gray[0],T);

    for(ID d = 0; d<3; d++)
      thermal_flux[d]->set(dof_indices_gray[0],heat_flux(d));

  }

}


void
Boltzmann::from_cell_to_nodal()
{
  const MeshBase& mesh = get_mesh();
  //Fourier System
  TiberLinearSystem& system_fourier = get_equation_systems().get_system<TiberLinearSystem>("fourier");
  DofMap& dof_map_fourier = system_fourier.get_dof_map();
  std::vector<unsigned int> dof_indices_fourier;
  //const NumericVector<double>& solution_fourier = *(system_fourier.solution);

  //Gray System
  TiberLinearSystem& system_gray = get_equation_systems().get_system<TiberLinearSystem>("gray");
  DofMap& dof_map_gray = system_gray.get_dof_map();
  std::vector<unsigned int> dof_indices_gray;
  // const NumericVector<double>& solution_gray = *(system_gray.solution);

  //--------------------------------------------

  system_fourier.solution->zero();
  for (ID d = 0; d< 3; d++)
    thermal_flux_nodal[d]->zero();
  //-----------------------------------------------------

  set<const Elem*>::iterator el = Domain.begin();
  const set<const Elem*>::iterator end_el = Domain.end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    dof_map_fourier.dof_indices (elem, dof_indices_fourier);
    dof_map_gray.dof_indices (elem, dof_indices_gray);

    double T = (*equilibrium_energy)(dof_indices_gray[0]);

    RealGradient heat_flux(0);
    for (ID d = 0; d<3; d++)
      heat_flux(d) = (*thermal_flux[d])(dof_indices_gray[0]);

    for (int n = 0; n < elem->n_nodes(); n++)
    {
      system_fourier.solution->add(dof_indices_fourier[n],T/node_conn[elem->node(n)]);

      for (ID d = 0; d<3; d++)
	thermal_flux_nodal[d]->add(dof_indices_fourier[n],heat_flux(d)/node_conn[elem->node(n)]);

    }

  }

}


//   //! Order the solution in correct mode
// void Boltzmann::build_elemental_results(const std::set<std::string>& variables,
// 				     std::vector<double>& results,
// 				     std::vector<std::string>& legend)
// {

//   EquationSystems& es = get_equation_systems();
//   const MeshBase& mesh = get_mesh();

//   const unsigned int nn  = mesh.n_active_elem();
//   const unsigned int dim = mesh.mesh_dimension();
//   legend.resize(variables.size());

//   unsigned int n_vars = 0;
//   int par = -1;
//   if (variables.count("partial"))
//   {
//     par = n_vars;
//     for (ID k = 0;k<  AngInt.n_slices;k++)
//     {

//       std::string label;
//       std::ostringstream i_str;
//       i_str << "energy_density" << k;
//       cout<<i_str.str()<<endl;
//       legend.push_back(i_str.str());
//       n_vars++;
//     }
//   }

//   //if (variables.count("ThermalFlux"))
//   // {
    

//   //}
  

//   results.resize(nn * n_vars,0.0);
//   legend.resize(n_vars);

//   TiberLinearSystem& my_system =
//     es.get_system<TiberLinearSystem>("gray");

//   const unsigned int  var = my_system.variable_number("T");

//   DofMap& dof_map =  my_system.get_dof_map();
//   std::vector<unsigned int> dof_indices;

//   FEType fe_type = dof_map.variable_type(var);
//   MeshBase::const_element_iterator it =    mesh.active_local_elements_begin();
//   const MeshBase::const_element_iterator end =     mesh.active_local_elements_end();

//   unsigned int elem_number = 0;
//   for ( ; it != end; ++it)
//   {
//     const Elem* elem = *it;

//     dof_map.dof_indices (elem, dof_indices);

//     unsigned int id = n_vars * elem_number;

//     if (par != -1)
//       for (ID k = 0;k<  AngInt.n_slices;k++)
// 	results[id+par+k] = (*sol_dir[k])(dof_indices[0]);
      

//     elem_number++;
//   }

//  results.resize(elem_number * n_vars);

// }

void
Boltzmann::solve_gray(void)
{


 EquationSystems& es = get_equation_systems();

  TiberLinearSystem& system =
    es.get_system<TiberLinearSystem>("gray");

  system.set_options(get_solver_options());


  //Fill the directional results
  for (ID k = 0; k<AngInt.n_slices ; k++ )
  {
    sol_dir[k]->zero();
    sol_dir[k]->add(*equilibrium_energy);
  }
  //-----Fill the boundary value (only for the outer boundary being in the gray domain------------------
  SD.clear();
  {

    SimulationEnvironment& se = get_environment();
    set<const Elem*>::iterator el = Domain.begin();
    const set<const Elem*>::iterator end_el = Domain.end();

    for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;
      ID dof = elem->dof_number(gray_sys_number,0,0);
      for (ID ns = 0; ns < elem->n_sides(); ns++)
      {
        const ElementSide elside(elem->top_parent(),ns);
        if (se.is_outer_boundary(elside))
        {
	  SD[elside].resize(AngInt.n_slices);
	  for (ID k = 0; k<AngInt.n_slices; k++ )
            SD[elside][k] = (*sol_dir[k])(dof);

	}
      }
    }
  }


  // if  (SimulationOptions::verbose() > 2)
  //  AngInt.print_info();

  double old_energy_norm = 0.0;
  equilibrium_energy->close();
  double energy_norm = equilibrium_energy->l2_norm();
  double err = 0.0;
  double E_err = 1.0;
  double J_err = 1.0;
  ID iter = 0;

  cout<<endl;
  cout<<"      GRAY..."<<endl;

  solid_angle_iter = 0;
  while (E_err > myopts.max_error & iter < myopts.max_iter)
  {

    for (ID k = 0; k < AngInt.n_slices; k++ )
    {

        solid_angle_iter = k;
        // vec_spec = AngInt.spec[k];
         d_omega = AngInt.d_omega[k];
         IntDir = AngInt.directions[k];
         dir = AngInt.dir[k];

       if  (SimulationOptions::verbose() > 2)
	AngInt.print_info(k);

       //system.init();
      //cout<<k<<endl;
      system.solve();

      sol_dir[k]->zero();
      sol_dir[k]->add(*(system.solution));

      (system.solution)->zero();
    }

    //----UPDATE BOUNDARY DATA-----------------------
    SideData::iterator it(SD.begin());
    SideData::const_iterator it_end(SD.end());

    for ( ; it != it_end; ++it)
    {
      for (ID k = 0; k<AngInt.n_slices ; k++ )
      {
	ID dof = ((it->first).elem())->dof_number(gray_sys_number,0,0);
      	double value =  (*sol_dir[k])(dof);
	SD[it->first][k] = value;
      }
    }

    //----UPDATE EQUILIBRIUM ENERGY FOR THE GRAY DOMAIN----------------------------
    double max_T = 0.0;
    set<const Elem*>::iterator el = Domain.begin();
    const set<const Elem*>::iterator end_el = Domain.end();
    for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;
      ID dof = elem->dof_number(gray_sys_number,0,0);
      equilibrium_energy->set(dof,0.0);
      for (ID k = 0; k<AngInt.n_slices ; k++ )
      {
	double value = (*sol_dir[k])(dof) * AngInt.d_omega[k]/AngInt.total_angle;
	equilibrium_energy->add(dof,value);
      }
      if (max_T < (*equilibrium_energy)(dof))
	max_T = (*equilibrium_energy)(dof);

    }
    old_energy_norm = energy_norm;
    equilibrium_energy->close();
    energy_norm = equilibrium_energy->l2_norm();
    //-----------------------------------------------------------------------

    //---------------------------------------------
    E_err = abs(energy_norm - old_energy_norm)/max(energy_norm,old_energy_norm);

    //  cout<<"MAX TEMP:  "<<max_T<<" K"<<endl;


   //----COMPUTE THERMAL FLUX--------------

    ofstream myfile;

    myfile.open ("energy.dat",ios_base::app);
    myfile << energy_norm <<"\n";
    myfile.close();

    //if  (SimulationOptions::verbose() > 1)
    //{
    //  double kappa = compute_effective_thermal_conductivity();
    //  cout<<"Effettive thermal conductivity:    "<<kappa* 100.0<<" W/(m K)"<<endl;
    //}
    //J_err = energy_conservation_check_traditional();
    ////J_err = energy_conservation_check();
    //--------------------------------------------------------------

    cout<<endl;

    iter +=1;
    if  (SimulationOptions::verbose() > 1)
      cout<<"          Iter: "<<iter<<" E_err: "<<E_err<<endl;


  }
  //COMPUTE THERMAL FLUX
          {

            for (ID i = 0; i<3 ; i++ )
            {
              thermal_flux[i]->close();
              thermal_flux[i]->zero();
            }

            set<const Elem*>::iterator el = Domain.begin();
            const set<const Elem*>::iterator end_el = Domain.end();
            for ( ; el != end_el ; ++el)
            {
              const Elem* elem = *el;
              ID dof = elem->dof_number(gray_sys_number,0,0);
              BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);
              mod.calculate(elem,elem->centroid());
              double vg = mod.get_sound_velocity();
              double cg = mod.get_heat_capacity();


              for (ID k = 0; k<AngInt.n_slices ; k++ )
              {
                double value = (*sol_dir[k])(dof) * cg * vg / AngInt.total_angle;

                for (ID i = 0; i<3 ; i++ )
                {
                  double flux = value * AngInt.directions[k](i);
                  thermal_flux[i]->add(dof,flux);

                }
              }
            }
          }

  cout<<endl;
  cout<<"      ...GRAY"<<endl;
  cout<<endl;

}

double
Boltzmann::energy_conservation_check()
{

  SimulationEnvironment& se = get_environment();

  //Fourier System
  EquationSystems& es = get_equation_systems();
  TiberLinearSystem& system =
    es.get_system<TiberLinearSystem>("fourier");
  DofMap& dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;
  //-----------------------------------------------

  FEType fe_type(FIRST,LAGRANGE);
  AutoPtr<FEBase>  fe(build_finite_element(dim,fe_type,true));

  QGauss qrule(dim, CONSTANT);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<RealGradient> >&  dphi_rstf = fe->get_dphi();
  const std::vector<Real>& JxW = fe->get_JxW();

  double check = 0.0;

  set<const Elem*>::iterator it = Domain.begin();
  const set<const Elem*>::iterator end = Domain.end();

  for ( ; it != end; ++it)
  {

   //  const Elem* elem = *it;
//     dof_map.dof_indices(elem, dof_indices);

//     //-------------------------PowerDissipated------------------------------------------------
//     bool has_node = false;
//     const ID num_sides = elem->n_sides();
//     for (ID ns = 0; ns<num_sides; ns++)
//     {
//        const ElementSide elside(elem->top_parent(),ns);
//        if (is_on_any_boundary(elside))
// 	 has_node = true;
//     }

//     if (!has_node)
//       continue;

//     fe->reinit(elem);

//     for (ID d = 0; d<dim; d++)
//       for (ID n = 0; n < elem->n_nodes() ;n ++)
// 	check += JxW[0] * (*thermal_flux_nodal[d])(dof_indices[0]) * dphi_rstf[n][0](d);

    //----------------------------------------------------------------------------------------

   }

  if (SimulationOptions::verbose() > 1)
    std::cout<<"Energy conservation: "<<check<<std::endl;

  return check;

}

double
Boltzmann::energy_conservation_check_traditional()
{

  //Gray System
  EquationSystems& es = get_equation_systems();
  TiberLinearSystem& system =
    es.get_system<TiberLinearSystem>("gray");
  DofMap& dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;
  //-----------------------------------------------


  const unsigned int tvar = system.variable_number("T");
  FEType fe_type = dof_map.variable_type(tvar);

  //------------BULK----------
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim,FIFTH);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  //--------------------------


  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  QGauss qrule_face(dim-1,CONSTANT);
  fe_face->attach_quadrature_rule(&qrule_face);

  // const std::vector<Point>& q_point = fe_face->get_xyz();
  //const std::vector<std::vector<RealGradient> >&  dphi = fe_face->get_dphi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& normal = fe_face->get_normals();

  double power_dissipated= 0.0;
  double check_abs = 0.0;

  set<const Elem*>::iterator it = Domain.begin();
  const set<const Elem*>::iterator end = Domain.end();

  Real total_heat_source = 0.0;
  for ( ; it != end; ++it)
  {

    const Elem* elem = *it;
    dof_map.dof_indices(elem, dof_indices);

    fe->reinit(elem);

    BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);

    //Energy emitted
    for (ID qp = 0; qp <  qrule.n_points(); qp++)
    {
      mod.calculate(elem,q_point[qp]);
      Real H = mod.get_total_heat_source();
      total_heat_source += H * JxW[qp];
    }
    //-------------------------PowerDissipated------------------------------------------------

    for (ID ns = 0; ns<elem->n_sides(); ns++)
    {
      const ElementSide elside(elem->top_parent(),ns);
      if (is_on_any_boundary(elside))
      {
	fe_face->reinit(elem,ns);

	for (ID d = 0; d<3; d++)
	{
          double value = JxW_face[0] *(*thermal_flux[d])(dof_indices[0]) * normal[0](d);
	  power_dissipated += value;
          check_abs += abs(value);
	}
      }

    }
  }

  double error = 0.0;
  if (total_heat_source>0)
    error = std::abs(1.0 - power_dissipated/total_heat_source) ;
  else
    error = 2.0 * std::abs(power_dissipated/check_abs);

  if (SimulationOptions::verbose() > 1)
  {
    std::cout<<"Power Emitted: "<<total_heat_source<<" W"<<std::endl;
    std::cout<<"Power Dissipated: "<<power_dissipated<<" W"<<std::endl;
  }

  return error;

}

double
Boltzmann::compute_power_dissipated()
{

  //Gray System
  EquationSystems& es = get_equation_systems();
  TiberLinearSystem& system =
    es.get_system<TiberLinearSystem>("fourier");

  const NumericVector<double>& solution = *(system.solution);

  DofMap& dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;
  //-----------------------------------------------


  const unsigned int tvar = system.variable_number("T");
  FEType fe_type = dof_map.variable_type(tvar);

  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  QGauss qrule_face(dim-1,FIFTH);
  fe_face->attach_quadrature_rule(&qrule_face);
  
  const std::vector<Point>& q_point_face = fe_face->get_xyz();
  const std::vector<std::vector<RealGradient> >&  dphi = fe_face->get_dphi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& normal = fe_face->get_normals();
  
  double power_dissipated = 0.0;

  set<const Elem*>::iterator it = Domain.begin();
  const set<const Elem*>::iterator end = Domain.end();

 
  for ( ; it != end; ++it)
  {

    const Elem* elem = *it;
    dof_map.dof_indices(elem, dof_indices);
    
    BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);
   
    const RealTensor& kappa = mod.get_total_thermal_conductivity();
 
    for (ID ns = 0; ns<elem->n_sides(); ns++)
    {
      const ElementSide elside(elem->top_parent(),ns);
      if (is_on_any_boundary(elside))
      {

	//BoltzmannBoundaryModel* mod_b =
	// get_interface_model<BoltzmannBoundaryModel>(elem, ns);
	
	//if (mod_b != NULL)
	//{
	fe_face->reinit(elem,ns);
	for (ID qp = 0; qp <  qrule_face.n_points(); qp++)
	{ 
	
	  for (ID alpha = 0; alpha<dof_indices.size() ;alpha ++)
	    power_dissipated -= JxW_face[qp] * solution(dof_indices[alpha]) * ((kappa * dphi[alpha][qp]) * normal[qp]);  
	  
	}
       
      }
    }
    
  }

  return power_dissipated;

}

double 
Boltzmann::compute_porosity()
{

  EquationSystems& es = get_equation_systems();
  TiberLinearSystem& system =
    es.get_system<TiberLinearSystem>("fourier");
  DofMap& dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;
  //-----------------------------------------------

  const unsigned int tvar = system.variable_number("T");
  FEType fe_type = dof_map.variable_type(tvar);
  
  //------------BULK----------
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim,FIFTH);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  //--------------------------

  set<const Elem*>::iterator it = Domain.begin();
  const set<const Elem*>::iterator end = Domain.end();

  Real total_volume = 0.0;
  for ( ; it != end; ++it)
  {

    const Elem* elem = *it;
    dof_map.dof_indices(elem, dof_indices);

    fe->reinit(elem);

    BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);

    for (ID qp = 0; qp <  qrule.n_points(); qp++)
    {
      total_volume += JxW[qp];
    }
  }

  total_volume *= get_scaling().get_calc_mesh_units();
  total_volume *= get_scaling().get_calc_mesh_units();
  total_volume *= get_scaling().get_calc_mesh_units();

  cout<<"Total Volume: "<<total_volume<<" m3";



 return  total_volume;

}


double
Boltzmann::compute_effective_thermal_conductivity()
{
 EquationSystems& es = get_equation_systems();


 if (is_gray)
 {

  //-----------------GRAY---------------
   TiberLinearSystem& gray_system =
     es.get_system<TiberLinearSystem>("gray");
   const NumericVector<Number>& gray_solution = gray_system.get_solution_vector();
   const DofMap& dof_map_gray= gray_system.get_dof_map();
   vector<unsigned int> dof_indices_gray;  
   const unsigned int tvar_gray = gray_system.variable_number("T");
   FEType fe_type_gray = dof_map_gray.variable_type(tvar_gray);

   //--------------surface-----------------
   AutoPtr<FEBase> fe_face_gray(build_finite_element(dim, fe_type_gray, true));
   QGauss qrule_face_gray(dim-1,CONSTANT);
   fe_face_gray->attach_quadrature_rule(&qrule_face_gray);
   const std::vector<Real>& JxW_face_gray = fe_face_gray->get_JxW();
   const std::vector<Point>& normal_gray = fe_face_gray->get_normals();
   //-----------------------------------------------------------------------------
  }
   //Fourier System
   TiberLinearSystem& system =
    es.get_system<TiberLinearSystem>("fourier");

  const NumericVector<double>& solution = *(system.solution);
  DofMap& dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;
  const unsigned int tvar = system.variable_number("T");
  FEType fe_type = dof_map.variable_type(tvar);

   //------------BULK----------
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim,FIFTH);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();

  //--------------SURFACE-----------------
  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  QGauss qrule_face(dim-1,FIFTH);
  fe_face->attach_quadrature_rule(&qrule_face);

  const std::vector<Point>& q_point_face = fe_face->get_xyz();
  const std::vector<std::vector<RealGradient> >&  dphi = fe_face->get_dphi();
  const std::vector<std::vector<double> >&  phi = fe_face->get_phi();

  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& normal = fe_face->get_normals();
  
  double abs_J = 0.0;
  double area = 0.0;
  double DeltaT = 0.0;


  double Phot = 0.0;
  double Pcold = 0.0;
  double Ahot = 0.0;
  double Acold = 0.0;
  double Thot = 0.0;
  double Tcold = 0.0;



  Point Point_hot(0);
  Point Point_cold(0);


  set<const Elem*>::iterator it = Domain.begin();
  const set<const Elem*>::iterator end = Domain.end();
  
  double total_heat_source = 0.0;
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    dof_map.dof_indices(elem, dof_indices);

    fe->reinit(elem);
    
    BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);
    const RealTensor& kappa = mod.get_total_thermal_conductivity();
    
    //Power emitted
    for (ID qp = 0; qp <  qrule.n_points(); qp++)
    {
      mod.calculate(elem,q_point[qp]);
      Real H = mod.get_total_heat_source();
      total_heat_source += H * JxW[qp];
    }

    for (ID ns = 0; ns<elem->n_sides(); ns++)
    {
      const ElementSide elside(elem->top_parent(),ns);
      if (is_on_any_boundary(elside))
      {
	BoltzmannBoundaryModel* mod_b =
	  get_interface_model<BoltzmannBoundaryModel>(elem, ns);
	
	fe_face->reinit(elem,ns);

	if (mod_b != NULL)
	{

	  if (mod_b->get_name() == myopts.cold_contact)
	  {

	    //For heat_reservoir contact we take the nominal value (the imposed one)
	    if (mod_b->get_type()=="boltzmann_bnd_heat_reservoir")
	    {
	      double a,b,c;
	      mod_b->get_coefficients(a,b,c);
	      Tcold = c;
	    }

	    //For periodic boundary condtion we take the naturally developed temperature
	    if (mod_b->get_type()=="boltzmann_bnd_periodic")
	    {
	      ID dof = elem->dof_number(gray_sys_number,0,0);
	      Tcold = (*equilibrium_energy)(dof);
	    }

	    Point_cold = q_point_face[0];
	    
	    for (ID qp = 0; qp <  qrule_face.n_points(); qp++)
	    { 
	      Acold += JxW_face[qp];
	      for (ID alpha = 0; alpha<dof_indices.size() ;alpha ++)
	      {
		double Px = (*thermal_flux_nodal[0])(dof_indices[alpha]) * normal[qp](0);
		double Py = (*thermal_flux_nodal[1])(dof_indices[alpha]) * normal[qp](1);
		double Pz = (*thermal_flux_nodal[2])(dof_indices[alpha]) * normal[qp](2);
                Pcold +=  JxW_face[qp] * phi[alpha][qp] * (Px + Py + Pz);
	      } 
	    }
	  }

	  if (mod_b->get_name() == myopts.hot_contact)
	  {
	    //For heat_reservoir contact we take the nominal value (the imposed one)
	    if (mod_b->get_type()=="boltzmann_bnd_heat_reservoir")
	    {
	      double a,b,c;
	      mod_b->get_coefficients(a,b,c);
	      Thot = c;
	    }

	    //For periodic boundary condtion we take the naturally developed temperature
	    if (mod_b->get_type()=="boltzmann_bnd_periodic")
	    {
	      ID dof = elem->dof_number(gray_sys_number,0,0);
	      Thot = (*equilibrium_energy)(dof);
	    }


	    Point_hot = q_point_face[0];

	    for (ID qp = 0; qp <  qrule_face.n_points(); qp++)
	    {
	      Ahot += JxW_face[qp];
              
	      for (ID alpha = 0; alpha<dof_indices.size() ;alpha ++)
	      {

		double Px = (*thermal_flux_nodal[0])(dof_indices[alpha]) * normal[qp](0);
		double Py = (*thermal_flux_nodal[1])(dof_indices[alpha]) * normal[qp](1);
		double Pz = (*thermal_flux_nodal[2])(dof_indices[alpha]) * normal[qp](2);
	        Phot +=  JxW_face[qp] * phi[alpha][qp] * (Px + Py + Pz);
	      } 
	    }
	  }
   
	}
      }
    }
  }


  //Guess of the direction of the heat transport
  //It simply takes the max coordinate displacement among contacts
  Point dist = Point_hot - Point_cold;

  double max_dist = - 1.0;

  for (ID d = 0; d<3; d++)
   if (abs(dist(d))>max_dist)
       max_dist = abs(dist(d));
  //-------------------------------------------

  //cout<<"Pcold_gray: "<<Pcold_gray<<" W"<<endl;  //! Order the solution in correct mode
  //cout<<"Phot_gray: "<<Phot_gray<<" W"<<endl;
  cout<<"Pcold: "<<Pcold<<" W"<<endl;  //! Order the solution in correct mode
  cout<<"Phot: "<<Phot<<" W"<<endl;
  cout<<"Tcold: "<<Tcold<<" K"<<endl;  //! Order the solution in correct mode
  cout<<"Thot: "<<Thot<<" K"<<endl;
  cout<<"Power emitted: "<<total_heat_source*1e-6<<" W"<<endl;

  cout<<abs(Thot - Tcold)<<endl;
 // cout<< max_dist * get_scaling().get_calc_mesh_units()<<endl;
  //cout<<Ahot<<endl;

  double effective_kappa = Pcold/Ahot/(abs(Thot - Tcold)) * max_dist * get_scaling().get_calc_mesh_units();

  ofstream myfile;

  myfile.open ("kappa.dat",ios_base::app);
  myfile << effective_kappa <<"\n";
  myfile.close();

  return effective_kappa;

}
double
Boltzmann::compute_power_emitted()
{

  //Gray System
  EquationSystems& es = get_equation_systems();
  TiberLinearSystem& system =
    es.get_system<TiberLinearSystem>("fourier");
  DofMap& dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;
  //-----------------------------------------------


  const unsigned int tvar = system.variable_number("T");
  FEType fe_type = dof_map.variable_type(tvar);

  //------------BULK----------
  AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
  QGauss qrule(dim,FIFTH);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  //--------------------------


  AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
  QGauss qrule_face(dim-1,CONSTANT);
  fe_face->attach_quadrature_rule(&qrule_face);

  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& normal = fe_face->get_normals();

  double power_dissipated= 0.0;
  double check_abs = 0.0;

  set<const Elem*>::iterator it = Domain.begin();
  const set<const Elem*>::iterator end = Domain.end();

  Real total_heat_source = 0.0;
  for ( ; it != end; ++it)
  {

    const Elem* elem = *it;
    dof_map.dof_indices(elem, dof_indices);

    fe->reinit(elem);

    BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);

    //Energy emitted
    for (ID qp = 0; qp <  qrule.n_points(); qp++)
    {
      mod.calculate(elem,q_point[qp]);
      Real H = mod.get_total_heat_source();
      total_heat_source += H * JxW[qp];
    }
    //-------------------------PowerDissipated------------------------------------------------

 //    for (ID ns = 0; ns<elem->n_sides(); ns++)
//     {
//       const ElementSide elside(elem->top_parent(),ns);
//       if (is_on_any_boundary(elside))
//       {
// 	fe_face->reinit(elem,ns);

// 	for (ID d = 0; d<dim; d++)
// 	{ solution(dof_indices[alpha]) * (dphi[alpha][qp] * normal[qp])
//           double value = JxW_face[0] *(*thermal_flux[d])(dof_indices[0]) * normal[0](d);
// 	  power_dissipated += value;
//           check_abs += abs(value);
// 	}
//       }

//     }
  }

  //double error = 0.0;
  //if (total_heat_source>0)
  //  error = std::abs(1.0 - power_dissipated/total_heat_source) ;
  //else
  //  error = 2.0 * std::abs(power_dissipated/check_abs);

  // if (SimulationOptions::verbose() > 1)
  // {
  //  std::cout<<"Power Emitted: "<<total_heat_source<<" W"<<std::endl;
  //  std::cout<<"Power Dissipated: "<<power_dissipated<<" W"<<std::endl;
  // }

  return  total_heat_source;

}

// double
// Boltzmann::compute_power_emitted()
// {

//   //Gray System
//   EquationSystems& es = get_equation_systems();
//   TiberLinearSystem& system =
//     es.get_system<TiberLinearSystem>("f");
//   DofMap& dof_map = system.get_dof_map();
//   std::vector<unsigned int> dof_indices;
//   //-----------------------------------------------


//   const unsigned int tvar = system.variable_number("T");
//   FEType fe_type = dof_map.variable_type(tvar);

//   //------------BULK----------
//   AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
//   QGauss qrule(dim,FIFTH);
//   fe->attach_quadrature_rule(&qrule);
//   const std::vector<Real>& JxW = fe->get_JxW();
//   const std::vector<Point>& q_point = fe->get_xyz();
//   //--------------------------


//   AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
//   QGauss qrule_face(dim-1,CONSTANT);
//   fe_face->attach_quadrature_rule(&qrule_face);

//   // const std::vector<Point>& q_point = fe_face->get_xyz();
//   //const std::vector<std::vector<RealGradient> >&  dphi = fe_face->get_dphi();
//   const std::vector<Real>& JxW_face = fe_face->get_JxW();
//   const std::vector<Point>& normal = fe_face->get_normals();

//   double power_dissipated= 0.0;
//   double check_abs = 0.0;

//   set<const Elem*>::iterator it = Domain.begin();
//   const set<const Elem*>::iterator end = Domain.end();

//   Real total_heat_source = 0.0;
//   for ( ; it != end; ++it)
//   {

//     const Elem* elem = *it;
//     dof_map.dof_indices(elem, dof_indices);

//     fe->reinit(elem);

//     BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);

//     for (ID qp = 0; qp <  qrule.n_points(); qp++)
//     {
//       mod.calculate(elem,q_point[qp]);
//       Real H = mod.get_total_heat_source();
//       total_heat_source += H * JxW[qp];
//     }
//     //-------------------------PowerDissipated------------------------------------------------

//     for (ID ns = 0; ns<elem->n_sides(); ns++)
//     {
//       const ElementSide elside(elem->top_parent(),ns);
//       if (is_on_any_boundary(elside))
//       {
// 	fe_face->reinit(elem,ns);

// 	for (ID d = 0; d<dim; d++)
// 	{
//           double value = JxW_face[0] *(*thermal_flux[d])(dof_indices[0]) * normal[0](d);
// 	  power_dissipated += value;
//           check_abs += abs(value);
// 	}
//       }

//     }
//   }

//   double error = 0.0;
//   if (total_heat_source>0)
//     error = std::abs(1.0 - power_dissipated/total_heat_source) ;
//   else
//     error = 2.0 * std::abs(power_dissipated/check_abs);

//   if (SimulationOptions::verbose() > 1)
//   {
//     std::cout<<"Power Emitted: "<<total_heat_source<<" W"<<std::endl;
//     std::cout<<"Power Dissipated: "<<power_dissipated<<" W"<<std::endl;
//   }

//   return error;

// }






void
Boltzmann::do_solve(void)
{
  _this = this;

  is_fourier_solved = false;
  is_gray_solved = false;
  first_guess = true;

  if (opts.fourier_guess)
  {
    Domain = GlobalDomain;
    cout<<"FOURIER over ALL...";
    solve_fourier();
    cout<<"...done."<<endl;
    is_fourier_solved = true;
    first_guess = false;
    EquationSystems& es = get_equation_systems();
    TiberLinearSystem& system =
      es.get_system<TiberLinearSystem>("fourier");
    initial_energy = (system.solution)->clone().release();
  }




  if (is_gray)
  {
    //! Partition of the simulation domain
    do_partition();
    
    bool do_multiscale = FourierDomain.size() > 0 && opts.do_fourier;
    //Only here the elemental results are filled
    if (GrayDomain.size() > 0 && opts.ms_iter > 0 )
    {
      if (do_multiscale)
      {
	cout<<"MULTISCALE LOOP...."<<endl;
	cout<<endl;
      }
      
      Domain = GlobalDomain;
      from_nodal_to_cell();
      
      double old_energy_norm = 0.0;
      equilibrium_energy->close();
      double energy_norm = equilibrium_energy->l2_norm();
      
      double err = 2.0 * opts.ms_error;
      ID iter = 0;
      
      while (err > opts.ms_error & iter < opts.ms_iter)
      {
	//-------Self Consistent Loop------------
	Domain = GrayDomain;
	
	if (iter>0) //The first guess for gray!
	{ is_gray_solved = false;
	  cout<<"FOURIER over GRAY"<<endl;
	  solve_fourier();
	  from_nodal_to_cell();
	}
	
	solve_gray();

	is_gray_solved = true;
	
	if (do_multiscale)
	{
	  cout<<"FOURIER over FOURIER"<<endl;
	  Domain = FourierDomain;
	  solve_fourier();
	  from_nodal_to_cell();
	}
	
	//---------------------------------------
	old_energy_norm = energy_norm;
	equilibrium_energy->close();
	energy_norm = equilibrium_energy->l2_norm();
	err = abs(energy_norm - old_energy_norm)/max(energy_norm,old_energy_norm);
	iter ++;

	Domain = GlobalDomain;
	from_cell_to_nodal();
	
	if (do_multiscale && iter> 0)
	  std::cout<<"      Iter: "<<iter<<" Error: "<<err<<std::endl;
	else
	  err  = opts.ms_error - 1e-6; //Quit soon
	
      }
      if (do_multiscale)
      {
	cout<<endl;
	cout<<"...MULTISCALE LOOP."<<endl;
      }
      
    }
  } //if is_gray

if  (SimulationOptions::verbose() > 1)
 {
  double kappa = compute_effective_thermal_conductivity();
  cout<<"Effettive thermal conductivity:    "<<kappa<<" W/(m K)"<<endl;
}
 //if  (SimulationOptions::verbose() > 1)
 //{
   //double pe = compute_power_emitted();
   //double pd = compute_power_dissipated();

   //cout<<"Power Emitted:    "<<pe<<" W"<<endl;
  // cout<<"Power Dissipated: "<<endl;


  // double  porosity = compute_porosity();


  //cout<<"Porosity:    "<<compute_porosity()<<endl;

 //}
  //J_err = energy_conservation_check();
}


void
Boltzmann::do_print_info(void)
{
  // Messages::info("THERMONEO");
}


PhysicalModel*
Boltzmann::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{

  return  BoltzmannModel::create(mat,options);

}



//! We need to create boundary condition model
PhysicalModel*
Boltzmann::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{
  return BoltzmannBoundaryModel::create(boundary, options);
}


//PhysicalModel*
//Boltzmann::create_boundary_model(const ModelOptions& options,
//   const Material* material_A, const Material* material_B) const
//{
//  return BoltzmannBoundaryModel::create(options);
//}


void
Boltzmann::get_solution_secure(std::map<ID, std::vector<double> >& values)
{
  if (values.count(MaxTemp))
  {
    TiberLinearSystem* system =
      &get_equation_systems().get_system<TiberLinearSystem>("fourier");
    const NumericVector<Number>& solution = system->get_solution_vector();
    double Tmax = solution.linfty_norm();
    values[MaxTemp] = std::vector<double>(1, Tmax);
  }
}



void
Boltzmann::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
   unsigned int np = p.size();

   TiberLinearSystem* system;
   system = &get_equation_systems().get_system<TiberLinearSystem>(
       "fourier");
   const NumericVector<Number>& solution = system->get_solution_vector();
   const DofMap& dof_map = system->get_dof_map();

   //Test Gray
 //   TiberLinearSystem* gray_system;
//    gray_system = &get_equation_systems().get_system<TiberLinearSystem>("gray");
//    const NumericVector<Number>& gray_solution = gray_system->get_solution_vector();
//    const DofMap& dof_map_gray= gray_system->get_dof_map();
//    vector<unsigned int> dof_indices_gray;

   //-------------
   const unsigned int u_var = system->variable_number("T");

   FEType fe_type = system->variable_type(u_var);
   AutoPtr<FEBase> fe(build_finite_element(dim, fe_type));

   vector<unsigned int> dof_indices;

   //element shape functions
   const vector<vector<Real> >& phi = fe->get_phi();
   const vector<vector<RealGradient> >& dphi = fe->get_dphi();
   const vector<Point>& real_pts = fe->get_xyz();

   ID subdomain = elem->subdomain_id();

   fe->reinit(elem, &p);

   dof_map.dof_indices(elem, dof_indices, u_var);

   const unsigned int n_dofs = dof_indices.size();
   BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);
   const RealTensor& kappa = mod.get_total_thermal_conductivity();

   for (unsigned int n = 0; n < np; n++)
    {
      mod.calculate(elem,real_pts[n]);
      
      if (values.count(LatticeTemp))
      {
	double T  = 0.0;
	for (unsigned int i = 0; i < n_dofs; i++)
	  T += phi[i][n] * solution(dof_indices[i]);
	
	values[LatticeTemp][n] = T;
      }
   
      if (values.count(ThermalFlux))  
      {
	
	RealGradient heat_flux(0);
	for (ID i = 0; i < n_dofs; i++)
	  for (ID d = 0; d < 3; d++)
	    heat_flux(d) += phi[i][n] * (*thermal_flux_nodal[d])(dof_indices[i]);
	
	
	for (ID d = 0; d < dim; d++)
	  values[ThermalFlux][d + 3 * n] = heat_flux(d);
	
      }

      
      if (values.count(ThermCond))
      {
	const RealTensor& kappa = mod.get_total_thermal_conductivity();
	values[ThermCond][0 + 3 * n] = kappa(0,0);
	values[ThermCond][1 + 3 * n] = kappa(1,1);
	values[ThermCond][2 + 3 * n] = kappa(2,2);
      }
      
      if (values.count(HeatSource))
      {
	
	Real H = mod.get_total_heat_source();
	values[HeatSource][n] = H;
	
      }
      
    }

   if (is_gray)
   {
     if (values.count(Partition))
     {
       double value = 0;
       if (GrayDomain.count(elem))
	 value = 1;
       
       values[Partition][0] = value;
     }
   }

   


//      if (values.count(GRAY)||
// 	 values.count(thermal)  )
//      {
//        dof_map_gray.dof_indices (elem, dof_indices_gray);

//        double H = (*equilibrium_energy)(dof_indices_gray[0]);

//        values[GRAY][0] = H;

//      }
   //if (values.count(FourierTemp))
      // {
      //  double T  = 0.0;
      //  for (unsigned int i = 0; i < n_dofs; i++)
      // 	 T += phi[i][n] * (*initial_energy)(dof_indices[i]);
      //  values[FourierTemp][n] = T;
      // }


}


void
Boltzmann::clear_system(const std::string& system_name)
{
  TiberLinearSystem& system = static_cast<TiberLinearSystem&>(
		      get_equation_systems().get_system(system_name));


  const MeshBase& mesh = get_mesh();

  vector<unsigned int> dof_indices;

  DenseMatrix<Number> Ke;

  DenseVector<Number> Fe;

  DofMap& dof_map =  system.get_dof_map();
 //ClearSystem
  MeshBase::const_element_iterator       el   = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  for ( ; el != end_el ; ++el)
  {

    const Elem* elem = *el;
    dof_map.dof_indices(elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();
    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);
    Fe.zero();
    Ke.zero();
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
  }

}

void
Boltzmann::do_assemble_global(EquationSystems& es, const std::string& system_name)
{

}

void
Boltzmann::do_assemble_gray(EquationSystems& es, const std::string& system_name)
{


 TiberLinearSystem& system = static_cast<TiberLinearSystem&>(
		       get_equation_systems().get_system("gray"));

   const MeshBase& mesh = get_mesh();
   DofMap& dof_map =  system.get_dof_map();
   const unsigned int tvar = system.variable_number("T");
   FEType fe_type = dof_map.variable_type(tvar);


   // the volume finite element
   AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
   //QGauss qrule(dim, FIFTH);
   QGauss qrule(dim, CONSTANT);
   fe->attach_quadrature_rule(&qrule);

   const vector<Real>& JxW = fe->get_JxW();
   const vector<Point>& q_point = fe->get_xyz();
   const vector<vector<Real> >& phi = fe->get_phi();
   const vector<vector<RealGradient> >& dphi = fe->get_dphi();


   // the surface finite element
   AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
   //QGauss qface(dim - 1, SIXTH);
   QGauss qface(dim - 1, CONSTANT);
   fe_face->attach_quadrature_rule(&qface);

   const vector<Real>& JxW_face = fe_face->get_JxW();
   const vector<Point>& qface_point = fe_face->get_xyz();
   const vector<vector<Real> >&  phi_face = fe_face->get_phi();
   const vector<vector<RealGradient> >& dphi_face = fe->get_dphi();
   const vector<Point>& normal = fe_face->get_normals();

   vector<unsigned int> dof_indices;

   DenseMatrix<Number> Ke;
   DenseVector<Number> Fe;

   //Clear the system (TO FIND ANOTHER METHOD)
   clear_system("gray");

   set<const Elem*>::iterator el = Domain.begin();
   const set<const Elem*>::iterator end_el = Domain.end();

  // system.matrix->clear();
  // system.rhs->clear();

   ID el_numb = 0;
   SimulationEnvironment& se = get_environment();
   for ( ; el != end_el ; ++el)
   {

     el_numb ++;
     const Elem* elem = *el;
     fe->reinit(elem);

     //cout<<"EL NUMBER: "<<el_numb<<endl;


     dof_map.dof_indices (elem,dof_indices);
     //const unsigned int n_dofs = dof_indices.size();


     //Update list of neighbors
     ID neighbor = elem->n_neighbors();

     for (ID k = 0; k < neighbor; k ++)
     {
       const Elem* elem_n = elem->neighbor(k);
       if (elem_n != NULL)
       {
	 std::vector<unsigned int> dof_indices_n;
	 dof_map.dof_indices (elem_n,dof_indices_n);
	 dof_indices.push_back(dof_indices_n[0]);
       }
       else
        dof_indices.push_back(0);
     }

     const unsigned int n_dofs = dof_indices.size();
     Fe.zero();

     Ke.zero();
     Ke.resize(n_dofs,n_dofs);
     Fe.resize(n_dofs);


     const unsigned int num_sides = elem->n_sides();

     BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);
     mod.calculate(elem,elem->centroid());

     double tg = mod.get_relaxation_time();
     double vg = mod.get_sound_velocity();
     double e_0 = (*equilibrium_energy)(dof_indices[0]);
     Real cg = mod.get_heat_capacity();
     Real heat_source = mod.get_total_heat_source() * 1e-6;

     //Assembly
     Ke(0,0) = JxW[0];

     ID nb = 0;
     for (ID ns = 0; ns < elem->n_sides(); ns++)
     {

       fe_face->reinit(elem,ns);
       const ElementSide elside(elem->top_parent(),ns);
       double in = dir * normal[0];
       double value = vg * tg * (IntDir * normal[0])/d_omega * JxW_face[0];

       if (in<0.0)
       {

	 if (se.is_outer_boundary(elside) || (is_on_GF_boundary(elside) && is_fourier_solved))
	 {
	   Fe(0) -= get_boundary_value(elside,normal[0],dir) * value;
	   nb ++;
	 }
	 else
	    Ke(0,ns + 1) += value;
       }
       else
         Ke(0,0) += value;

     }

     //This term is a mess
     //if (nb>0)
      // Fe(0) /=nb;


    Fe(0) += (e_0 +  heat_source * tg / cg) * JxW[0];

    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices);

   }

   system.matrix->close();
   system.rhs->close();
   //system.matrix->print_matlab("K.m");

   //system.rhs->print_matlab("F.m");

}


double
Boltzmann::get_boundary_value(ElementSide elside, Point normal, Point dir)
{

  BoltzmannBoundaryModel* mod =
    get_interface_model<BoltzmannBoundaryModel>(elside.elem(), elside.side());


  SimulationEnvironment& se = get_environment();

  double value = 0.0;
  if (mod != NULL)
  {
    mod->calculate(elside.elem(), elside.side(), elside.elem()->centroid());

    if (mod->get_type() == "boltzmann_bnd_heat_reservoir")
    {
      double a, b, c;
      //c is the temperature
      mod->get_coefficients(a, b, c);
      value = c;
    }

    if (mod->get_type() == "boltzmann_bnd_periodic")
      {
         //Get the periodicity vector
         RealGradient periodicity = mod->get_periodicity();
         double DeltaT = mod->get_deltaT();

         //Get the point to search for
         AutoPtr< Elem > elem_face =  elside.elem()->build_side(elside.side());
         Point point_side = elem_face->point(0);

         //The 1 - delta term has been included in order to be sure that the point is in the element
         Point pt_to_search = (point_side + periodicity);


         //Get the element of that point
          SideData::iterator it(SD.begin());
          SideData::const_iterator it_end(SD.end());

         value = 0.0;
         for ( ; it != it_end; ++it)
            if ((it->first).elem() -> contains_point(pt_to_search))
            {

              ID local_dof = elside.elem()->dof_number(gray_sys_number,0,0);
              ID per_dof = (it->first).elem()->dof_number(gray_sys_number,0,0);
              double DeltaT = (*equilibrium_energy)(local_dof)-(*equilibrium_energy)(per_dof);

                //value = SD[it->first][solid_angle_iter] + DeltaT;
              value = SD[it->first][solid_angle_iter] + DeltaT;

            }
         //cout<<point_side<<endl;
        //cout<<pt_to_search<<endl;
         if (value == 0)
           cout<<"ERROR: no ELEMENT FOUND!!!"<<endl;



      }
  }
  else //Internal boundary or wall
  {
    //Internal boundary
    if (is_on_GF_boundary(elside))
    {

      ID dof = (elside.elem()->neighbor(elside.side()))->dof_number(gray_sys_number,0,0);
      value = (*equilibrium_energy)(dof);

    }
    else //Wall
    {

      //Diffusive
      if (myopts.diffusive)
      {
	ID dof = elside.elem()->dof_number(gray_sys_number,0,0);
	value = (*equilibrium_energy)(dof);

      }
      else
      {//SPECULAR

        ID dof = elside.elem()->dof_number(gray_sys_number,0,0);


        //Check which direction is the specular for thi side
        for (ID k = 0; k<AngInt.n_slices; k++)
        {
          Point check = (AngInt.dir[k] - dir)/(AngInt.dir[k] - dir).size() - normal;

          if (check.size()<1e-5)
            value = SD[elside][k];
      }

    }

   }
  }

  return value;


}

void
Boltzmann::do_assemble_fourier(EquationSystems& es, const std::string& system_name)
{
  TiberLinearSystem& system_fourier = static_cast<TiberLinearSystem&>(
								      get_equation_systems().get_system("fourier"));

   const MeshBase& mesh = get_mesh();

   DofMap& dof_map =  system_fourier.get_dof_map();

   const unsigned int tvar = system_fourier.variable_number("T");

   FEType fe_type = dof_map.variable_type(tvar);

   SimulationEnvironment& se = get_environment();
   // the volume finite element
   AutoPtr<FEBase> fe(build_finite_element(dim, fe_type, true));
   QGauss qrule(dim, FIFTH);
   fe->attach_quadrature_rule(&qrule);

   const vector<Real>& JxW = fe->get_JxW();
   const vector<Point>& q_point = fe->get_xyz();
   const vector<vector<Real> >& phi = fe->get_phi();
   const vector<vector<RealGradient> >& dphi = fe->get_dphi();


   // the surface finite element
   AutoPtr<FEBase> fe_face(build_finite_element(dim, fe_type, true));
   QGauss qface(dim - 1, SIXTH);
   fe_face->attach_quadrature_rule(&qface);

   const vector<Real>& JxW_face = fe_face->get_JxW();
   const vector<Point>& qface_point = fe_face->get_xyz();
   const vector<vector<Real> >&  phi_face = fe_face->get_phi();
   const vector<vector<RealGradient> >& dphi_face = fe->get_dphi();
   const vector<Point>& normal = fe_face->get_normals();

   vector<unsigned int> dof_indices;

   DenseMatrix<Number> Ke;
   DenseVector<Number> Fe;

   //Clear the system (TO FIND ANOTHER METHOD)
   clear_system("fourier");


   set<const Elem*>::iterator el = Domain.begin();
   const set<const Elem*>::iterator end_el = Domain.end();

   for ( ; el != end_el ; ++el)
   {
     const Elem* elem = *el;

     dof_map.dof_indices(elem, dof_indices);

     const unsigned int n_dofs = dof_indices.size();

     //resize the element matrix/rhs (does also zero them out)
     Ke.resize(n_dofs, n_dofs);

     Fe.resize(n_dofs);
     fe->reinit(elem);

     BoltzmannModel& mod = *get_bulk_model<BoltzmannModel>(elem);

     // loop over the quadrature points
     for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
     {

       mod.calculate(elem,q_point[qp]);

       const RealTensor& kappa = mod.get_total_thermal_conductivity();
       double heat_source = mod.get_total_heat_source();
      
       for (unsigned int i = 0; i < n_dofs; i++)
       {
         for (unsigned int j = 0; j < n_dofs; j++)
           Ke(i, j) += JxW[qp] * dphi[i][qp] * (kappa * dphi[j][qp]);

	   Fe(i) += JxW[qp] * heat_source * phi[i][qp];
       }

     }

     // the sides
     for (unsigned int s = 0; s < elem->n_sides(); s++)
     {
       const ElementSide elside(elem->top_parent(),s);
       if (se.is_outer_boundary(elside) ||
	   is_on_GF_boundary(elside))
       {
	 
	 BoltzmannBoundaryModel* mod =
	   get_interface_model<BoltzmannBoundaryModel>(elem, s);
	 
	 bool do_boundary = false;
	 double a, b, c;
	 if (mod != NULL)
	 {
	   mod->calculate(elem,s,elem->centroid() );
	   mod->get_coefficients(a, b, c);
	   do_boundary = true;  
	 }
	 else
	 {             
	   if (is_on_GF_boundary(elside) && !first_guess) //If this an Gray/Fourier Boundary
	   {
 	      ID dof = (elem->neighbor(s))->dof_number(gray_sys_number,0,0);
 	      do_boundary = true;

 	      if (is_gray_solved)  //-----------Get the Gray flux-------------------------------
 	      {
		RealGradient heat_flux(0);
		for (ID i = 0; i < 3; i++)
		  heat_flux(i) = (*thermal_flux[i])(dof);

		//Put here something
	        double normal_flux =  heat_flux * normal[0];
		
	        a = 0;
		b = 1;
		c = -normal_flux;
	      }
	      else //Impose the computed temperature 
	      {
	        double T = (*equilibrium_energy)(dof);
		
	        a = 1;
	      	b = 0;
	        c = T;
		
	      }
	   }
	 }//if (is_on_GF_boundary(elside) && !first_guess
	 
	 if (do_boundary)
	 {
	   fe_face->reinit(elem, s);
	   if ((b < 1e-10) && (b >= 0)) b = 1e-20;
	   else if ((b > -1e-10) && (b<= 0)) b = -1e-20;
	   a /= b;
	   c /= b;
	   
	   for (unsigned int qp = 0; qp < qface.n_points(); qp++)
	     for (unsigned int i = 0; i < n_dofs; i++)
	     {
	       for (unsigned int j = 0; j < n_dofs; j++)
		 Ke(i, j) += a * JxW_face[qp] * (phi_face[i][qp] * phi_face[j][qp]);
	       
	       Fe(i) += c * JxW_face[qp] * phi_face[i][qp];
	     } 
	 } 
       }
       
     }
     

     dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
     system_fourier.matrix->add_matrix(Ke, dof_indices);
     system_fourier.rhs->add_vector(Fe, dof_indices);

   }//Elem
   //cout<<n_elem<<endl;
   system_fourier.matrix->close();
   system_fourier.matrix->print_matlab("K.m");
   system_fourier.rhs->close();
   system_fourier.rhs->print_matlab("F.m");
}

bool
Boltzmann::get_fourier_boundary(ElementSide elside,double a, double b, double c)
{

  bool boundary = false;

  BoltzmannBoundaryModel* mod =
    get_interface_model<BoltzmannBoundaryModel>(elside.elem(), elside.side());

  if (mod != NULL)
  {
    mod->calculate(elside.elem(), elside.side(),elside.elem()->centroid() );
    mod->get_coefficients(a, b, c);
    boundary = true;
  }
  else
  {
    // if (domain_boundary.count(elside.elem())) //If this an Gray/Fourier Boundary
    // {


    // }

  }


}


void
Boltzmann::AngularIntegrator::compute_custom_direction(std::vector<Point> custom_dir)
{

  //Omega is not theta dependent but simply uniform.
  weight = 1.0;
  //total_angle = 4.0 * M_PI;


  n_slices = custom_dir.size();

  total_angle = n_slices;

  directions.resize(n_slices);
  d_omega.resize(n_slices,1);
  dir.resize(n_slices);
  theta_vec.resize(n_slices,0);
  phi_vec.resize(n_slices,0);
  spec.resize(n_slices);
  dir = custom_dir;
  directions = custom_dir;

  // directions = custom_dir * 2*M_PI;

  //for (ID k = 0; k<n_slices; k++)
    // {

    //directions[k] = custom_dir[k];// * M_PI;
    //d_omega[k] = 4.0 * M_PI/n_slices;
  // d_omega[k] = 1;
    //directions[k] = custom_dir[k] * d_omega[k];
  //  theta_vec[k] = 0.0;
  //  phi_vec[k] = 0.0;

  // }

  //Spec vectors
  //for (ID k1 = 0;k1<n_slices; k1++)
  //  for (ID k2 = 0;k2<n_slices; k2++)

  //  {
  //   Point sum = directions[k1] + directions[k2];
  //   if (sum.size() < 1e-4)
  //spec[k1]=k2;
  // }


}
void
Boltzmann::AngularIntegrator::compute_directions()
{


  //double min_theta, max_theta, min_phi, max_phi;
      //double d_phi,theta, phi, d_theta;
    /*ID k = 0;

    switch (dim)
    {

    case 1 :

     theta_slices = 1;
      if (phi_slices == 0)
           phi_slices = 4;

      min_theta = 0.0;
      max_theta = M_PI;

      min_phi = M_PI * 0.5;
      max_phi = M_PI * 2.0 + M_PI * 0.5;


      weight = 1.0;
      n_slices = theta_slices * phi_slices;
      spec.resize(n_slices);
      total_angle = 2 * M_PI;

      directions.resize(n_slices);
      integrate_directions.resize(n_slices);
      d_omega.resize(n_slices);
      dir.resize(n_slices);
      theta_vec.resize(n_slices);
      phi_vec.resize(n_slices);

      d_phi =  (max_phi - min_phi) / phi_slices;

      k = 0;
      for (ID n_phi = 0; n_phi < phi_slices; n_phi++)
      {

       phi = min_phi + d_phi * n_phi + phi_zero * M_PI/180.0;

       d_omega[k] = d_phi;

       dir[k](0) = sin(phi);
       dir[k](1) = cos(phi);
       dir[k](2) = 0.0;

       directions[k](0) =  2.0 * sin(phi) * sin(0.5 * d_phi);
       directions[k](1) =  2.0 * cos(phi) * sin(0.5 * d_phi);
       directions[k](2) =  0.0;

       theta_vec[k] = M_PI/2.0;
       phi_vec[k] = phi;

       k++;

      }

      break;

    case 2 :

      //theta_slices = 1;
      //if (phi_slices == 0)

      min_theta = 0.0;
      max_theta = M_PI * 0.5;

      min_phi = M_PI * 0.5;
      max_phi = M_PI * 2.0 + M_PI * 0.5;

      weight =  2.0;

      n_slices = theta_slices * phi_slices;
      spec.resize(n_slices);

      directions.resize(n_slices);
      integrate_directions.resize(n_slices);
      d_omega.resize(n_slices);
      dir.resize(n_slices);
      theta_vec.resize(n_slices);
      phi_vec.resize(n_slices);

      d_theta =  (max_theta - min_theta) / theta_slices;
      d_phi =  (max_phi - min_phi) / phi_slices;

      k = 0;
      for (ID n_phi = 0; n_phi < phi_slices; n_phi++)
      {
        // phi = min_phi + d_phi * 0.5 + d_phi * n_phi;

        phi = min_phi + d_phi * n_phi + phi_zero * M_PI/180.0 ;

        for (ID n_theta = 0; n_theta < theta_slices; n_theta++)
        {
          theta = min_theta + d_theta * 0.5 + d_theta * n_theta;

          d_omega[k] =  weight * 2.0 * sin(theta) * sin (0.5 * d_theta) * d_phi;
          //d_omega[k] = d_phi;

          dir[k](0) = sin(theta) * sin(phi);
          dir[k](1) = sin(theta) * cos(phi);
          dir[k](2) = cos(theta);

          directions[k](0) =  weight * sin(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
          directions[k](1) =  weight * cos(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
          directions[k](2) =  weight * 0.5 * d_phi * sin(2.0 * theta) * sin(d_theta);

          theta_vec[k] = theta;
          phi_vec[k] = phi;

          k++;
        }
      }


      //Specular vectors
    //  for (ID k1 = 0; k1<n_slices; k1++)
      //    for (ID k2 = 0; k2<n_slices; k2++)
        //    {
          //    Point sum = directions[k1] + directions[k2];
           //   cout<<sum<<endl;
            //  if ((sum(0) < 1e-4) && (sum(1) < 1e-4) && (sum(2) < 1e-4) )
            //     spec[k1]=k2;
           // }


      total_angle = 4.0 * M_PI;
      break;

    case 3 :

      std::cout<<"3D"<<std::endl;

      if (theta_slices == 0)
	theta_slices = 2;

      if (phi_slices == 0)
	phi_slices = 4;


      min_theta = 0.0;
      max_theta = M_PI;

      min_phi = 0.0;
      max_phi = M_PI * 2.0;

      weight = 1.0;


      n_slices = theta_slices * phi_slices;

      spec.resize(n_slices);

      total_angle = 4.0 * M_PI;



      directions.resize(n_slices);
      integrate_directions.resize(n_slices);
      d_omega.resize(n_slices);
      dir.resize(n_slices);
      theta_vec.resize(n_slices);
      phi_vec.resize(n_slices);

      d_theta =  (max_theta - min_theta) / theta_slices;
      d_phi =  (max_phi - min_phi) / phi_slices;

      k = 0;
      for (ID n_phi = 0; n_phi < phi_slices; n_phi++)
      {
        // phi = min_phi + d_phi * 0.5 + d_phi * n_phi;
        phi = min_phi + d_phi * n_phi + phi_zero * M_PI/180.0 ;

        for (ID n_theta = 0; n_theta < theta_slices; n_theta++)
        {
          theta = min_theta + d_theta * 0.5 + d_theta * n_theta;

          d_omega[k] =  2.0 * sin(theta) * sin (0.5 * d_theta) * d_phi;
          //d_omega[k] = d_phi;

          dir[k](0) = sin(theta) * sin(phi);
          dir[k](1) = sin(theta) * cos(phi);
          dir[k](2) = cos(theta);

          directions[k](0) =  weight * sin(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
          directions[k](1) =  weight * cos(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
          directions[k](2) =  weight * 0.5 * d_phi * sin(2.0 * theta) * sin(d_theta);

          theta_vec[k] = theta;
          phi_vec[k] = phi;

          k++;
        }
      }

      break;
    }*/

  double delta = 0.0;

  double  min_theta = 0.0 + delta;
  double  max_theta = M_PI+ delta;

  double  min_phi = M_PI * 0.5 + delta;
  double  max_phi = M_PI * 2.0 + M_PI * 0.5+ delta;

    weight =  1.0;

    n_slices = theta_slices * phi_slices;
    //spec.resize(n_slices);

     directions.resize(n_slices);
     integrate_directions.resize(n_slices);
     d_omega.resize(n_slices);
     dir.resize(n_slices);
     theta_vec.resize(n_slices);
     phi_vec.resize(n_slices);

     double d_theta =  (max_theta - min_theta) / theta_slices;
     double d_phi =  (max_phi - min_phi) / phi_slices;

     ID k = 0;
     for (ID n_phi = 0; n_phi < phi_slices; n_phi++)
     {
       double phi = min_phi + d_phi * n_phi + phi_zero * M_PI/180.0 ;

       for (ID n_theta = 0; n_theta < theta_slices; n_theta++)
       {
         double theta = min_theta + d_theta * 0.5 + d_theta * n_theta;

         d_omega[k] =  weight * 2.0 * sin(theta) * sin (0.5 * d_theta) * d_phi;

         dir[k](0) = sin(theta) * sin(phi);
         dir[k](1) = sin(theta) * cos(phi);
         dir[k](2) = cos(theta);

         directions[k](0) =  weight * sin(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
         directions[k](1) =  weight * cos(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
         directions[k](2) =  weight * 0.5 * d_phi * sin(2.0 * theta) * sin(d_theta);

         theta_vec[k] = theta;
         phi_vec[k] = phi;

         k++;
       }
     }

    total_angle = 4.0 * M_PI;
}



void
Boltzmann::AngularIntegrator::print_info(void)
{


  std::cout<<"General:  "<<"theta slices: "<<theta_slices<<" phi slices: "<<phi_slices<<std::endl;
  //std::cout<<"d_theta:  "<<d_theta / ( M_PI) * 180.0<<"  d_phi:  "<<d_phi / ( M_PI) * 180.0<<std::endl;
  //std::cout<<"Weight factor:  "<<weight<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"  "<<std::endl;


  for (ID k =0; k<n_slices; k ++)
  {

    std::cout<<"Direction: "<<k+1<<std::endl;
    std::cout<<"  Theta: "<<theta_vec[k] / ( M_PI) * 180.0 <<"  phi: "<<phi_vec[k] / ( M_PI) * 180.0<<std::endl;

    std::cout<<"d_omega:  "<<  d_omega[k]<<std::endl;


    std::cout<<"sx:  "<<   dir[k](0)  <<"  sy:  "<<   dir[k](1) <<"  sz:  "<<   dir[k](2)<<std::endl;
    std::cout<<"six:  "<<   directions[k](0)  <<"  siy:  "<<   directions[k](1) <<"  siz:  "<<   directions[k](2)<<std::endl;
    std::cout<<"  "<<std::endl;

  }

  std::cout<<"  "<<std::endl;

  //check 1
  double solid_angle = 0.0;
  double s0 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  for (ID n = 0; n < n_slices; n++)
  {
    solid_angle += d_omega[n] ;
    s0 += directions[n](0);
    s1 += directions[n](1);
    s2 += directions[n](2);

  }
  std::cout<<"Check:  "<<std::endl;
  std::cout<<"Solid angle: "<<solid_angle/total_angle<<std::endl;


  std::cout<<"Directions:  "<<std::endl;
  std::cout<<s0<<std::endl;
  std::cout<<s1<<std::endl;
  std::cout<<s2<<std::endl;



  //std::cout<<"  "<<std::endl;
  //std::cout<<"TotalAngle:  "<<total_angle<<std::endl;
  //std::cout<<"  "<<std::endl;
  //std::cout<<"Specular vector: "<<std::endl;
  //for (ID k =0; k<n_slices;k++)
  //  std::cout<<"Dir "<<k<<" : "<<spec[k]<<std::endl;
  //std::cout<<"  "<<std::endl;
}


void
Boltzmann::AngularIntegrator::print_info(ID k)
{



    std::cout<<"Direction: "<<k<<std::endl;
    std::cout<<"  Theta: "<<theta_vec[k] / ( M_PI) * 180.0 <<"  phi: "<<phi_vec[k] / ( M_PI) * 180.0<<std::endl;

    std::cout<<"d_omega:  "<<  d_omega[k]/(4.0 * M_PI) <<std::endl;


    std::cout<<"sx:  "<<   dir[k](0)  <<"  sy:  "<<   dir[k](1) <<"  sz:  "<<   dir[k](2)<<std::endl;
    std::cout<<"six:  "<<   directions[k](0)  <<"  siy:  "<<   directions[k](1) <<"  siz:  "<<   directions[k](2)<<std::endl;
    std::cout<<"  "<<std::endl;



  std::cout<<"  "<<std::endl;

}


  //cout<<Fe(0)<<endl;

    // END ASSEMBLY1


     //ASEEMBLY FINAL-------------

  //    ID OB = 0;
//      double f_value = 0.0;
//      for (ID ns = 0; ns < elem->n_sides(); ns++)
//      {
//        fe_face->reinit(elem,ns);
//        double in = dir * normal[0];
//        const ElementSide elside(elem->top_parent(),ns);

//        BoltzmannBoundaryModel* b_mod =
// 	 get_surface_model<BoltzmannBoundaryModel>(elem,ns);

//        if (in<0.0)
// 	 if (se.is_outer_boundary(elside) && (b_mod == NULL))
// 	 {
// 	   f_value = get_boundary_value(elside);
// 	   OB = 1;
// 	 }
//      }

//      if (OB==1)
//      {
//        Ke(0,0) = 1.0;
//        Fe(0) = f_value;

//      }
//      else
//      {
//        Ke(0,0) = 1.0/tg * JxW[0];

//        ID nb = 0;
//        for (ID ns = 0; ns < elem->n_sides(); ns++)
//        {

// 	 fe_face->reinit(elem,ns);
// 	 BoltzmannBoundaryModel* b_mod =
//     	   get_surface_model<BoltzmannBoundaryModel>(elem,ns);
// 	 const ElementSide elside(elem->top_parent(),ns);
// 	 double in = dir * normal[0];
// 	 double value = vg * (IntDir * normal[0]) * JxW_face[0];

// 	 if (in<0.0)
// 	 {
// 	   // if (se.is_outer_boundary(elside) || (is_on_GF_boundary(elside) && is_fourier_solved))
// 	   if (is_on_GF_boundary(elside) && is_fourier_solved || (b_mod != NULL))
// 	   {
// 	     //   nb ++;
// 	     Fe(0) -= get_boundary_value(elside) * value;
// 	   }
// 	   else
// 	     Ke(0,ns + 1) += value;


// 	 }
// 	 else
// 	 Ke(0,0) += value;

//        }
//        // if (nb>0)
//        // Fe(0) /=nb;

//        Fe(0)   += (e_0/tg  +  heat_source * d_omega) * JxW[0];

//      }

     //-------------------------------



//       //Assembly1

//      //Check if there are outer boundary on in negative direction

//       ID OB = 0;
//      double f_value = 0.0;
//      for (ID ns = 0; ns < elem->n_sides(); ns++)
//      {
//        fe_face->reinit(elem,ns);
//        double in = dir * normal[0];
//        const ElementSide elside(elem->top_parent(),ns);

//        BoltzmannBoundaryModel* b_mod =
// 	 get_surface_model<BoltzmannBoundaryModel>(elem,ns);

//        if (in<0.0)
// 	 if (se.is_outer_boundary(elside) && (b_mod == NULL))
// 	 {
// 	   f_value = get_boundary_value(elside);
// 	   OB = 1;
// 	 }
//      }
//      //----------------------

//       if (OB == 1)
//       {
//       Ke(0,0) = 1;
//       Fe(0) = f_value;
//       }
//      else
//      {

//        Ke(0,0) = 1.0/tg * JxW[0];
//        ID nb = 0;
//        for (ID ns = 0; ns < elem->n_sides(); ns++)
//        {

// 	 fe_face->reinit(elem,ns);
// 	 BoltzmannBoundaryModel* b_mod =
// 	   get_surface_model<BoltzmannBoundaryModel>(elem,ns);

// 	 const ElementSide elside(elem->top_parent(),ns);
// 	 double in = dir * normal[0];
// 	 double value = (IntDir * normal[0]) * vg * tg * JxW_face[0]/d_omega;

// 	 if (in<0.0)
// 	 {
// 	   //if (is_on_GF_boundary(elside) && is_fourier_solved || (b_mod != NULL))
// 	   if (se.is_outer_boundary(elside) || (is_on_GF_boundary(elside) && is_fourier_solved))
// 	     Fe(0) -= get_boundary_value(elside) * value;
// 	   else
// 	     Ke(0,ns + 1) += value;

// 	   // nb ++;
// 	 }
// 	 else
// 	   Ke(0,0) += value;

// 	 if (nb>0)
// 	   Fe(0) /=nb;

// 	 Fe(0)  += (e_0  +  heat_source * tg / cg) * JxW[0];

//        }
//         }

    // END ASSEMBLY1
