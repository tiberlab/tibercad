// $Id: ThermalBalance.C 2457 2011-03-06 23:52:12Z gromano $

#include "Thermal.h"
#include "ThermalModel.h"
#include "ThermalBoundaryModel.h"
#include "TiberLinearSystem.h"
#include "Messages.h"
#include "SimulationOptions.h"
#include "ModelOptions.h"

#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/fe_interface.h"

#include "TiberModule.h"


using namespace std;


Thermal*
Thermal::_this = NULL;



Thermal::Thermal(const ModelOptions& options) :
  SimulationInterface(options)
{
  // there's nothing to be done
}



Thermal*
Thermal::create(const ModelOptions& options)
{
  // we could use the options to create different implementations
  // or something like that.
  return new Thermal(options);
}



void
Thermal::do_init(void)
{
  parse_options();

  ID dim = get_mesh().mesh_dimension();

  create_equation_system("linear");
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();


  system.add_variable("T", libMesh::FIRST, &(this->get_region_ids()));
  system.attach_assemble_function(assemble);
  system.init();

  libMesh::NumericVector<double>& solution = system.get_local_solution_vector();
  solution.zero();
  solution.add(SimulationOptions::T);
  solution.close();
  system.update();

}





//-------------------------------------------------------------------------//
Thermal::~Thermal()
{
   

}


void
Thermal::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(LatticeTemp, REAL, NODES, "K");
  declare_solution(MaxTemp, REAL, GLOBAL, "K");
  declare_solution(ThermalFlux, VECTOR, NODES, "W/m^2");
  declare_solution(ThermCond, VECTOR, NODES, "W/m K");
  declare_solution(HeatSource, REAL, NODES, "W/m^3");

}



double
Thermal::compute_power_dissipated()
{

  //Gray System
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const libMesh::NumericVector<double>& solution = *(system.solution);

  libMesh::DofMap& dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;
  //-----------------------------------------------


  const unsigned int tvar = system.variable_number("T");
  libMesh::FEType fe_type = dof_map.variable_type(tvar);

  
  ID  dim = get_mesh().mesh_dimension();

  libMesh::UniquePtr<libMesh::FEBase> fe_face(build_finite_element(dim, fe_type, true));
  libMesh::QGauss qrule_face(dim-1,libMesh::FIFTH);
  fe_face->attach_quadrature_rule(&qrule_face);
  
  const std::vector<Point>& q_point_face = fe_face->get_xyz();
  const std::vector<std::vector<libMesh::RealGradient> >&  dphi = fe_face->get_dphi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& normal = fe_face->get_normals();
  
  double power_dissipated = 0.0;
 

  MeshBase::const_element_iterator       it  = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end = this->active_local_elements_end();

  for ( ; it != end; ++it)
  {

    const Elem* elem = *it;
    dof_map.dof_indices(elem, dof_indices);
    
    ThermalModel& mod = *get_bulk_model<ThermalModel>(elem);
   
    const libMesh::RealTensor& kappa = mod.get_total_thermal_conductivity();
 
    for (ID ns = 0; ns<elem->n_sides(); ns++)
    {
      const ElementSide elside(elem->top_parent(),ns);
      if ( get_environment().is_on_boundary(elside))
      {

	fe_face->reinit(elem,ns);
	for (ID qp = 0; qp <  qrule_face.n_points(); qp++)
	{ 
	
	  for (ID alpha = 0; alpha<dof_indices.size() ;alpha ++)
	  {
	    power_dissipated -= JxW_face[qp] * solution(dof_indices[alpha]) * ((kappa * dphi[alpha][qp]) * normal[qp]);  

	  }
	}
       
      }
    }
    
  }

  this->get_solver_communicator().sum(power_dissipated);

  return(power_dissipated);
}

double
Thermal::compute_power_emitted()
{

  //Gray System
  libMesh::EquationSystems& es = get_equation_systems();
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();
  const libMesh::NumericVector<Number>& solution = system.get_solution_vector();
  libMesh::DofMap& dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;
  //-----------------------------------------------
  ID  dim = get_mesh().mesh_dimension();

  const unsigned int tvar = system.variable_number("T");
  libMesh::FEType fe_type = dof_map.variable_type(tvar);

  //------------BULK----------
  libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type, true));
  libMesh::QGauss qrule(dim,libMesh::FIFTH);
  fe->attach_quadrature_rule(&qrule);
  const vector<vector<Real> >& phi = fe->get_phi();
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  //--------------------------


  libMesh::UniquePtr<libMesh::FEBase> fe_face(build_finite_element(dim, fe_type, true));
  libMesh::QGauss qrule_face(dim-1,libMesh::CONSTANT);
  fe_face->attach_quadrature_rule(&qrule_face);

  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<Point>& normal = fe_face->get_normals();

  MeshBase::const_element_iterator       it  = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end = this->active_local_elements_end();

  Real total_heat_source = 0.0;
  for ( ; it != end; ++it)
  {

    const Elem* elem = *it;
    dof_map.dof_indices(elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();

    fe->reinit(elem);

    ThermalModel& mod = *get_bulk_model<ThermalModel>(elem);

    //Energy emitted
    for (ID qp = 0; qp <  qrule.n_points(); qp++)
    {
      double T  = 0.0;
      for (unsigned int i = 0; i < n_dofs; i++)
        T += phi[i][qp] * solution(dof_indices[i]);

      mod.calculate(elem, q_point[qp], T);
      Real H = mod.get_total_heat_source();
      total_heat_source += H * JxW[qp];
    }
 
  }

  this->get_solver_communicator().sum(total_heat_source);

  return(total_heat_source);
}



void
Thermal::do_solve(void)
{
  _this = this;

  libMesh::EquationSystems& es = get_equation_systems();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  system.set_options(get_solver_options());
  system.solve();

  const libMesh::NumericVector<Number>& solution = system.get_solution_vector();
  _max_temperature = solution.linfty_norm();
  this->get_solver_communicator().max(_max_temperature);

  ostringstream os;
  os << "Maximum temperature : " << _max_temperature << " K";
  Messages::info(os.str());
}


void
Thermal::do_print_info(void)
{
  // Messages::info("THERMONEO");
}


PhysicalModel*
Thermal::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{

  return  ThermalModel::create(mat, options);

}



PhysicalModel*
Thermal::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{
  return ThermalBoundaryModel::create(boundary, options);
}



void
Thermal::get_solution_secure(std::map<ID, std::vector<double> >& values)
{
  if (values.count(MaxTemp))
  {
    values[MaxTemp] = std::vector<double>(1, _max_temperature);
  }
}



void
Thermal::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
   unsigned int np = p.size();

   TiberLinearSystem* system = &get_equation_system<TiberLinearSystem>();
   const libMesh::NumericVector<Number>& solution = system->get_solution_vector();
   const libMesh::DofMap& dof_map = system->get_dof_map();

   ID  dim = get_mesh().mesh_dimension();
     
   //-------------
   const unsigned int u_var = system->variable_number("T");

   libMesh::FEType fe_type = system->variable_type(u_var);
   libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type));

   vector<unsigned int> dof_indices;

   //element shape functions
   const vector<vector<Real> >& phi = fe->get_phi();
   const vector<vector<libMesh::RealGradient> >& dphi = fe->get_dphi();
   const vector<Point>& real_pts = fe->get_xyz();

   ID subdomain = elem->subdomain_id();

   fe->reinit(elem, &p);

   dof_map.dof_indices(elem, dof_indices, u_var);

   const unsigned int n_dofs = dof_indices.size();
   ThermalModel& mod = *get_bulk_model<ThermalModel>(elem);
   const libMesh::RealTensor& kappa = mod.get_total_thermal_conductivity();

   for (unsigned int n = 0; n < np; n++)
    {
      double T  = 0.0;
      for (unsigned int i = 0; i < n_dofs; i++)
        T += phi[i][n] * solution(dof_indices[i]);

      if (values.count(LatticeTemp))
	values[LatticeTemp][n] = T;
   
      mod.calculate(elem,real_pts[n], T);

      if (values.count(ThermalFlux))  
      {

        libMesh::RealGradient heat_flux(0);
	for (ID alpha = 0; alpha<dof_indices.size() ;alpha ++)
	  heat_flux -= solution(dof_indices[alpha]) * (kappa * dphi[alpha][n]);
	
	

	for (ID d = 0; d < dim; d++)
	  values[ThermalFlux][d + 3 * n] = heat_flux(d);
	
      }

      
      if (values.count(ThermCond))
      {
	const libMesh::RealTensor& kappa = mod.get_total_thermal_conductivity();
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


}


void
Thermal::do_assemble(libMesh::EquationSystems& es, const std::string& system_name)
{
  TiberLinearSystem& system_fourier = get_equation_system<TiberLinearSystem>();
   const libMesh::NumericVector<Number>& solution = system_fourier.get_solution_vector();

   const MeshBase& mesh = get_mesh();

   libMesh::DofMap& dof_map =  system_fourier.get_dof_map();

   const unsigned int tvar = system_fourier.variable_number("T");

   libMesh::FEType fe_type = dof_map.variable_type(tvar);
   ID  dim = get_mesh().mesh_dimension();
   SimulationEnvironment& se = get_environment();
   // the volume finite element
   libMesh::UniquePtr<libMesh::FEBase> fe(build_finite_element(dim, fe_type, true));
   libMesh::QGauss qrule(dim, libMesh::FIFTH);
   fe->attach_quadrature_rule(&qrule);

   const vector<Real>& JxW = fe->get_JxW();
   const vector<Point>& q_point = fe->get_xyz();
   const vector<vector<Real> >& phi = fe->get_phi();
   const vector<vector<libMesh::RealGradient> >& dphi = fe->get_dphi();


   // the surface finite element
   libMesh::UniquePtr<libMesh::FEBase> fe_face(build_finite_element(dim, fe_type, true));
   libMesh::QGauss qface(dim - 1, libMesh::SIXTH);
   fe_face->attach_quadrature_rule(&qface);

   const vector<Real>& JxW_face = fe_face->get_JxW();
   const vector<Point>& qface_point = fe_face->get_xyz();
   const vector<vector<Real> >&  phi_face = fe_face->get_phi();
   const vector<vector<libMesh::RealGradient> >& dphi_face = fe->get_dphi();
   const vector<Point>& normal = fe_face->get_normals();

   vector<unsigned int> dof_indices;

   libMesh::DenseMatrix<Number> Ke;
   libMesh::DenseVector<Number> Fe;

   MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
   for ( ; el != end_el ; ++el)
   {

     const Elem* elem = *el;

     dof_map.dof_indices(elem, dof_indices);

     const unsigned int n_dofs = dof_indices.size();

     //resize the element matrix/rhs (does also zero them out)
     Ke.resize(n_dofs, n_dofs);

     Fe.resize(n_dofs);
     fe->reinit(elem);

     ThermalModel& mod = *get_bulk_model<ThermalModel>(elem);

     // loop over the quadrature points
     for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
     {
       double T  = 0.0;
       for (unsigned int i = 0; i < n_dofs; i++)
         T += phi[i][qp] * solution(dof_indices[i]);

       mod.calculate(elem, q_point[qp], T);
       
       const libMesh::RealTensor& kappa = mod.get_total_thermal_conductivity();
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
      ThermalBoundaryModel* mod_int =
        get_interface_model<ThermalBoundaryModel>(elem, s);

      if (mod_int != NULL)
      {
        fe_face->reinit(elem, s);

        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          mod_int->calculate(elem, s, qface_point[qp]);

          double a, b, c;
          mod_int->get_coefficients(a, b, c);

          // we use a penalty approach here for its simplicity
          if ((b < 1e-10) && (b >= 0)) b = 1e-10;
          else if ((b > -1e-10) && (b<= 0)) b = -1e-10;

          a /= b;
          c /= b;

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
   
   system_fourier.matrix->close();
   //system_fourier.matrix->print_matlab("K.m");
   system_fourier.rhs->close();
   //system_fourier.rhs->print_matlab("f.m");
}



