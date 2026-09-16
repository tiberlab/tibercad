/*  
 * This file is part of the tiberCAD module degradation.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file Degradation.cpp
 * \brief tiberCAD degradation module implementation.
 *
 * \note This file is part of module degradation.
 */


#include "Degradation.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/solver/TiberTransientSystem.h"
#include "tibercad/io/Messages.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "tibercad/module/TiberModule.h" // This is needed in order to create the shared module library
#include <cmath>


using namespace std;
using namespace libMesh;


Degradation::Degradation(const ModelOptions& options) :
  SimulationInterface(options),
  _my_assembly(this)
{
  // there's nothing to be done
}


Degradation::~Degradation(void)
{
  // there's nothing to be done
}


void
Degradation::do_init(void)
{
  parse_options();

  std::string water_ingress_sim = get_option("p_h2o_sim", "wi"); // get p_h2o from WI simulation 
  _humidity_model = SimulationInterface::find_solution_provider(water_ingress_sim, "PartialPressure");

  std::string oxygen_ingress_sim = get_option("p_o2_sim", "oi"); // get p_o2 from OI simulation 
  _oxygen_model = SimulationInterface::find_solution_provider(oxygen_ingress_sim, "PartialPressure");

  // create a transient linear equation system 
  create_equation_system("transient");

  // get the reference to it
  TiberTransientLinSystem& system = get_equation_system<TiberTransientLinSystem>();

  // we use the fraction of active perovskite as variable "f"
  // all other quantities are then related to "f"
  system.add_variable("f", CONSTANT, MONOMIAL, &(this->get_region_ids())); 
  system.attach_assemble_object(_my_assembly);
  system.init();
  system.get_local_solution_vector().zero();  
  system.get_local_solution_vector().add(1.0);  
  system.update();
}


void
Degradation::parse_options(void)
{
  // temperature, and photon flux ( -> generated carriers) from input file
  _cell_temp = get_option("cell_temperature", SimulationOptions::temperature);  // [K]
  _el_den = get_option("n", 1.5e21); // [(ph. m^-2 s^-1)^0.72] // 1 Sun equivalent photon flux, ^0.72 when computing reaction rate
}


double
Degradation::do_eval_deg_rate(const libMesh::Elem* elem,
    const libMesh::Point& p)
{
  //physical constants
  const double k_b = 8.61733e-5; //[eV/K]

  //list of parameters from J Mater Chem A 2025 13 38436
  //DPO
  const double k_0_DPO = 1.5e-16; //[mol m^-2 s^-1 kPa^-1 (ph. m^-2 s^-1)^-0.72]
  const double E_A_DPO = 0.56; //[eV]
  const double K_2 = 9.5e-4; //[kPa^-1]
  const double K_3 = 2.3e-14; //[(ph. m^-2 s^-1)^-0.72]

  //WPO
  const double k_0_WPO = 1.3e-27; //[mol m^-2 s^-1 kPa^-1 (ph. m^-2 s^-1)^-0.72]
  const double E_A_WPO = -0.12; //+- 0.12 [eV]

  //hydration
  const double f_hyd = 0.55; //[unitless]
  const double k_0_hyd = 3.2e22; //[kPa^-1]
  const double E_A_hyd = 1.5; //[eV]

  //thermally activated processes
  double K_DPO = k_0_DPO*exp(-E_A_DPO/(k_b*_cell_temp)); //[mol m^-2 s^-1 kPa^-1 (ph. m^-2 s^-1)^-0.72]
  double K_WPO = k_0_WPO*exp(-E_A_WPO/(k_b*_cell_temp)); //[mol m^-2 s^-1 kPa^-1 (ph. m^-2 s^-1)^-0.72]
  double K_hyd = k_0_hyd*exp(-E_A_hyd/(k_b*_cell_temp)); //[kPa^-1]
  
  //expression of degradation rate as a function of humidity (WPO and hydration factor) from water ingress simulation
  //and oxygen pressure from oxygen ingress simulation
  double rate = 0;
  double den = 0; //mixed order denominator
  double p_o2 = 0;
  double r_DPO = 0; //DPO
  double humidity = 0;
  double r_WPO = 0; //WPO
  
  if (_humidity_model.first != nullptr)
  {
    if (_humidity_model.first->get_solution(elem, _humidity_model.second, humidity, p, false))       
    {
      if (_oxygen_model.first != nullptr)
      {
        if (_oxygen_model.first->get_solution(elem, _oxygen_model.second, p_o2, p, false))
        {
         den = 1 + K_2*p_o2*1e-3*(1 + K_3*pow(_el_den,0.72)); //[unitless]
         r_DPO = K_DPO*p_o2*1e-3*pow(_el_den,0.72)/den; //[mol m^-2 s^-1]
         r_WPO = K_WPO*p_o2*1e-3*pow(_el_den,0.72)*humidity*1e-3/pow(den,2); //[mol m^-2 s^-1]
         rate = (r_DPO + r_WPO)*(1 - f_hyd*K_hyd*humidity*1e-3/(den + K_hyd*humidity*1e-3)); //[mol m^-2 s^-1]
        }
      }
    }
  }

  return rate; //[mol/m^2/s]

}


void
Degradation::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(PerovskiteFraction, REAL, CELL, ""); //unitless
}


void
Degradation::do_solve(void)
{

  //we need to sweep over time
  double current_time = TiberCad::get_global_time(); //[s]
  //same as WI simulation
  TiberTransientLinSystem& system = get_equation_system<TiberTransientLinSystem>();
  system.set_options(get_solver_options());
  system.set_target_time(current_time);
  system.solve();

}


void
Degradation::do_print_info(void)
{
  Messages::info("Kinetic model for the perovskite degradation");
}


void
Degradation::assemble(void)
{

  const double rho = 4201; // [kg/m^3] perovskite density
  const double w = 0.627; // [kg/mol] perovskite molar weight
  const double t_psk = 280e-9; // [m] perovskite layer thickness
  const double c = rho*t_psk/w; // [mol/m^2] perovskite concentration

  // Here, (1+delta_t*r')f_n+1 = f_n, where
  // r' = r/c, delta_t = time step, f = fraction of active perovskite

  // In WI, (1+delta_t*A)f_n+1 = f_n, where
  // (1+delta_t*A) = K -> f_n+1 = K^-1 f_n

  TiberTransientLinSystem& system = get_equation_system<TiberTransientLinSystem>();

  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();
  
  // this vector contains the diagonal of the A matrix
  NumericVector<libMesh::Number>& t_weight = system.get_vector("t_weight");
  t_weight.zero(); //see transient system def

  // the time relevant for assembly
  double current_time = system.time;

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  DofMap& dof_map = system.get_dof_map();

  const unsigned int uvar = system.variable_number("f");

  FEType fe_type = dof_map.variable_type(uvar);

  // the finite element
  unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));
  unique_ptr<QBase> qrule(QBase::build(QGAUSS, dim, CONSTANT));
  fe->attach_quadrature_rule(qrule.get());
  
  // as we have a CONSTANT MONOMIAL, fe->get_JxW()[0] is the cell area (if dim = 2)
  const vector<Real>& JxW = fe->get_JxW(); 

  // as we have a CONSTANT MONOMIAL, fe->get_phi()[0][0] is always 1.0
  const vector<vector<Real> >& phi = fe->get_phi(); 

  // add this if you need the physical coordinates of the quadrature points
  const vector<Point>& q_point = fe->get_xyz(); 

  vector<unsigned int> dof_indices;

  DenseMatrix<Number> Ke;
  DenseVector<Number> A;
  DenseVector<Number> sol;

  MeshBase::const_element_iterator el(this->active_local_elements_begin());
  const MeshBase::const_element_iterator end_el(this->active_local_elements_end());

  for ( ; el != end_el ; ++el) // loop over elements
  {
    const Elem* elem = *el; 

    dof_map.dof_indices(elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();

    fe->reinit(elem);

    // resize the element matrix/rhs (does also zero them out)
    Ke.resize(n_dofs, n_dofs);
    A.resize(n_dofs);
    sol.resize(n_dofs);

    // get the current solution
    dof_map.extract_local_vector(solution, dof_indices, sol); 

    // we get P_H2O and P_O2 and evaluate the degradation rate for each centroid
    double rate = do_eval_deg_rate(elem,elem->vertex_average()); 

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {
      for (unsigned int i = 0; i < n_dofs; i++)
      {
        A(i) += 1.0;
        Ke(i,i) += rate/c;
      }
    }

    system.matrix->add_matrix(Ke, dof_indices);
    t_weight.add_vector(A, dof_indices);

  }

  system.matrix->close();
  t_weight.close();

}

void
Degradation::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{

  TiberTransientLinSystem& system = get_equation_system<TiberTransientLinSystem>();
  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();
  const DofMap& dof_map = system.get_dof_map();
  const unsigned int u_var = system.variable_number("f");
  vector<unsigned int> dof_indices;
  dof_map.dof_indices(elem, dof_indices, u_var);

  // cell data variable
  double f = solution(dof_indices[0]);
  // we assure f >= 0
  if (f < 0) f = 0;
  values[PerovskiteFraction][0] = f;
  
}







