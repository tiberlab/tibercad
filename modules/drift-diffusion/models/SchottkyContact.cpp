/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file SchottkyContact.C
 * \brief tiberCAD driftdiffusion module implementation.
 *
 * \note This file is part of module driftdiffusion.
 */


#include "SchottkyContact.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/base/ModelErrorException.h"

#include "tibercad/module/TiberModule.h"



SchottkyContact::SchottkyContact(const ModelOptions& options)
  : ElectricalContact(options),
    _band('c'),
    _fixed_barrier(true),
    _thermionic_emission(true),
    _metal_Ef(0)
{
  set_type(0, DIRICHLET);
  set_type(1, DIRICHLET);
  set_type(2, DIRICHLET);
}


void
SchottkyContact::do_init(void)
{

  ElectricalContact::do_init();

  double barrier = get_option("barrier_height", 0.8);
  barrier = get_option("barrier", barrier);

  std::string band = get_options().get_option("band", "c");
  if (band == "c")
    _metal_Ef = get_conduction_band_edge() - barrier;
  else if (band == "v")
    _metal_Ef = get_valence_band_edge() + barrier;
  else
    throw ModelErrorException("SchottkyContact: undefined band specified");

  _band = band[0];

  // we override when work function is given explicitly
  if (has_option("work_function"))
  {
    // if the workfunction is set then we assume that this one should
    // be fixed and not the resulting barrier
    _fixed_barrier = false;
    _metal_Ef = -get_option("work_function", _metal_Ef);
  }

  if (has_option("metal_fermilevel"))
  {
    // sanity check
    if (has_option("work_function"))
      throw ModelErrorException("You may not specify both work_function "
          "and metal_fermilevel for Schottky contact.");

    // if the workfunction is set then we assume that this one should
    // be fixed and not the resulting barrier
    _fixed_barrier = false;
    _metal_Ef = get_option("metal_fermilevel", _metal_Ef);
  }

  _fixed_barrier = get_option("fixed_barrier", _fixed_barrier);

  if (_fixed_barrier)
  {
    // this is a dirty trick, but assures the barrier in strained case
    if (_band == 'c')
      _metal_Ef -= get_conduction_band_edge();
    else
      _metal_Ef -= get_valence_band_edge();
  }

  _thermionic_emission = get_option("thermionic_emission", true);
  if (_thermionic_emission)
  {
    double temp = Constants::k_B * SimulationOptions::T;
    double vth_n = get_bulk_dd_properties()->get_conduction_band().get_thermal_velocity(temp);
    double vth_p = get_bulk_dd_properties()->get_valence_band().get_thermal_velocity(temp);

    // the effective emission velocity is fac * v_thermal

    const double fac = 0.23032943; // std::sqrt(1.0 / (6.0 * M_PI));
    set_recombination_velocities(fac * vth_n, fac * vth_p);
    
     // use the Scott and Malliaras formula for recombination velocities;
    _scott_malliaras = get_option("scott_malliaras", false);
    _barrier_lowering = get_option("barrier_lowering", false);
   }

  //_tunneling = get_option("tunneling", true);

}


void
SchottkyContact::do_compute(void)
{
  double val = _metal_Ef;
  if (_fixed_barrier)
  {
    if (_band == 'c')
      val += get_conduction_band_edge();
    else
      val += get_valence_band_edge();
  }
 
  if (_thermionic_emission)
  {
    const DriftDiffusionProperties* dd = get_bulk_dd_properties();
    const PointData& pd = get_point_data();

    double vth_n = dd->get_conduction_band().get_thermal_velocity(pd.electron_vt);
    double vth_p = dd->get_valence_band().get_thermal_velocity(pd.hole_vt);

    const double fac = 0.23032943; // std::sqrt(1.0 / (6.0 * M_PI));
    //std::cout<<"ve_th="<<vth_n<<" vh_th="<<vth_p<<std::endl;
    set_recombination_velocities(fac * vth_n, fac * vth_p);
 
    const libMesh::RealGradient& e_field = dd->get_electric_field();
    const double pi = 3.14159265358979323846;
    const libMesh::RealTensor& permittivity  = dd->get_relative_permittivity();
    const libMesh::RealGradient& eps_e_field = permittivity * e_field;
    double epsilon = Constants::e0 / (Constants::e * 100);

    if (e_field.norm() != 0 && !std::isnan(e_field.norm()) )
    {
      epsilon = (Constants::e0 / (Constants::e * 100)) * eps_e_field.norm() /  e_field.norm();   //this MUST be changed, it works only for homogeneous materials
    }

    if (_barrier_lowering )
    {
      if (_band == 'c')
      {
        double delta = sqrt(e_field.norm() / (4 * pi * epsilon));
        val += delta;
        //std::cout<<"efield="<<e_field.size()<<" epsilon="<<epsilon<<std::endl;
        //std::cout<<sqrt(e_field.size() / (4 * pi * epsilon))<<std::endl;
      }
      else
      {
        double delta = sqrt(e_field.norm() / (4 * pi * epsilon));
        val -= delta;
      }
    }

    if (_scott_malliaras )
    {
      double e_temp = pd.electron_vt;
      double h_temp = pd.hole_vt;
      double vn0 = 16 * pi * epsilon * e_temp * e_temp * get_bulk_dd_properties()->get_electron_mobility(); 
      double vp0 = 16 * pi * epsilon * h_temp * h_temp * get_bulk_dd_properties()->get_hole_mobility();

      double e_rC = 1 / (4 * pi * epsilon * e_temp);
      double e_f = e_rC * e_field.norm() / e_temp;
      //std::cout<<e_rC / e_temp<<std::endl;
      double e_psi = (1 + sqrt(e_f) - sqrt(1 + 2*sqrt(e_f))) / e_f;
      double h_rC = 1 / (4 * pi * epsilon * h_temp);
      double h_f = h_rC * e_field.norm() / h_temp;
      double h_psi = (1 + sqrt(h_f) - sqrt(1 + 2*sqrt(h_f))) / h_f;
      double vnE = vn0 * (1/(e_psi * e_psi) - e_f) / 4;
      double vpE = vp0 * (1/(h_psi * h_psi) - h_f) / 4;
      set_recombination_velocities(vnE, vpE);
      //std::cout<<"temp="<<e_temp<<" vnE="<<vnE<<" vpE="<<vpE<<std::endl;
    }
  }
  set_contact_fermilevel(val); 
  ElectricalContact::do_compute();
}



