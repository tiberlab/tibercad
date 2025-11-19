// $Id: SchottkyContact.C 1912 2010-04-16 10:13:35Z maufder $

#include "SchottkyContact.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/base/ModelErrorException.h"

#include "tibercad/module/TiberModule.h"



SchottkyContact::SchottkyContact(const ModelOptions& options)
  : ElectricalContact(options)
{

}


void
SchottkyContact::do_init(void)
{
  ElectricalContact::do_init();

  set_type(n_known_carriers(), DIRICHLET);

  _band = get_option("band", "");
  _band = get_option("carrier", _band);

  if (!has_option("metal_fermilevel") &&
      !has_option("work_function") &&
      !has_option("barrier") &&
      !has_option("barrier_height"))
    throw InitFailedException("Schottky Contact: One of 'work_function' or 'metal_fermilevel' "
        "or 'barrier' must be set");

  _fixed_barrier = true;
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

  //we set the barrier
  if (has_option("barrier") || has_option("barrier_height") || _fixed_barrier)
  {

    if (_band == "")
      throw InitFailedException("Schottky Contact: in order to fix the barrier a "
          "carrier must be provided with 'band' keyword");

    if (get_carrier_properties(_band) == nullptr)
      throw InitFailedException("Schottky Contact: carrier '" + _band + "' not found'");
  }

  const DriftDiffusionProperties* dd = get_bulk_dd_properties();
  ID band_id = dd->get_carrier_id(_band);

  if (!has_option("metal_fermilevel") && !has_option("work_function"))
  {
    double barrier = get_option("barrier_height", 0.8);
    barrier = get_option("barrier", barrier);

    double charge = get_carrier_properties(band_id)->get_charge();

    _metal_Ef = (charge < 0) ? get_carrier_band_edge(band_id) - barrier :
        get_carrier_band_edge(band_id) + barrier;
  }

  if (_fixed_barrier)
  {
    // this is a dirty trick, but assures the barrier in strained case
    _metal_Ef -=  get_carrier_band_edge(band_id);
  }



  _correction_factor = get_option("correction_factor", _correction_factor);

  _thermionic_emission = get_option("thermionic_emission", true);
  if (_thermionic_emission)
  {
    for (unsigned int i = 0; i < this->n_known_carriers(); i++)
      set_type(i, NEUMANN);

    double temp = Constants::k_B * SimulationOptions::T;
    const double fac = 0.23032943 * _correction_factor;

    for ( auto&& it : get_bulk_dd_properties()->get_carrier_properties())
    {
      double v_rec = it.second->get_thermal_velocity(temp);
      set_recombination_velocities(it.first, fac*v_rec);
    }

     // use the Scott and Malliaras formula for recombination velocities;
    _scott_malliaras = get_option("scott_malliaras", false);
    _barrier_lowering = get_option("barrier_lowering", false);
  }
  else
  {
    for (unsigned int i = 0; i < this->n_known_carriers(); i++)
     set_type(i, DIRICHLET);
  }
}


void
SchottkyContact::do_compute(void)
{
  const DriftDiffusionProperties* dd = get_bulk_dd_properties();
  ID band_id = dd->get_carrier_id(_band);

  double val = _metal_Ef;
  if (_fixed_barrier)
  {
    val +=  get_carrier_band_edge(band_id);
  }

  const PointData& pd = get_point_data();

  if (_thermionic_emission)
  {
    const double fac = 0.23032943 * _correction_factor;
    const double pi = 3.14159265358979323846;
    double epsilon = Constants::e0 / (Constants::e * 100);
    libMesh::RealGradient e_field = dd->get_electric_field();

    const libMesh::RealTensor& permittivity  = dd->get_relative_permittivity();
    const libMesh::RealGradient& eps_e_field = permittivity * e_field;

    if (e_field.norm() != 0) // && !std::isnan(e_field.norm()) )
    {
      epsilon = (Constants::e0 / (Constants::e * 100)) * eps_e_field.norm() /  e_field.norm();   //this MUST be changed, it works only for homogeneous materials
    }


    if (_scott_malliaras )
    {
      const libMesh::RealTensor& permittivity  = dd->get_relative_permittivity();
      libMesh::RealGradient eps_e_field = permittivity * e_field;

      for ( auto&& it : get_carrier_properties())
      {

        //double temp = pd.carrier_vt[it.first];
        double temp = it.second->get_temperature();
        double v0 = 16 * pi * epsilon * temp * temp * dd->get_q_mobility(it.first); 

        double rC = 1 / (4 * pi * epsilon * temp);
        double f = rC * e_field.norm() / temp;
        double psi = (1 + sqrt(f) - sqrt(1 + 2*sqrt(f))) / f;
        double vE = v0 * (1/(psi * psi) - f) / 4;
        set_recombination_velocities(it.first, vE);
      }
    }
    else
    {
      for ( auto&& it : get_carrier_properties())
      {
        //double temp = pd.carrier_vt[it.first];
        double temp = it.second->get_temperature();
        double v_rec = it.second->get_thermal_velocity(temp);
        set_recombination_velocities(it.first, fac*v_rec);
      }
    }


    if (_barrier_lowering)
    {
      double F = e_field * this->get_face_normal();
      if (dd->get_carrier_properties(band_id)->get_charge() <= 0)
      {
        if (F > 0)
        {
          double delta = sqrt(F / (4 * pi * epsilon));
          val += delta;
        }
      }
      else
      {
        if (F < 0)
        {
          double delta = sqrt(-F / (4 * pi * epsilon));
          val -= delta;
        }
      }
    }

  }

  set_contact_fermilevel(val); 
  ElectricalContact::do_compute();
}



