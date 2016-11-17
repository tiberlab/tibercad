// $Id: ElectricalContact.C 4135 2015-09-25 10:19:38Z maufder $

#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"
#include "ParticleDensity.h"
#include "Initializer.h"


ElectricalContact::ElectricalContact(const ModelOptions& options)
  : DDInterfaceModel(options),
    _voltage(0.0),
    _surfres(0.0),
    _contact_fermilevel(0.0)
{
  // this is a real contact
  has_current(true);
}


void
ElectricalContact::do_init(void)
{
  DDInterfaceModel::do_init();
  double vrec = -1.0;
  _vrec.resize(n_carriers(), 0);
  _fixed_vrec.resize(n_carriers(), 0);

  if (get_option("zero_field", false))
    set_type(0, NEUMANN);

  get_parameter("rec_velocity", vrec);

  for (unsigned int i=0; i < n_carriers(); i++)
  {
    _vrec[i] = vrec;

    if (get_option("zero_grad_fermi", false) || (_vrec[i] >= 0))
    {
      set_type(i+1, NEUMANN);
      _fixed_vrec[i] = true;
    }
    else
      _vrec[i] = 0.0;
  }

  get_parameter("contact_resistance", _surfres);

  get_parameter("voltage", _voltage);

}


void
ElectricalContact::do_compute(void)
{
  
  if (get_type(0) != NEUMANN)
    coeff_g(0) = _contact_fermilevel + get_inner_voltage();

  for (auto&& it : qids())
  {
    if (get_type(it.second) != NEUMANN)
      coeff_g(it.second) = get_inner_voltage();
    else
    {
      const PointData& pd = get_point_data();

      // we have to take the correct equilibrium density!
      CarrierProperties* cp = get_bulk_dd_properties()->get_carrier_properties(it.first);
      //cp->calculate(pd.carrier_vt[it.first]);
      cp->calculate(pd.carrier_vt.at(it.first));
      std::pair<double, double> dens(cp->get_density_and_derivative(0, _contact_fermilevel));
      double dens0 = dens.first; 

      //double R = _vrec[it.second-1] * (pd.q_density[it.first] - dens0);
      //double dR = _vrec[it.second-1] * pd.q_density_derivative[it.first];
      double R = _vrec[it.second-1] * (pd.q_density.at(it.first) - dens0);
      double dR = _vrec[it.second-1] * pd.q_density_derivative.at(it.first);

      coeff_g(it.second) += R;
      jacobian(it.second, 0) += dR;
      jacobian(it.second, it.second) -= dR;
    }
  }
/*
  if (get_type(2) != NEUMANN)
    coeff_g(2) = get_inner_voltage();
  else
  {
    const PointData& pd = get_point_data();

    // we have to take the correct equilibrium density!
    CarrierProperties& vb = get_bulk_dd_properties()->get_valence_band();
    vb.calculate(pd.hole_vt);
    std::pair<double, double> hldens(vb.get_density_and_derivative(0, _contact_fermilevel));
    double p0 = hldens.first;
    //ParticleDensity& hl = get_dd_properties()->get_holes();
    //vb.set_temperature(pd.hole_vt);
    //hl.set_classical_parameters(vb.get_effective_DOS(),
    //    _contact_fermilevel - get_dd_properties()->get_valence_band_edge(), 0,
    //    pd.hole_vt);
    //double p0 = hl.get_particle_density();

    double Rp = _vrec_p * (pd.hole_density - p0);
    double dRp = _vrec_p * pd.hole_density_derivative;

    coeff_g(2) += Rp;
    jacobian(2, 0) += dRp;
    jacobian(2, 2) -= dRp;
  }
 */
}

/*
double
ElectricalContact::get_contact_voltage_drop(void) const
{
  double j = Constants::e * (get_normal_hole_flux() -
      get_normal_electron_flux());

  // a negative current means inflowing current
  return -_surfres * j;
}
*/



