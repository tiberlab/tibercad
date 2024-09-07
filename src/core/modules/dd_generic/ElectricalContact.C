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
  _vrec.resize(n_known_carriers(), -1.0);
  _fixed_vrec.resize(n_known_carriers(), 0);

  if (get_option("zero_field", false))
    set_type(n_known_carriers(), NEUMANN);

  get_parameter("rec_velocity", _vrec);
  get_parameter("recombination_velocity", _vrec);
  if (_vrec.size() == 1)
    _vrec.resize(n_known_carriers(), _vrec[0]);
  if (_vrec.size() < n_known_carriers())
    throw InitFailedException("Need to provide contact recombination "
        "velocity for all carriers");

  for (unsigned int i = 0; i < n_known_carriers(); i++)
  {

    if (get_option("zero_grad_fermi", false) || (_vrec[i] >= 0))
    {
      set_type(i, NEUMANN);
      _fixed_vrec[i] = true;
      if (get_option("zero_grad_fermi", false))
        _vrec[i] = 0.0;
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
  
  if (get_type(n_known_carriers()) != NEUMANN)
    coeff_g(n_known_carriers()) = _contact_fermilevel + get_inner_voltage();

  for (auto&& it : this->get_carrier_properties())
  {
    ID carrier = it.first;
    if (get_type(carrier) != NEUMANN)
      coeff_g(carrier) = get_inner_voltage();
    else
    {
      const PointData& pd = get_point_data();

      // we have to take the correct equilibrium density!
      CarrierProperties* cp = get_bulk_dd_properties()->get_carrier_properties(carrier);
      cp->calculate(pd.carrier_vt[carrier]);
      std::pair<double, double> dens(cp->get_density_and_derivative(0, _contact_fermilevel));
      double dens0 = dens.first; 

      double R = _vrec[carrier] * (pd.q_density[carrier] - dens0);
      double dR = _vrec[carrier] * pd.q_density_derivative[carrier];

      coeff_g(carrier) += R;
      jacobian(carrier, n_known_carriers()) += dR;
      jacobian(carrier, carrier) -= dR;
    }
  }
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



