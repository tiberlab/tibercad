// $Id$

#include "ElectricalContact.h"



ElectricalContact::ElectricalContact(const ModelOptions& options)
  : DDInterfaceModel(options),
    _voltage(0.0),
    _surfres(0.0),
    _barrier(0.0)
{
  // this is a real contact
  has_current(true);
}


void
ElectricalContact::do_init(void)
{

  if (get_option("zero_field", false))
    set_type(0, NEUMANN);

  DDInterfaceModel::do_init();

  if (get_option("zero_grad_fermi_e", false))
    set_type(1, NEUMANN);
  if (get_option("zero_grad_fermi_h", false))
    set_type(2, NEUMANN);

  get_parameter("contact_resistance", _surfres);

  get_parameter("voltage", _voltage);
}


void
ElectricalContact::do_compute(void)
{
  if (get_type(0) != NEUMANN)
    coeff_g(0) = _barrier + get_inner_voltage();

  if (get_type(1) != NEUMANN)
    coeff_g(1) = get_inner_voltage();

  if (get_type(2) != NEUMANN)
    coeff_g(2) = get_inner_voltage();

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



