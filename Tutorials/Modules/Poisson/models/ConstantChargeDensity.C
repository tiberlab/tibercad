// $Id$

#include "ConstantChargeDensity.h"



// The first string is the class name, the second one
// is the type of the model (here it is a charge density model),
// the third one is the specific model implementation.
// The library name will then be charge_density_constant.so
TIBER_MODULE(ConstantChargeDensity, charge_density, constant)


void
ConstantChargeDensity::do_init(void)
{
  get_parameter("charge_density", _charge_density);
}


double
ConstantChargeDensity::calculate_charge_density(const Elem* elem, const Point& point)
{
  return _charge_density;
}
