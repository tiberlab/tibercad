// $Id$

#include "ConstantChargeDensity.h"

#include "TiberModule.h"


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
