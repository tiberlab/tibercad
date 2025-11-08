// $Id: ConstantChargeDensity.C 3864 2014-06-13 10:24:26Z maufder $

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
