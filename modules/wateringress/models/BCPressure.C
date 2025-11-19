// $Id$

#include "BCPressure.h"
#include "WIUtils.h"
#include "tibercad/base/SimulationOptions.h"

#include "tibercad/module/TiberModule.h"

using namespace libMesh;


void
BCPressure::do_init(void)
{
  _relative_humidity = get_option("relative_humidity", _relative_humidity);
}


void
BCPressure::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  double temp = SimulationOptions::temperature;
  double pressure = WIUtils::goff_gratch(temp);
  pressure *= _relative_humidity / 100;

  set_coefficients(1, 0, pressure);
}

