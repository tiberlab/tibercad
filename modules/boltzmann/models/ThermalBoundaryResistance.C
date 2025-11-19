
// $Id: ThermalBoundaryResistance.C 2362 2011-02-21 01:02:31Z gromano $

#include "ThermalBoundaryResistance.h"
#include "tibercad/base/SimulationOptions.h"

#include "tibercad/module/TiberModule.h"

void
ThermalBoundaryResistance::do_init(void)
{

  _resistance = SimulationOptions::temperature;
  get_parameter("r_surf", _resistance);
  get_parameter("temperature", _temperature);

}

void
ThermalBoundaryResistance::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  set_coefficients(1.0/_resistance, 1.0, _temperature/_resistance);
}

