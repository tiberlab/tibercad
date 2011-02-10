
// $Id$

#include "ThermalBoundaryResistance.h"
#include "SimulationOptions.h"

TIBER_MODULE(ThermalBoundaryResistance,boundary,surface_resistance)
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

