 // $Id$

#include "HeatReservoir.h"

TIBER_MODULE(HeatReservoir,Boundary,heat_reservoir)

void
HeatReservoir::do_init(void)
{
  get_parameter("temperature", _temperature);

}

void
HeatReservoir::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  set_coefficients(1.0, 0.0, _temperature);
}

